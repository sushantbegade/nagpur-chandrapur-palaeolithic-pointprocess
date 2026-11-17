# =====================================================================
# 40_stage1_recovery_runner.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file; REWRITES R/recovery_cells.R
#           after backing it up to outputs/recovery_cells_before_40.R)
# Purpose : Stage 1 of the recovery study: 19 cells x 100 reps, nsim 99, 4 workers,
#           checkpointed batches (4 reps x 19 cells), hierarchy M1 -> M3 -> M4.
#           Adds a floored-sd z-max (primary) with raw z-max, argmax radius and curve
#           stored; stage-labelled failures. BLIND monitoring: the report prints
#           only timing and failure counts, never classes or p-values.
# Outputs : outputs/40_report.txt, 40_s1_batch_NN.csv, 40_stage1_all.csv,
#           40_sessionInfo.txt, protocol_deviation_log_before_40.csv,
#           recovery_cells_before_40.R; R/recovery_cells.R (rewritten);
#           docs/protocol_deviation_log.csv (2 rows appended)
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           Resumable: finished batches are skipped. Restart R before running.
#           Keep the machine awake. Workers load the cache only (no terra).
#           Stage 2 (reps 101-200) is decided in script 41 by the pre-set rule.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))
suppressPackageStartupMessages(library(parallel))

OUT_DIR    <- dir_of("10_simulation_recovery", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
TESTS_FILE <- file.path(REPO_DIR, "R", "recovery_tests.R"); stopifnot(file.exists(TESTS_FILE))
CELLS_FILE <- file.path(REPO_DIR, "R", "recovery_cells.R"); stopifnot(file.exists(CELLS_FILE))
OBS_FILE   <- file.path(REPO_DIR, "data_derived", "recovery_obs_cache.rds"); stopifnot(file.exists(OBS_FILE))
N_WORK <- 4; N_REP <- 100; BATCH <- 4; NSIM <- 99

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULTS <- setNames(gm$mult, gm$cell)
OBS <- readRDS(OBS_FILE)
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")
source(TESTS_FILE)
source(CELLS_FILE)
file.copy(CELLS_FILE, file.path(OUT_DIR, "recovery_cells_before_40.R"), overwrite = TRUE)

# ---- replacement functions ----
REC_NUMF <- c("sec", "n_raw", "n", "n_dropped", "p1", "p3", "p4", "zmax3", "zmax4", "zraw3", "zraw4",
              "zarg_r3", "zarg_r4", "zarg_cv3", "zarg_cv4", "usable1", "usable3", "usable4",
              "sigma2", "alpha", "gamma", "r_sel", "sat_sel", "sec1", "sec3", "sec4")
num1 <- function(x) if (is.null(x)) NA_real_ else as.numeric(x)
zsum <- function(obs, S, floorf = 0.25) {
  m <- rowMeans(S); s <- apply(S, 1, sd)
  med <- median(s[is.finite(s) & s > 0]); if (!is.finite(med)) med <- 1
  sf <- pmax(s, floorf * med)
  list(raw = (obs - m) / s, flo = (obs - m) / sf)
}
curves_test_r <- function(Xobs, sims, bw, rGv) {
  f0 <- fast_fit(Xobs); ob <- fast_curves_r(Xobs, f0, bw, rGv)
  sL <- matrix(NA_real_, length(rL), length(sims)); sG <- matrix(NA_real_, length(rGv), length(sims))
  dsim <- integer(length(sims))
  for (k in seq_along(sims)) {
    Xr <- to_poly(sims[[k]]); Xs <- clean_pat(Xr); dsim[k] <- npoints(Xr) - npoints(Xs)
    fs <- fast_fit(Xs); cs <- fast_curves_r(Xs, fs, bw, rGv); sL[, k] <- cs$L; sG[, k] <- cs$g
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rGv, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  zmax <- NA_real_; zraw <- NA_real_; zr <- NA_real_; zcv <- NA_real_
  if (sum(okc) >= 2) {
    zL <- zsum(ob$L, sL[, okc, drop = FALSE]); zG <- zsum(ob$g, sG[, okc, drop = FALSE])
    fl <- abs(c(zL$flo, zG$flo)); rw <- abs(c(zL$raw, zG$raw)); rr <- c(rL, rGv)
    fl[!is.finite(fl)] <- NA; rw[!is.finite(rw)] <- NA
    if (any(!is.na(fl))) zmax <- max(fl, na.rm = TRUE)
    if (any(!is.na(rw))) { ia <- which.max(rw); zraw <- rw[ia]; zr <- rr[ia]; zcv <- if (ia <= length(rL)) 1 else 2 }
  }
  list(p = erl_p(csL, csG, 2), usable = sum(okc), zmax = zmax, zraw = zraw, zarg_r = zr, zarg_cv = zcv, drop_sim = sum(dsim))
}
run_hierarchy <- function(X, nsim, bw, rGv, seed, force = FALSE) {
  R <- list(p1 = NA_real_, p3 = NA_real_, p4 = NA_real_, zmax3 = NA_real_, zmax4 = NA_real_,
            zraw3 = NA_real_, zraw4 = NA_real_, zarg_r3 = NA_real_, zarg_r4 = NA_real_, zarg_cv3 = NA_real_, zarg_cv4 = NA_real_,
            usable1 = NA_real_, usable3 = NA_real_, usable4 = NA_real_, sigma2 = NA_real_, alpha = NA_real_,
            gamma = NA_real_, r_sel = NA_real_, sat_sel = NA_real_, sec1 = NA_real_, sec3 = NA_real_, sec4 = NA_real_,
            class = "failed_M1")
  t1 <- unname(proc.time()["elapsed"])
  a <- tryCatch(m1_test_r(X, nsim, bw, rGv, seed + 1), error = function(e) list(p = NA_real_, usable = 0))
  R$sec1 <- unname(proc.time()["elapsed"]) - t1; R$p1 <- num1(a$p); R$usable1 <- num1(a$usable)
  if (is.na(R$p1)) return(R)
  if (!force && R$p1 >= TESTS$alpha) { R$class <- "Poisson"; return(R) }
  R$class <- "failed_M3"
  t3 <- unname(proc.time()["elapsed"])
  b <- tryCatch(m3_test_r(X, nsim, bw, rGv, seed + 2), error = function(e) list(p = NA_real_, usable = 0))
  R$sec3 <- unname(proc.time()["elapsed"]) - t3; R$p3 <- num1(b$p); R$usable3 <- num1(b$usable)
  R$zmax3 <- num1(b$zmax); R$zraw3 <- num1(b$zraw); R$zarg_r3 <- num1(b$zarg_r); R$zarg_cv3 <- num1(b$zarg_cv)
  R$sigma2 <- num1(b$sigma2); R$alpha <- num1(b$alpha)
  if (is.na(R$p3)) return(R)
  if (!force && R$p3 >= TESTS$alpha) { R$class <- "LGCP"; return(R) }
  R$class <- "failed_M4"
  t4 <- unname(proc.time()["elapsed"])
  d <- tryCatch(m4_search_test(X, nsim, bw, rGv, seed + 3), error = function(e) list(p = NA_real_, usable = 0))
  R$sec4 <- unname(proc.time()["elapsed"]) - t4; R$p4 <- num1(d$p); R$usable4 <- num1(d$usable)
  R$zmax4 <- num1(d$zmax); R$zraw4 <- num1(d$zraw); R$zarg_r4 <- num1(d$zarg_r); R$zarg_cv4 <- num1(d$zarg_cv)
  R$gamma <- num1(d$gamma); R$r_sel <- num1(d$r_sel); R$sat_sel <- num1(d$sat_sel)
  if (is.na(R$p4) || is.na(R$zmax4) || is.na(R$zmax3)) return(R)
  R$class <- if (R$p1 >= TESTS$alpha) "Poisson" else if (R$p3 >= TESTS$alpha) "LGCP" else
    if (R$p4 >= TESTS$alpha && R$zmax4 <= 0.5 * R$zmax3) "Geyer" else "inadequate"
  R
}
rec_task <- function(job) {
  t0 <- unname(proc.time()["elapsed"])
  cc <- job$cell; i <- job$i; ci <- match(cc, RECOVERY_CELLS$cell)
  sb <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 100000
  n_t <- RECOVERY_CELLS$n[ci]; bw <- bw_for_n(n_t); rGv <- rG_for(n_t)
  out <- tryCatch({
    X0 <- gen_cell(cc, sb + 1000 * ci + i); X <- clean_pat(X0)
    res <- run_hierarchy(X, job$nsim, bw, rGv, sb + 500000 + 1000 * ci + i, force = isTRUE(job$force))
    c(list(n_raw = npoints(X0), n = npoints(X), n_dropped = npoints(X0) - npoints(X), err = NA_character_), res)
  }, error = function(e) list(class = "failed_gen", err = conditionMessage(e)))
  c(list(cell = cc, i = i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
dump(c("REC_NUMF", "num1", "zsum", "env_eta", "b0_for", "jitter_pat", "gen_cell", "rG_for", "fast_curves_r", "curves_test_r",
       "m1_test_r", "m3_test_r", "m4_search_test", "run_hierarchy", "rec_task", "rec_rows_df", "init_worker2"), file = CELLS_FILE)

# ---- deviation log (before sink; 2 clarification rows) ----
D <- "2026-11-11"; AB <- "author-delegated (assistant decision; clarification, no frozen item changed)"
ROWS <- data.frame(
  date = D,
  item = c("Clarification: max standardised deviation uses an sd floor (0.25 x median simulated sd per curve)",
           "Clarification: Stage 1 recovery run - blind monitoring, batches, failure handling"),
  reason = c("Script 39 forced smoke (single replicates, nsim 29) showed pointwise |z| up to 15 alongside non-rejecting ERL p (G3: p4 0.133, zmax4 11.2; N1: p4 0.067, zmax4 9.3), suggesting that z-max can be driven by a single radius with tiny simulated sd (INFERENCE). The floor was chosen from synthetic smoke output only, before any Stage 1 result. Unfloored z-max, argmax radius and curve are stored for every replicate.",
             "Stage 1 = 19 cells x 100 reps, nsim 99, 25 batches of 4 reps x 19 cells. The report prints only timing and failure counts (no classes or p-values) so that the staged-extension rule is not subject to peeking. Failed replicates (stage-labelled) are excluded from denominators and reported; a cell with more than 10% failures is flagged non-evaluable."),
  effect_on_inference = c("Geyer classification depends on the floored z-max; the unfloored version is reported as sensitivity. No archaeological model has been fitted.",
                          "Interim results are not inspected; Gate 6 analysis (script 41) is run once after Stage 1. No archaeological model has been fitted."),
  approved_by = AB, stringsAsFactors = FALSE)
if (identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_40.csv"), overwrite = TRUE)
  NEW <- ROWS[!(ROWS$item %in% LOG0$item), , drop = FALSE]
  if (nrow(NEW)) write.csv(rbind(LOG0, NEW), f_log, row.names = FALSE)
  LOG_MSG <- sprintf("Appended %d row(s) to docs/protocol_deviation_log.csv (backup in outputs).", nrow(NEW))
} else {
  write.csv(ROWS, file.path(OUT_DIR, "40_log_rows_for_manual_append.csv"), row.names = FALSE)
  LOG_MSG <- "Log columns differ; rows saved to 40_log_rows_for_manual_append.csv (append manually)."
}

WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }

REPORT <- file.path(OUT_DIR, "40_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("40_stage1_recovery_runner | ", format(Sys.time()), "\n", sep = "")
    cat(sprintf("cells %d x %d reps | nsim %d | batch %d reps (%d jobs) | workers %d\n",
                nrow(RECOVERY_CELLS), N_REP, NSIM, BATCH, BATCH * nrow(RECOVERY_CELLS), N_WORK))
    cat(LOG_MSG, "\nrecovery_cells.R rewritten; backup outputs/recovery_cells_before_40.R\n")
    cl <- makeCluster(N_WORK)
    ini <- clusterCall(cl, init_worker2, REPO_DIR, CACHE_FILE, TESTS_FILE, CELLS_FILE, OBS_FILE, MULTS)
    cat("workers initialised:", all(unlist(ini)), "\n")
    
    cat("\n== Pre-flight: one forced G3 replicate (i = 902, nsim 29) ==\n")
    pf <- rec_rows_df(parLapply(cl[1], list(list(cell = "G3", i = 902, nsim = 29, force = TRUE)), rec_task))
    print(pf[, c("cell", "class", "zmax3", "zraw3", "zarg_r3", "zarg_cv3", "zmax4", "zraw4", "zarg_r4", "zarg_cv4", "err")],
          digits = 4, row.names = FALSE)
    need <- c("zmax3", "zraw3", "zarg_r3", "zarg_cv3", "zmax4", "zraw4", "zarg_r4", "zarg_cv4", "p1", "p3", "p4")
    if (any(is.na(pf[, need]))) stop("pre-flight: new fields contain NA; stop and paste output")
    cat("pre-flight OK\n")
    
    cat("\n== Stage 1 batches (blind: timing and failures only) ==\n")
    nb <- ceiling(N_REP / BATCH)
    for (b in seq_len(nb)) {
      f_b <- file.path(OUT_DIR, sprintf("40_s1_batch_%02d.csv", b))
      if (file.exists(f_b)) { cat(sprintf("  batch %2d: found, skipped\n", b)); next }
      ii <- ((b - 1) * BATCH + 1):min(b * BATCH, N_REP)
      jobs <- unlist(lapply(ii, function(i) lapply(RECOVERY_CELLS$cell, function(cc)
        list(cell = cc, i = i, nsim = NSIM, force = FALSE))), recursive = FALSE)
      t0 <- unname(proc.time()["elapsed"])
      res <- parLapplyLB(cl, jobs, rec_task)
      db <- rec_rows_df(res); db$batch <- b; db$batch_wall_sec <- unname(proc.time()["elapsed"]) - t0
      write.csv(db, f_b, row.names = FALSE)
      cat(sprintf("  batch %2d (reps %d-%d): %d replicates, %.0f s wall, mean %.0f s per replicate, failed %d  [%s]\n",
                  b, min(ii), max(ii), nrow(db), unname(proc.time()["elapsed"]) - t0, mean(db$sec),
                  sum(grepl("^failed", db$class)), format(Sys.time(), "%m-%d %H:%M")))
    }
    
    fs <- file.path(OUT_DIR, sprintf("40_s1_batch_%02d.csv", seq_len(nb)))
    ALL <- do.call(rbind, lapply(fs[file.exists(fs)], read.csv, stringsAsFactors = FALSE))
    write.csv(ALL, file.path(OUT_DIR, "40_stage1_all.csv"), row.names = FALSE)
    cat("\n== Stage 1 completeness (no classes printed) ==\n")
    TB <- do.call(rbind, lapply(split(ALL, ALL$cell), function(d) data.frame(
      cell = d$cell[1], replicates = nrow(d), failed = sum(grepl("^failed", d$class)),
      mean_n = mean(d$n, na.rm = TRUE), mean_sec = mean(d$sec))))
    print(TB[match(RECOVERY_CELLS$cell, TB$cell), ], digits = 4, row.names = FALSE)
    cat(sprintf("\ntotal replicates %d of %d | failed %d | wall mean per batch %.0f s\n", nrow(ALL),
                nrow(RECOVERY_CELLS) * N_REP, sum(grepl("^failed", ALL$class)), mean(ALL$batch_wall_sec)))
    cat("\n== Warning tally (main session only; worker warnings not captured) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink(); if (exists("cl")) try(stopCluster(cl), silent = TRUE) })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "40_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")