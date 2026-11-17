# =====================================================================
# 39_forced_hierarchy_smoke.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file; REWRITES R/recovery_cells.R
#           after backing it up to outputs/recovery_cells_before_39.R)
# Purpose : exercise M3 and M4 code paths that script 38 could not reach
#           (nsim 19 -> minimum ERL p 0.05, rejection p < 0.05 impossible).
#           run_hierarchy gets a 'force' argument (all three stages run) and
#           stage-labelled failures; 19 cells x 1 rep, nsim 29 (minimum p 0.033),
#           forced path, 4 workers.
# Outputs : outputs/39_report.txt, 39_smoke_forced.csv, 39_sessionInfo.txt,
#           recovery_cells_before_39.R; R/recovery_cells.R (rewritten)
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           Timings feed the Stage 1 projection. p-values and classes of single
#           replicates at nsim 29 carry no inferential information.
#           Restart R before running. Workers load the cache only (no terra).
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
N_WORK     <- 4
NSIM_SMOKE <- 29
I_SMOKE    <- 901

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULTS <- setNames(gm$mult, gm$cell)
OBS <- readRDS(OBS_FILE)
source(TESTS_FILE)
source(CELLS_FILE)          # defines the script-38 versions in this session
file.copy(CELLS_FILE, file.path(OUT_DIR, "recovery_cells_before_39.R"), overwrite = TRUE)

# ---- replacement functions (force argument, stage-labelled failures) ----
run_hierarchy <- function(X, nsim, bw, rGv, seed, force = FALSE) {
  R <- list(p1 = NA_real_, p3 = NA_real_, p4 = NA_real_, zmax3 = NA_real_, zmax4 = NA_real_,
            usable1 = NA_real_, usable3 = NA_real_, usable4 = NA_real_, sigma2 = NA_real_, alpha = NA_real_,
            gamma = NA_real_, r_sel = NA_real_, sat_sel = NA_real_, sec1 = NA_real_, sec3 = NA_real_, sec4 = NA_real_,
            class = "failed_M1")
  t1 <- unname(proc.time()["elapsed"])
  a <- tryCatch(m1_test_r(X, nsim, bw, rGv, seed + 1), error = function(e) list(p = NA_real_, usable = 0))
  R$sec1 <- unname(proc.time()["elapsed"]) - t1; R$p1 <- a$p; R$usable1 <- a$usable
  if (is.na(a$p)) return(R)
  if (!force && a$p >= TESTS$alpha) { R$class <- "Poisson"; return(R) }
  R$class <- "failed_M3"
  t3 <- unname(proc.time()["elapsed"])
  b <- tryCatch(m3_test_r(X, nsim, bw, rGv, seed + 2),
                error = function(e) list(p = NA_real_, usable = 0, zmax = NA_real_, sigma2 = NA_real_, alpha = NA_real_))
  R$sec3 <- unname(proc.time()["elapsed"]) - t3; R$p3 <- b$p; R$usable3 <- b$usable; R$zmax3 <- b$zmax
  R$sigma2 <- b$sigma2; R$alpha <- b$alpha
  if (is.na(b$p)) return(R)
  if (!force && b$p >= TESTS$alpha) { R$class <- "LGCP"; return(R) }
  R$class <- "failed_M4"
  t4 <- unname(proc.time()["elapsed"])
  d <- tryCatch(m4_search_test(X, nsim, bw, rGv, seed + 3),
                error = function(e) list(p = NA_real_, usable = 0, zmax = NA_real_, gamma = NA_real_, r_sel = NA_real_, sat_sel = NA_real_))
  R$sec4 <- unname(proc.time()["elapsed"]) - t4; R$p4 <- d$p; R$usable4 <- d$usable; R$zmax4 <- d$zmax
  R$gamma <- d$gamma; R$r_sel <- d$r_sel; R$sat_sel <- d$sat_sel
  if (is.na(d$p) || is.na(d$zmax) || is.na(b$zmax)) return(R)
  R$class <- if (R$p1 >= TESTS$alpha) "Poisson" else if (b$p >= TESTS$alpha) "LGCP" else
    if (d$p >= TESTS$alpha && d$zmax <= 0.5 * b$zmax) "Geyer" else "inadequate"
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
dump(c("REC_NUMF", "env_eta", "b0_for", "jitter_pat", "gen_cell", "rG_for", "fast_curves_r", "curves_test_r", "m1_test_r",
       "m3_test_r", "m4_search_test", "run_hierarchy", "rec_task", "rec_rows_df", "init_worker2"), file = CELLS_FILE)

WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }

REPORT <- file.path(OUT_DIR, "39_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("39_forced_hierarchy_smoke | ", format(Sys.time()), "\n", sep = "")
    cat(sprintf("cells %d x 1 rep (i = %d) | nsim %d (min ERL p %.3f) | forced M1 -> M3 -> M4 | workers %d\n",
                nrow(RECOVERY_CELLS), I_SMOKE, NSIM_SMOKE, 1 / (NSIM_SMOKE + 1), N_WORK))
    cat("recovery_cells.R rewritten; backup outputs/recovery_cells_before_39.R\n")
    cl <- makeCluster(N_WORK)
    ini <- clusterCall(cl, init_worker2, REPO_DIR, CACHE_FILE, TESTS_FILE, CELLS_FILE, OBS_FILE, MULTS)
    cat("workers initialised:", all(unlist(ini)), "\n")
    jobs <- lapply(RECOVERY_CELLS$cell, function(cc) list(cell = cc, i = I_SMOKE, nsim = NSIM_SMOKE, force = TRUE))
    t0 <- unname(proc.time()["elapsed"])
    res <- parLapplyLB(cl, jobs, rec_task)
    wall <- unname(proc.time()["elapsed"]) - t0
    SM <- rec_rows_df(res)
    write.csv(SM, file.path(OUT_DIR, "39_smoke_forced.csv"), row.names = FALSE)
    cat(sprintf("wall %.0f s for %d forced replicates\n\n", wall, nrow(SM)))
    print(SM[, c("cell", "class", "n", "p1", "p3", "p4", "zmax3", "zmax4", "sigma2", "alpha", "gamma",
                 "r_sel", "sat_sel", "usable1", "usable3", "usable4", "sec1", "sec3", "sec4", "err")],
          digits = 3, row.names = FALSE)
    cat("\n-- Checks --\n")
    cat("failed (any stage):", sum(grepl("^failed", SM$class)), "of", nrow(SM), "\n")
    cat("all three p-values numeric:", sum(!is.na(SM$p1) & !is.na(SM$p3) & !is.na(SM$p4)), "of", nrow(SM), "\n")
    cat("both z-max finite:", sum(is.finite(SM$zmax3) & is.finite(SM$zmax4)), "of", nrow(SM), "\n")
    cat("distinct selected (r, sat) pairs:", nrow(unique(SM[!is.na(SM$r_sel), c("r_sel", "sat_sel")])), "\n")
    print(table(r_sel = SM$r_sel, sat_sel = SM$sat_sel, useNA = "ifany"))
    cat(sprintf("seconds per test at nsim %d in workers: M1 mean %.1f | M3 mean %.1f | M4 incl. profile fits mean %.1f\n",
                NSIM_SMOKE, mean(SM$sec1, na.rm = TRUE), mean(SM$sec3, na.rm = TRUE), mean(SM$sec4, na.rm = TRUE)))
    cat("classes (natural rule, single forced replicates; NOT inference):\n"); print(table(SM$class))
    cat("\n== Warning tally (main session only; worker warnings not captured) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink(); if (exists("cl")) try(stopCluster(cl), silent = TRUE) })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "39_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")