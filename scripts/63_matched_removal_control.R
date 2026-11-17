# =====================================================================
# 63_matched_removal_control.R
# Stage   : 13_sensitivity
# Folder  : 13_sensitivity/scripts/   (NEW; writes R/ctl_tasks.R)
# Purpose : EXPLORATORY matched control for script 60: is the loss of the M1 deficiency when the 60
#           GIS_MAP sites are removed (n = 137, p = 0.127 at nsim 999) a sample-size effect?
#           Group G: the fixed GIS_MAP-excluded pattern (137 sites), 80 independent simulation seeds.
#           Group R: 80 random removals of 60 sites (n = 137 each), one test per subset.
#           Full M1 test (E3, joint ERL, nstep 2) at nsim 199 for every replicate (same generic code
#           as scripts 49-59), 4 workers, batches of 20. Retention = share with p < 0.05, compared
#           with the script-50 baseline (197 sites, nsim 199: 0.92).
# Outputs : outputs/63_report.txt, 63_batch_NN.csv, 63_matched_control_replicates.csv,
#           63_matched_control_summary.csv, 63_sessionInfo.txt; R/ctl_tasks.R
# Notes   : RUN ONCE (resumable by batch). Reading rules (set before running): GIS-SPECIFIC if
#           retention(R) >= 0.80 and retention(G) <= 0.60 and Fisher exact p (G vs R) < 0.01 (the loss is
#           not explained by sample size); SAMPLE-SIZE if retention(G) is within 0.15 of retention(R) and
#           Fisher p >= 0.05; otherwise INTERMEDIATE. nsim 199 p-scale differs from nsim 999 (scripts 51-52):
#           only same-nsim comparisons are used. Restart R before running.
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

OUT_DIR <- dir_of("13_sensitivity", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
f_res <- file.path(OUT_DIR, "63_matched_control_summary.csv")
if (file.exists(f_res)) stop("63_matched_control_summary.csv exists: this analysis runs once. Stop; do not rerun.")
GENERIC_FILE <- file.path(REPO_DIR, "R", "robust_generic.R"); JIT_FILE <- file.path(REPO_DIR, "R", "jitter_tasks.R")
stopifnot(file.exists(GENERIC_FILE), file.exists(JIT_FILE))
CTL_FILE <- file.path(REPO_DIR, "R", "ctl_tasks.R")
N_G <- 80; N_R <- 80; NSIM <- TESTS$nsim_recovery; BATCH <- 20; N_WORK <- 4
SEED_C <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 170000

# ---- open inputs BEFORE sink ----
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
SITES_G <- read.csv(f_sites, stringsAsFactors = FALSE)
stopifnot(nrow(SITES_G) == 197, max(abs(X_obs_km$x - SITES_G$x_m / 1000), abs(X_obs_km$y - SITES_G$y_m / 1000)) < 1e-3)
f_s50 <- file.path(REPO_DIR, "11_generative_gof", "outputs", "50_jitter_summary.csv"); stopifnot(file.exists(f_s50))
f_s60 <- file.path(OUT_DIR, "60_gis_excluded_result.csv"); stopifnot(file.exists(f_s60))
S50 <- read.csv(f_s50, stringsAsFactors = FALSE); S60 <- read.csv(f_s60, stringsAsFactors = FALSE)
source(GENERIC_FILE); source(JIT_FILE)

# ---- control task (dumped for workers) ----
ctl_task <- function(job) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch({
    X <- X_obs_km[job$keep]; r <- jit_test(X, job$nsim, job$seed)
    list(n = npoints(X), p = r$p, usable = r$usable, po2 = r$pobs[1], po3 = r$pobs[2], po5 = r$pobs[3],
         ps2 = r$psim[1], ps3 = r$psim[2], ps5 = r$psim[3], err = NA_character_)
  }, error = function(e) list(err = conditionMessage(e)))
  c(list(grp = job$grp, i = job$i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
crow <- function(L) do.call(rbind, lapply(L, function(x) {
  g <- function(f) if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  data.frame(grp = x$grp, i = x$i, sec = x$sec, n = g("n"), p = g("p"), usable = g("usable"), po2 = g("po2"), po3 = g("po3"), po5 = g("po5"),
             ps2 = g("ps2"), ps3 = g("ps3"), ps5 = g("ps5"), err = if (is.null(x$err)) NA_character_ else as.character(x$err), stringsAsFactors = FALSE)
}))
init_worker7 <- function(repo, cache, generic, jit, ctl, sites_path) {
  assign("REPO_DIR", repo, envir = globalenv()); Sys.setenv(RECOVERY_CACHE = cache)
  source(file.path(repo, "R", "recovery_worker.R"), local = globalenv())
  source(file.path(repo, "R", "recovery_functions.R"), local = globalenv())
  assign("SITES_G", read.csv(sites_path, stringsAsFactors = FALSE), envir = globalenv())
  source(generic, local = globalenv()); source(jit, local = globalenv()); source(ctl, local = globalenv())
  TRUE
}
dump(c("ctl_task", "crow"), file = CTL_FILE)
wilson <- function(k, n) as.numeric(suppressWarnings(prop.test(k, n, correct = FALSE)$conf.int))

REPORT <- file.path(OUT_DIR, "63_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("63_matched_removal_control (EXPLORATORY) | ", format(Sys.time()), "\n", sep = "")
  cat(sprintf("group G: GIS_MAP-excluded subset x %d seeds | group R: %d random removals of 60 sites | nsim %d | batch %d | workers %d | seed base %d\n", N_G, N_R, NSIM, BATCH, N_WORK, SEED_C))
  cat(sprintf("reference: script 60 (GIS-excluded, nsim 999) p %.4f | script 50 baseline retention (197 sites, nsim 199) %.3f\n", S60$p, S50$retention[S50$type == "base"]))
  gis <- SITES_G$coordinate_source == "GIS_MAP"; keepG <- !gis
  set.seed(SEED_C); KEEPR <- lapply(seq_len(N_R), function(i) { k <- rep(TRUE, 197); k[sample.int(197, 60)] <- FALSE; k })
  cat(sprintf("subset sizes: G %d | R %d-%d | random removals containing k GIS_MAP sites: mean %.1f (expected 18.3)\n", sum(keepG),
              min(vapply(KEEPR, sum, 1L)), max(vapply(KEEPR, sum, 1L)), mean(vapply(KEEPR, function(k) sum(!k & gis), 1L))))
  
  cl <- makeCluster(N_WORK)
  ini <- clusterCall(cl, init_worker7, REPO_DIR, CACHE_FILE, GENERIC_FILE, JIT_FILE, CTL_FILE, f_sites)
  cat("workers initialised:", all(unlist(ini)), "\n")
  ex <- clusterCall(cl, function() c(exists("ctl_task"), exists("jit_test"), exists("SITES_G"), "terra" %in% loadedNamespaces()))
  cat("worker check (ctl_task, jit_test, SITES_G present; terra loaded):", paste(ex[[1]], collapse = " "), "\n")
  
  cat("\n== Worker-vs-main identity check (one random-removal job, nsim 29) ==\n")
  jt <- list(grp = "R", i = 901, keep = KEEPR[[1]], nsim = 29, seed = SEED_C + 901)
  m <- crow(list(ctl_task(jt))); w <- crow(parLapply(cl[1], list(jt), ctl_task))
  cat(sprintf("main: n %d p %.4f pairs %.0f/%.0f/%.0f sims %.4f/%.4f/%.4f | worker: n %d p %.4f pairs %.0f/%.0f/%.0f sims %.4f/%.4f/%.4f | worker time %.0f s\n",
              m$n, m$p, m$po2, m$po3, m$po5, m$ps2, m$ps3, m$ps5, w$n, w$p, w$po2, w$po3, w$po5, w$ps2, w$ps3, w$ps5, w$sec))
  if (!(isTRUE(all.equal(m$p, w$p)) && isTRUE(all.equal(c(m$ps2, m$ps3, m$ps5), c(w$ps2, w$ps3, w$ps5))) && is.na(m$err) && is.na(w$err))) stop("worker path differs from main path; stop and paste output")
  
  cat("\n== Batches ==\n")
  JOBS <- list()
  for (i in seq_len(max(N_G, N_R))) {
    if (i <= N_G) JOBS[[length(JOBS) + 1]] <- list(grp = "G", i = i, keep = keepG, nsim = NSIM, seed = SEED_C + 1000 + i)
    if (i <= N_R) JOBS[[length(JOBS) + 1]] <- list(grp = "R", i = i, keep = KEEPR[[i]], nsim = NSIM, seed = SEED_C + 5000 + i)
  }
  nb <- ceiling(length(JOBS) / BATCH)
  for (b in seq_len(nb)) {
    f_b <- file.path(OUT_DIR, sprintf("63_batch_%02d.csv", b))
    if (file.exists(f_b)) { cat(sprintf("  batch %d: found, skipped\n", b)); next }
    grp <- JOBS[((b - 1) * BATCH + 1):min(b * BATCH, length(JOBS))]; t0 <- unname(proc.time()["elapsed"])
    db <- crow(parLapplyLB(cl, grp, ctl_task)); db$batch <- b; write.csv(db, f_b, row.names = FALSE)
    cat(sprintf("  batch %d of %d: %d replicates, %.0f s wall, failed %d [%s]\n", b, nb, nrow(db), unname(proc.time()["elapsed"]) - t0, sum(!is.na(db$err)), format(Sys.time(), "%H:%M")))
  }
  
  cat("\n== Retention ==\n")
  fs <- file.path(OUT_DIR, sprintf("63_batch_%02d.csv", seq_len(nb)))
  ALL <- do.call(rbind, lapply(fs[file.exists(fs)], read.csv, stringsAsFactors = FALSE)); write.csv(ALL, file.path(OUT_DIR, "63_matched_control_replicates.csv"), row.names = FALSE)
  SM <- do.call(rbind, lapply(c("G", "R"), function(g) { d <- ALL[ALL$grp == g & !is.na(ALL$p), ]; k <- sum(d$p < TESTS$alpha); ci <- wilson(k, nrow(d)); q <- quantile(d$p, c(.1, .25, .5, .75, .9))
  data.frame(group = g, description = if (g == "G") "GIS_MAP-excluded (fixed subset), seeds vary" else "60 random sites removed, subset varies", reps = nrow(d), failed = sum(ALL$grp == g) - nrow(d),
             n_median = median(d$n), retained = k, retention = k / nrow(d), ci_lo = ci[1], ci_hi = ci[2], p_q10 = q[1], p_q25 = q[2], p_median = q[3], p_q75 = q[4], p_q90 = q[5],
             pair_ratio_2 = median(d$po2 / d$ps2), pair_ratio_3 = median(d$po3 / d$ps3), pair_ratio_5 = median(d$po5 / d$ps5), min_usable = min(d$usable)) }))
  print(SM, digits = 3, row.names = FALSE)
  g <- SM[SM$group == "G", ]; r <- SM[SM$group == "R", ]
  ft <- fisher.test(matrix(c(g$retained, g$reps - g$retained, r$retained, r$reps - r$retained), 2))
  cat(sprintf("\nretention G %.3f (%.3f-%.3f) | R %.3f (%.3f-%.3f) | baseline (197 sites) %.3f | Fisher exact p (G vs R) %.4f\n", g$retention, g$ci_lo, g$ci_hi, r$retention, r$ci_lo, r$ci_hi, S50$retention[S50$type == "base"], ft$p.value))
  verdict <- if (r$retention >= 0.80 && g$retention <= 0.60 && ft$p.value < 0.01) "GIS-SPECIFIC: the loss is not explained by sample size" else
    if (abs(g$retention - r$retention) <= 0.15 && ft$p.value >= 0.05) "SAMPLE-SIZE: random removal of 60 sites loses retention similarly" else "INTERMEDIATE: partly sample size, partly the removed sites"
  cat("verdict (pre-set rule):", verdict, "\n")
  SM$fisher_p <- ft$p.value; SM$verdict <- verdict; SM$baseline_retention <- S50$retention[S50$type == "base"]; SM$date <- format(Sys.time()); write.csv(SM, f_res, row.names = FALSE)
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink(); if (exists("cl")) try(stopCluster(cl), silent = TRUE) })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "63_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")