# =====================================================================
# 59_jitter_sensitivity.R
# Stage   : 13_sensitivity
# Folder  : 13_sensitivity/scripts/   (NEW; writes R/jitter_sens_tasks.R)
# Purpose : coordinate-jitter SENSITIVITY scenarios for the frozen M1 deficiency (script 50 was the
#           primary scenario): FIELD fixed at 10 m; scenarios S1 (AUTHOR 10, GIS 100 m),
#           S2 (AUTHOR 100, GIS 500 m), S3 (AUTHOR 250, GIS 1000 m); 100 replicates each, full M1 test
#           per jittered pattern at nsim 199 (same test and code as script 50), 4 workers, batches of 20.
#           Draws outside the window or in NA-covariate pixels are redrawn (<= 50 tries; fallbacks counted).
# Outputs : outputs/59_report.txt, 59_batch_NN.csv, 59_jitter_sens_replicates.csv,
#           59_jitter_sens_summary.csv, 59_sessionInfo.txt; R/jitter_sens_tasks.R
# Notes   : RUN ONCE (resumable by batch). Retention = share of replicates with p < 0.05; frozen
#           criterion >= 0.80, read against the script-50 baseline retention (0.92 at nsim 199).
#           Scenario sigmas are sensitivity values from the freeze; they are NOT estimates of the
#           actual coordinate error (GIS-map datum unknown). Restart R before running.
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
f_res <- file.path(OUT_DIR, "59_jitter_sens_summary.csv")
if (file.exists(f_res)) stop("59_jitter_sens_summary.csv exists: this analysis runs once. Stop; do not rerun.")
GENERIC_FILE <- file.path(REPO_DIR, "R", "robust_generic.R"); JIT_FILE <- file.path(REPO_DIR, "R", "jitter_tasks.R")
stopifnot(file.exists(GENERIC_FILE), file.exists(JIT_FILE))
JIT2_FILE <- file.path(REPO_DIR, "R", "jitter_sens_tasks.R")
N_R <- 100; NSIM <- TESTS$nsim_recovery; BATCH <- 20; N_WORK <- 4
SEED_J2 <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 140000
SC <- list(S1 = c(FIELD = 10, AUTHOR_REPORTED = 10, GIS_MAP = 100),
           S2 = c(FIELD = 10, AUTHOR_REPORTED = 100, GIS_MAP = 500),
           S3 = c(FIELD = 10, AUTHOR_REPORTED = 250, GIS_MAP = 1000))

# ---- open inputs BEFORE sink ----
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
SITES_G <- read.csv(f_sites, stringsAsFactors = FALSE)
stopifnot(nrow(SITES_G) == 197, all(SITES_G$coordinate_source %in% names(SC$S1)))
f_s50 <- file.path(REPO_DIR, "11_generative_gof", "outputs", "50_jitter_summary.csv"); stopifnot(file.exists(f_s50))
S50 <- read.csv(f_s50, stringsAsFactors = FALSE)
source(GENERIC_FILE); source(JIT_FILE)

# ---- scenario jitter + task (dumped for workers) ----
jitter_obs2 <- function(seed, sig_class_m) {
  set.seed(seed)
  sig <- unname(sig_class_m[SITES_G$coordinate_source]) / 1000
  x <- X_obs_km$x; y <- X_obs_km$y; fb <- 0L; disp <- numeric(length(x))
  for (i in seq_along(x)) {
    ok <- FALSE
    for (tr in seq_len(50)) {
      xn <- x[i] + rnorm(1, 0, sig[i]); yn <- y[i] + rnorm(1, 0, sig[i])
      if (inside.owin(xn, yn, W_km)) {
        z <- vapply(covs, function(im) lookup.im(im, xn, yn, naok = TRUE), numeric(1))
        if (!anyNA(z)) { disp[i] <- sqrt((xn - x[i])^2 + (yn - y[i])^2); x[i] <- xn; y[i] <- yn; ok <- TRUE; break }
      }
    }
    if (!ok) fb <- fb + 1L
  }
  list(X = ppp(x, y, window = W_km, check = FALSE), fallbacks = fb, disp = disp)
}
jit_task2 <- function(job) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch({
    jj <- jitter_obs2(job$seed_j, job$sig); r <- jit_test(jj$X, job$nsim, job$seed_t)
    list(n = npoints(jj$X), p = r$p, usable = r$usable, n_sim = r$n_sim, po2 = r$pobs[1], po3 = r$pobs[2], po5 = r$pobs[3],
         ps2 = r$psim[1], ps3 = r$psim[2], ps5 = r$psim[3], fallbacks = jj$fallbacks, med_disp_km = median(jj$disp), err = NA_character_)
  }, error = function(e) list(err = conditionMessage(e)))
  c(list(scn = job$scn, i = job$i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
jrow2 <- function(L) do.call(rbind, lapply(L, function(x) {
  g <- function(f) if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  data.frame(scn = x$scn, i = x$i, sec = x$sec, n = g("n"), p = g("p"), usable = g("usable"), n_sim = g("n_sim"),
             po2 = g("po2"), po3 = g("po3"), po5 = g("po5"), ps2 = g("ps2"), ps3 = g("ps3"), ps5 = g("ps5"),
             fallbacks = g("fallbacks"), med_disp_km = g("med_disp_km"), err = if (is.null(x$err)) NA_character_ else as.character(x$err), stringsAsFactors = FALSE)
}))
init_worker6 <- function(repo, cache, generic, jit, jit2, sites_path) {
  assign("REPO_DIR", repo, envir = globalenv()); Sys.setenv(RECOVERY_CACHE = cache)
  source(file.path(repo, "R", "recovery_worker.R"), local = globalenv())
  source(file.path(repo, "R", "recovery_functions.R"), local = globalenv())
  assign("SITES_G", read.csv(sites_path, stringsAsFactors = FALSE), envir = globalenv())
  source(generic, local = globalenv()); source(jit, local = globalenv()); source(jit2, local = globalenv())
  TRUE
}
dump(c("jitter_obs2", "jit_task2", "jrow2"), file = JIT2_FILE)
wilson <- function(k, n) as.numeric(suppressWarnings(prop.test(k, n, correct = FALSE)$conf.int))

REPORT <- file.path(OUT_DIR, "59_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("59_jitter_sensitivity | ", format(Sys.time()), "\n", sep = "")
  cat(sprintf("scenarios %d x %d replicates | nsim %d | batch %d | workers %d | seed base %d\n", length(SC), N_R, NSIM, BATCH, N_WORK, SEED_J2))
  cat("classes:", paste(sprintf("%s %d", names(table(SITES_G$coordinate_source)), as.integer(table(SITES_G$coordinate_source))), collapse = ", "), "\n")
  
  cat("\n== A. Jitter geometry (100 main-session draws per scenario) ==\n")
  cls <- SITES_G$coordinate_source; GE <- data.frame()
  for (sn in names(SC)) {
    D <- lapply(seq_len(100), function(i) jitter_obs2(SEED_J2 + 77 * match(sn, names(SC)) + i, SC[[sn]])); disp <- do.call(rbind, lapply(D, function(d) d$disp))
    for (cl in names(SC[[sn]])) GE <- rbind(GE, data.frame(scenario = sn, class = cl, sites = sum(cls == cl), sigma_m = SC[[sn]][[cl]],
                                                           expected_median_m = 1.1774 * SC[[sn]][[cl]], observed_median_m = 1000 * median(disp[, cls == cl]), observed_p95_m = 1000 * quantile(disp[, cls == cl], 0.95),
                                                           fallback_sites_per_draw = mean(vapply(D, function(d) d$fallbacks, 1L))))
  }
  print(GE, digits = 4, row.names = FALSE)
  
  cl <- makeCluster(N_WORK)
  ini <- clusterCall(cl, init_worker6, REPO_DIR, CACHE_FILE, GENERIC_FILE, JIT_FILE, JIT2_FILE, f_sites)
  cat("\nworkers initialised:", all(unlist(ini)), "\n")
  ex <- clusterCall(cl, function() c(exists("jit_task2"), exists("jit_test"), exists("SITES_G"), "terra" %in% loadedNamespaces()))
  cat("worker check (jit_task2, jit_test, SITES_G present; terra loaded):", paste(ex[[1]], collapse = " "), "\n")
  
  cat("\n== B. Worker smoke (1 per scenario, nsim 19; completion check only) ==\n")
  sm <- parLapply(cl, lapply(names(SC), function(sn) list(scn = sn, i = 901, nsim = 19, sig = SC[[sn]], seed_j = SEED_J2 + 901, seed_t = SEED_J2 + 5901)), jit_task2)
  smd <- jrow2(sm); print(smd[, c("scn", "n", "usable", "fallbacks", "med_disp_km", "sec", "err")], digits = 3, row.names = FALSE)
  if (any(!is.na(smd$err)) || any(is.na(smd$p))) stop("worker smoke failed; paste output")
  
  cat("\n== C. Batches ==\n")
  JOBS <- list()
  for (i in seq_len(N_R)) for (sn in names(SC))
    JOBS[[length(JOBS) + 1]] <- list(scn = sn, i = i, nsim = NSIM, sig = SC[[sn]], seed_j = SEED_J2 + 1000 + 1000 * match(sn, names(SC)) + i,
                                     seed_t = SEED_J2 + 500000 + 1000 * match(sn, names(SC)) + i)
  nb <- ceiling(length(JOBS) / BATCH)
  for (b in seq_len(nb)) {
    f_b <- file.path(OUT_DIR, sprintf("59_batch_%02d.csv", b))
    if (file.exists(f_b)) { cat(sprintf("  batch %2d: found, skipped\n", b)); next }
    grp <- JOBS[((b - 1) * BATCH + 1):min(b * BATCH, length(JOBS))]; t0 <- unname(proc.time()["elapsed"])
    res <- parLapplyLB(cl, grp, jit_task2); db <- jrow2(res); db$batch <- b
    write.csv(db, f_b, row.names = FALSE)
    cat(sprintf("  batch %2d of %d: %d replicates, %.0f s wall, failed %d [%s]\n", b, nb, nrow(db), unname(proc.time()["elapsed"]) - t0, sum(!is.na(db$err)), format(Sys.time(), "%H:%M")))
  }
  
  cat("\n== D. Retention ==\n")
  fs <- file.path(OUT_DIR, sprintf("59_batch_%02d.csv", seq_len(nb)))
  ALL <- do.call(rbind, lapply(fs[file.exists(fs)], read.csv, stringsAsFactors = FALSE))
  write.csv(ALL, file.path(OUT_DIR, "59_jitter_sens_replicates.csv"), row.names = FALSE)
  base_ret <- S50$retention[S50$type == "base"]; prim_ret <- S50$retention[S50$type == "jitter"]
  SUMR <- do.call(rbind, lapply(names(SC), function(sn) {
    d <- ALL[ALL$scn == sn & !is.na(ALL$p), ]; k <- sum(d$p < TESTS$alpha); n <- nrow(d); ci <- wilson(k, n); q <- quantile(d$p, c(0.1, 0.25, 0.5, 0.75, 0.9))
    data.frame(scenario = sn, author_m = SC[[sn]][["AUTHOR_REPORTED"]], gis_m = SC[[sn]][["GIS_MAP"]], reps = n, failed = sum(ALL$scn == sn) - n, retained = k,
               retention = k / n, ci_lo = ci[1], ci_hi = ci[2], p_q10 = q[1], p_q25 = q[2], p_median = q[3], p_q75 = q[4], p_q90 = q[5], min_usable = min(d$usable),
               pair_ratio_2 = median(d$po2 / d$ps2), pair_ratio_3 = median(d$po3 / d$ps3), pair_ratio_5 = median(d$po5 / d$ps5),
               mean_fallbacks = mean(d$fallbacks), verdict = if (k / n >= 0.80) "frozen criterion MET (>= 0.80)" else "frozen criterion NOT MET (< 0.80)")
  }))
  print(SUMR, digits = 3, row.names = FALSE)
  cat(sprintf("\nreference (script 50, nsim 199): baseline retention %.3f | primary jitter (AUTHOR 25 m, GIS 250 m) retention %.3f\n", base_ret, prim_ret))
  SUMR$baseline_retention <- base_ret; SUMR$primary_retention <- prim_ret; SUMR$date <- format(Sys.time())
  write.csv(SUMR, f_res, row.names = FALSE)
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink(); if (exists("cl")) try(stopCluster(cl), silent = TRUE) })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "59_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")