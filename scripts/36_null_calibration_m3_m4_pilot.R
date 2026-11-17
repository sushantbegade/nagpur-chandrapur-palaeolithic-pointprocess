# =====================================================================
# 36_null_calibration_m3_m4_pilot.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing;
#           writes R/recovery_tests.R, which later recovery scripts source)
# Purpose : (0) test whether NA-covariate points arise in boundary pixels
#               (1e6 uniform points in W);
#           (1) worker-path smoke test (nsim 19);
#           (2) null-calibration pilot with 4 workers:
#               M3 test on L4 truth (sigma2 1.5, scale 8 km, rLGCP eps 0.5, no Gibbs);
#               M4 test on G3 truth fitted at the true (r = 3, sat = 3), rbord 0;
#               N_REP reps each, nsim = NSIM, checkpointed in batches;
#           (3) size summary and amendment of log row "NA-covariate points removed".
# Outputs : outputs/36_report.txt, 36_batch_NN.csv, 36_null_pilot.csv,
#           36_na_boundary_probe.csv, 36_sessionInfo.txt,
#           protocol_deviation_log_before_36.csv; docs/protocol_deviation_log.csv (amended)
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           Pilot resolution: 40 reps detect only gross size problems.
#           M4 is calibrated at the true (r, sat); the profile search over
#           the (r, sat) grid is NOT part of this test.
#           Pre-set rule (before running): per test, ANTICONSERVATIVE if the
#           lower 95% binomial bound of the rejection rate at 0.05 exceeds 0.05;
#           CONSERVATIVE if mean p > 0.60; FAILURES if > 10% of reps give no p.
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

OUT_DIR <- dir_of("10_simulation_recovery", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
TESTS_FILE <- file.path(REPO_DIR, "R", "recovery_tests.R")

SEED9      <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 8000
N_REP      <- 40
BATCH      <- 8
NSIM       <- 99
N_WORK     <- 4
LGCP_VAR   <- 1.5
LGCP_SCALE <- 8
EPS_LGCP   <- 0.5
NREP_RMH   <- 5e5
N_UNI      <- 1e6
N_CHUNK    <- 5

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULT_G3 <- gm$mult[gm$cell == "G3"]; stopifnot(length(MULT_G3) == 1)
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")

# ---- test functions (dumped to R/recovery_tests.R and sourced by workers) ----
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
fill_im <- function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
trend_im <- function(b) fill_im(Reduce(`+`, Map(function(a, nmc) a * unname(b[nmc]), covs, names(covs))) + unname(b["(Intercept)"]))
lgcp_mu <- function() fill_im(Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta)) + truth$b0)
na_cov_pts <- function(P) {
  if (npoints(P) == 0) return(logical(0))
  Z <- do.call(cbind, lapply(covs, function(im) lookup.im(im, P$x, P$y, naok = TRUE)))
  rowSums(is.na(Z)) > 0
}
clean_pat <- function(P) { bad <- na_cov_pts(P); if (any(bad)) P[!bad] else P }
sim_geyer <- function(mult, gamma, r, sat, n_start, nsim, seed, nrep = NREP_RMH) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}
curves_test <- function(Xobs, sims, bw) {
  Xo <- clean_pat(Xobs)
  f0 <- fast_fit(Xo); ob <- fast_curves(Xo, f0, bw)
  sL <- matrix(NA_real_, length(rL), length(sims)); sG <- matrix(NA_real_, length(rG), length(sims))
  dsim <- integer(length(sims))
  for (k in seq_along(sims)) {
    Xr <- to_poly(sims[[k]]); Xs <- clean_pat(Xr); dsim[k] <- npoints(Xr) - npoints(Xs)
    fs <- fast_fit(Xs); cs <- fast_curves(Xs, fs, bw)
    sL[, k] <- cs$L; sG[, k] <- cs$g
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc), drop_sim = sum(dsim), sims_dropped = sum(dsim > 0))
}
m3_test <- function(Xobs, nsim, bw, seed) {
  if (nsim < 19) stop("ERL test at alpha 0.05 needs nsim >= 19")
  set.seed(seed)
  Xc <- clean_pat(Xobs)
  fit <- tryCatch(kppm(Xc, trend = formula_E3, clusters = "LGCP", model = "exponential",
                       covariates = covs, method = "clik2"), error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA_real_, usable = 0L, drop_sim = NA_integer_, sims_dropped = NA_integer_,
                                sigma2 = NA_real_, alpha = NA_real_))
  s2 <- unname(fit$par["sigma2"]); al <- unname(fit$par["alpha"])
  eta <- trend_im(coef(fit)); n <- npoints(Xc)
  S <- rLGCP(model = "exp", mu = eta - s2 / 2, param = list(var = s2, scale = al), win = W_km,
             eps = EPS_LGCP, saveLambda = TRUE, nsim = nsim, drop = TRUE)
  sims <- lapply(as_list(S), function(P) rpoint(n, f = attr(P, "Lambda"), win = W_km, forcewin = TRUE))
  r <- curves_test(Xc, sims, bw)
  c(r, list(sigma2 = s2, alpha = al))
}
m4_test <- function(Xobs, r_km, sat, nsim, bw, seed) {
  if (nsim < 19) stop("ERL test at alpha 0.05 needs nsim >= 19")
  set.seed(seed)
  Xc <- clean_pat(Xobs)
  fit <- tryCatch(geyer_fit(Xc, r = r_km, sat = sat, rbord = 0), error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA_real_, usable = 0L, drop_sim = NA_integer_, sims_dropped = NA_integer_, gamma = NA_real_))
  gh <- exp(unname(tail(coef(fit), 1))); n <- npoints(Xc)
  S <- simulate(fit, nsim = nsim, drop = TRUE, verbose = FALSE, progress = FALSE,
                control = rmhcontrol(nrep = NREP_RMH, p = 1, expand = 1), start = rmhstart(n.start = n))
  r <- curves_test(Xc, as_list(S), bw)
  c(r, list(gamma = gh))
}
null_task <- function(job) {
  t0 <- unname(proc.time()["elapsed"])
  base <- SEEDS$master + unname(SEEDS$offsets["empirical"])
  bw <- bw_for_n(RECOVERY_ENV$expected_n)
  out <- tryCatch({
    if (job$test == "M3") {
      set.seed(base + 3000 + job$i)
      X <- clean_pat(to_poly(rLGCP(model = "exp", mu = lgcp_mu() - LGCP_VAR / 2,
                                   param = list(var = LGCP_VAR, scale = LGCP_SCALE), win = W_km,
                                   eps = EPS_LGCP, nsim = 1, drop = TRUE)))
      r <- m3_test(X, job$nsim, bw, base + 4000 + job$i)
      list(n = npoints(X), p = r$p, usable = r$usable, e1 = r$sigma2, e2 = r$alpha, dsim = r$drop_sim, err = NA_character_)
    } else {
      X <- to_poly(sim_geyer(job$mult, 2.5, 3, 3, 197, 1, base + 5000 + job$i))
      r <- m4_test(X, 3, 3, job$nsim, bw, base + 6000 + job$i)
      list(n = npoints(X), p = r$p, usable = r$usable, e1 = r$gamma, e2 = NA_real_, dsim = r$drop_sim, err = NA_character_)
    }
  }, error = function(e) list(n = NA_real_, p = NA_real_, usable = 0L, e1 = NA_real_, e2 = NA_real_,
                              dsim = NA_integer_, err = conditionMessage(e)))
  c(list(test = job$test, i = job$i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
init_worker <- function(repo, cache, tests) {
  assign("REPO_DIR", repo, envir = globalenv())
  Sys.setenv(RECOVERY_CACHE = cache)
  source(file.path(repo, "R", "recovery_worker.R"), local = globalenv())
  source(file.path(repo, "R", "recovery_functions.R"), local = globalenv())
  source(tests, local = globalenv())
  TRUE
}
rows_df <- function(L) do.call(rbind, lapply(L, function(x) data.frame(
  test = x$test, i = x$i, sec = as.numeric(x$sec), n = as.numeric(x$n), p = as.numeric(x$p),
  usable = as.numeric(x$usable), e1 = as.numeric(x$e1), e2 = as.numeric(x$e2),
  dsim = as.numeric(x$dsim), err = as.character(x$err), stringsAsFactors = FALSE)))
summ_test <- function(d, label) {
  ok <- !is.na(d$p); p <- d$p[ok]; k <- length(p)
  cat(sprintf("\n-- %s --\n", label))
  cat(sprintf("reps %d | valid p %d | failed %d%s\n", nrow(d), k, nrow(d) - k,
              if (any(!is.na(d$err))) paste0(" | errors: ", paste(unique(d$err[!is.na(d$err)]), collapse = " ; ")) else ""))
  if (k < 5) { cat("too few valid p-values\n"); return(invisible("NO VERDICT")) }
  r05 <- sum(p < 0.05); r10 <- sum(p < 0.10)
  c05 <- binom.test(r05, k)$conf.int; c10 <- binom.test(r10, k)$conf.int
  ksp <- suppressWarnings(ks.test(p, "punif")$p.value)
  cat(sprintf("mean p %.3f (uniform: 0.5) | share p > 0.5: %.2f | KS p %.3f (ties, indicative)\n", mean(p), mean(p > 0.5), ksp))
  cat(sprintf("rejections at 0.05: %d/%d = %.3f (95%% CI %.3f-%.3f) | at 0.10: %d/%d = %.3f (CI %.3f-%.3f)\n",
              r05, k, r05 / k, c05[1], c05[2], r10, k, r10 / k, c10[1], c10[2]))
  cat(sprintf("seconds per test: median %.0f, max %.0f | n: median %.0f (range %.0f-%.0f) | min usable sims %d | synthetic points dropped (cleaning): %d\n",
              median(d$sec), max(d$sec), median(d$n, na.rm = TRUE), min(d$n, na.rm = TRUE), max(d$n, na.rm = TRUE),
              as.integer(min(d$usable)), as.integer(sum(d$dsim, na.rm = TRUE))))
  cat(sprintf("estimates: e1 median %.3f (IQR %.3f-%.3f) | e2 median %s\n", median(d$e1, na.rm = TRUE),
              quantile(d$e1, 0.25, na.rm = TRUE), quantile(d$e1, 0.75, na.rm = TRUE),
              if (all(is.na(d$e2))) "NA" else sprintf("%.2f (IQR %.2f-%.2f)", median(d$e2, na.rm = TRUE),
                                                      quantile(d$e2, 0.25, na.rm = TRUE), quantile(d$e2, 0.75, na.rm = TRUE))))
  v <- if ((nrow(d) - k) / nrow(d) > 0.10) "FAILURES > 10%" else if (c05[1] > 0.05) "ANTICONSERVATIVE"
  else if (mean(p) > 0.60) "CONSERVATIVE (power cost)" else "no gross size problem detected"
  cat("verdict (pre-set rule):", v, "\n")
  invisible(v)
}

# ---- write the tests file for later scripts and workers ----
dump(c("LGCP_VAR", "LGCP_SCALE", "EPS_LGCP", "NREP_RMH", "as_list", "fill_im", "trend_im", "lgcp_mu", "na_cov_pts",
       "clean_pat", "sim_geyer", "curves_test", "m3_test", "m4_test", "null_task"), file = TESTS_FILE)

WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }

REPORT <- file.path(OUT_DIR, "36_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("36_null_calibration_m3_m4_pilot | ", format(Sys.time()), "\n", sep = "")
    cat(sprintf("N_REP %d per test | batch %d | nsim %d | workers %d | tests file %s\n", N_REP, BATCH, NSIM, N_WORK, TESTS_FILE))
    
    ## ---------------- 0. NA-covariate boundary-pixel probe ----------------
    cat("\n== 0. NA-covariate points among uniform points in W (boundary-pixel hypothesis) ==\n")
    set.seed(SEED9 + 1)
    tot <- 0; nb <- 0; bdv <- numeric(0)
    for (ch in seq_len(N_CHUNK)) {
      U <- runifpoint(N_UNI / N_CHUNK, win = W_km)
      bad <- na_cov_pts(U)
      tot <- tot + npoints(U); nb <- nb + sum(bad)
      if (any(bad)) bdv <- c(bdv, bdist.points(U[bad]))
    }
    cat(sprintf("uniform points %d | NA-covariate points %d (rate %.2e per point; observed in script 35 rerun: about 1 per 7e4)\n", tot, nb, nb / tot))
    if (nb > 0) cat(sprintf("boundary distance (km) of those points: min %.3f median %.3f max %.3f -> hypothesis %s\n",
                            min(bdv), median(bdv), max(bdv), if (median(bdv) < 1) "SUPPORTED (boundary pixels)" else "NOT supported (interior)"))
    else cat("none found: boundary-pixel hypothesis not confirmed by this probe (rate below about 1e-6)\n")
    write.csv(data.frame(uniform_points = tot, na_points = nb, rate = nb / tot,
                         bd_median_km = if (nb) median(bdv) else NA_real_), file.path(OUT_DIR, "36_na_boundary_probe.csv"), row.names = FALSE)
    
    ## ---------------- cluster ----------------
    cl <- makeCluster(N_WORK)
    ini <- clusterCall(cl, init_worker, REPO_DIR, CACHE_FILE, TESTS_FILE)
    cat("\nworkers initialised:", all(unlist(ini)), "\n")
    ex <- clusterCall(cl, function() c(exists("m3_test"), exists("m4_test"), exists("truth"), "terra" %in% loadedNamespaces()))
    cat("worker check (m3_test, m4_test, truth present; terra loaded):", paste(ex[[1]], collapse = " "), "\n")
    
    ## ---------------- 1. smoke test ----------------
    cat("\n== 1. Smoke test through the worker path (nsim 19, i = 900) ==\n")
    sm <- parLapply(cl, list(list(test = "M3", i = 900, nsim = 19, mult = NA_real_),
                             list(test = "M4", i = 900, nsim = 19, mult = MULT_G3)), null_task)
    smd <- rows_df(sm); print(smd, digits = 3, row.names = FALSE)
    if (all(!is.na(smd$err))) stop("smoke test failed for both tests; paste output")
    
    ## ---------------- 2. batches ----------------
    cat("\n== 2. Null-calibration batches ==\n")
    nb_batches <- ceiling(N_REP / BATCH)
    for (b in seq_len(nb_batches)) {
      f_b <- file.path(OUT_DIR, sprintf("36_batch_%02d.csv", b))
      if (file.exists(f_b)) { cat(sprintf("  batch %d: found %s, skipped\n", b, basename(f_b))); next }
      ii <- ((b - 1) * BATCH + 1):min(b * BATCH, N_REP)
      jobs <- c(lapply(ii, function(i) list(test = "M3", i = i, nsim = NSIM, mult = NA_real_)),
                lapply(ii, function(i) list(test = "M4", i = i, nsim = NSIM, mult = MULT_G3)))
      t0 <- unname(proc.time()["elapsed"])
      res <- parLapplyLB(cl, jobs, null_task)
      db <- rows_df(res); db$batch <- b; db$batch_wall_sec <- unname(proc.time()["elapsed"]) - t0
      write.csv(db, f_b, row.names = FALSE)
      cat(sprintf("  batch %d: %d tests in %.0f s wall | M3 p: %s | M4 p: %s\n", b, nrow(db), unname(proc.time()["elapsed"]) - t0,
                  paste(sprintf("%.2f", db$p[db$test == "M3"]), collapse = " "),
                  paste(sprintf("%.2f", db$p[db$test == "M4"]), collapse = " ")))
    }
    
    ## ---------------- 3. summary ----------------
    fs <- file.path(OUT_DIR, sprintf("36_batch_%02d.csv", seq_len(nb_batches)))
    ALL <- do.call(rbind, lapply(fs[file.exists(fs)], read.csv, stringsAsFactors = FALSE))
    write.csv(ALL, file.path(OUT_DIR, "36_null_pilot.csv"), row.names = FALSE)
    cat("\n== 3. Size summary (pilot; resolution coarse) ==\n")
    summ_test(ALL[ALL$test == "M3", ], "M3 test, L4 truth (e1 = sigma2_hat, e2 = alpha_hat)")
    summ_test(ALL[ALL$test == "M4", ], "M4 test, G3 truth at true (r, sat) (e1 = gamma_hat)")
    wall <- tapply(ALL$batch_wall_sec, ALL$batch, function(x) x[1])
    cat(sprintf("\nwall seconds per batch (8 M3 + 8 M4 tests): %s | mean %.0f s => %.1f s wall per test, 4 workers\n",
                paste(sprintf("%.0f", wall), collapse = " "), mean(wall), mean(wall) / (2 * BATCH)))
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink(); if (exists("cl")) try(stopCluster(cl), silent = TRUE) })

## ---------------- log amendment (after sink closed) ----------------
idx <- grep("^Clarification: NA-covariate points removed", LOG0$item)
if (length(idx) == 1 && identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_36.csv"), overwrite = TRUE)
  LOG1 <- LOG0
  LOG1$reason[idx] <- paste("Script 34 M3 rep 4 failed: fast_fit dropped a data point (214 of 215). Script 35 found 0 NA-covariate pixel centres",
                            "inside the window and 0 NA-covariate points in 24,722 synthetic points, so the mechanism is NOT established",
                            "(hypothesis, tested in script 36: points in boundary pixels whose centres lie outside the window).",
                            "Rerun with cleaning completed 4 of 4; one simulated point was dropped in 1 of 396 simulations.")
  LOG1$effect_on_inference[idx] <- paste("Rare synthetic points (about 1 per 7e4 simulated points) are removed before fitting; counts are logged per test.",
                                         "Observed archaeological sites have 0 NA covariates (script 35). No archaeological model has been fitted.")
  write.csv(LOG1, f_log, row.names = FALSE)
  cat("Amended log row", idx, "(backup protocol_deviation_log_before_36.csv)\n")
} else cat("Log row not amended (match count", length(idx), "or columns differ); amend manually.\n")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "36_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")