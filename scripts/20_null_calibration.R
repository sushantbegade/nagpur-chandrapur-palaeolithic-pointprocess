# =============================================================================
# 20_null_calibration.R   (v3: workers load a cache, no terra on workers)
# Stage: 10_simulation_recovery | validates the adopted fast pipeline and the M1 test
# Purpose: (1) write data_derived/recovery_cache.rds, R/recovery_worker.R and
#   R/recovery_functions.R; (2) fast-fit re-check and pipeline timing; (3) one-worker
#   diagnostic (19 simulations); (4) null calibration of the M1 test (parallel);
#   (5) cost re-estimate and batch schedule. SYNTHETIC patterns only.
# Requires R/recovery_setup.R (script 17).
# Outputs: data_derived/recovery_cache.rds, R/recovery_worker.R, R/recovery_functions.R,
#   10_simulation_recovery/outputs/null_calibration_report.txt, null_calibration.csv,
#   cost_estimate_v3.csv, sessionInfo_20.txt ; docs/protocol_deviation_log.csv (note row)
# Not executed by assistant. RESTART R FIRST, then run top to bottom (20-30 min).
# WARNING: first line clears the objects in your R workspace (not files).
# =============================================================================

rm(list = ls(envir = globalenv()), envir = globalenv()); invisible(gc())

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
N_NULL_REPS <- 24L
NSIM_NULL   <- 99L
DIAG_NSIM   <- 19L        # minimum for a GET test at alpha 0.05

WORKER_TEXT <- r"---(
# recovery_worker.R -- written by 20_null_calibration.R. Light loader for parallel workers (no terra).
# Requires REPO_DIR from R/config.R and data_derived/recovery_cache.rds.
params_f <- file.path(REPO_DIR, "R", "freeze_v1_params.R"); stopifnot(file.exists(params_f)); source(params_f)
if (!grepl("^FROZEN", FREEZE_STATUS)) stop("Freeze v1 is not FROZEN. Aborting.")
suppressPackageStartupMessages({ library(spatstat.geom); library(spatstat.explore); library(spatstat.model)
                                 library(spatstat.random); library(GET) })
cache <- readRDS(file.path(REPO_DIR, "data_derived", "recovery_cache.rds"))
covs <- cache$covs; W_km <- cache$W_km; X_obs_km <- cache$X_obs_km; formula_E3 <- cache$formula_E3
QUAD_KM <- QUADRATURE_KM$primary
RGRID <- seq(0, R_MAX_DECISION_KM, by = R_STEP_KM)
rL <- RGRID[RGRID >= R_MIN_KM["A"] & RGRID <= R_MAX_DECISION_KM]
rG <- RGRID[RGRID >= GINHOM_LOWER_KM["A"] & RGRID <= R_MAX_DECISION_KM]
to_poly <- function(P, label = "") {
  Q <- P[W_km]
  if (Window(Q)$type != "polygonal") stop("pattern window is not polygonal after restriction")
  Q
}
col_of <- function(f, pref) { for (nm in pref) if (nm %in% names(f)) return(f[[nm]]); stop("no column among: ", paste(names(f), collapse = ",")) }
make_truth <- function(env_mult = 1, n_target = RECOVERY_ENV$expected_n) {
  set.seed(SEEDS$master + RECOVERY_ENV$seed_offset)
  beta <- setNames(rnorm(length(covs), 0, RECOVERY_ENV$coef_sd), names(covs)) * env_mult
  eta <- Reduce(`+`, Map(function(a, b) a * b, covs, beta))
  lam0 <- exp(eta)[W_km, drop = FALSE]
  b0 <- log(n_target / integral(lam0))
  list(beta = beta, b0 = b0, lamW = lam0 * exp(b0))
}
)---"

FUNCTIONS_TEXT <- r"---(
# recovery_functions.R -- written by 20_null_calibration.R; sourced by recovery scripts and workers.
# Requires covs, W_km, X_obs_km, QUAD_KM, RGRID, rL, rG, col_of, to_poly, make_truth, formula_E3
# (from recovery_setup.R in the main session, or recovery_worker.R on workers).
covsW <- lapply(covs, function(im) im[W_km, drop = FALSE])
area_km2 <- area.owin(W_km)
bw_for_n <- function(n) (PCF_STOYAN_C / sqrt(n / area_km2)) / sqrt(5)   # kernel sd; the Stoyan value is the half-width
L_ARGS    <- list(renormalise = TRUE, normpower = 1)
G_ARGS_NS <- list(kernel = "epanechnikov", zerocor = "convolution", divisor = "r", renormalise = TRUE, normpower = 2)
DUMMY_FIXED <- default.dummy(X_obs_km, eps = QUAD_KM)
truth <- make_truth()
make_quad <- function(Xp) quadscheme(Xp, dummy = DUMMY_FIXED, method = "grid")
fast_fit <- function(Xp) {
  Q <- make_quad(Xp)
  w <- w.quad(Q); isd <- is.data(Q); xq <- x.quad(Q); yq <- y.quad(Q)
  Z <- vapply(covs, function(im) lookup.im(im, xq, yq, naok = TRUE), numeric(length(xq)))
  ok <- stats::complete.cases(Z) & w > 0
  X <- cbind(`(Intercept)` = 1, Z)
  gfit <- suppressWarnings(stats::glm.fit(X[ok, , drop = FALSE], (isd / w)[ok], weights = w[ok], family = poisson(), intercept = FALSE))
  beta <- gfit$coefficients
  Xd <- X[isd & ok, , drop = FALSE]
  if (nrow(Xd) != npoints(Xp)) stop("data points dropped in fast_fit: ", nrow(Xd), " of ", npoints(Xp))
  list(beta = beta, lam_data = as.numeric(exp(Xd %*% beta)))
}
geyer_fit <- function(Xp, r, sat, rbord = 10)
  ppm(make_quad(Xp), trend = formula_E3, interaction = Geyer(r = r, sat = sat), covariates = covs, rbord = rbord)
lam_im_of <- function(beta) exp(beta[1] + Reduce(`+`, Map(function(a, b) a * b, covsW, beta[-1])))
fast_curves <- function(Xp, ftb, bw) {
  Lh <- do.call(Linhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "isotropic"), L_ARGS))
  gh <- do.call(pcfinhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "Ripley", bw = bw), G_ARGS_NS))
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(rG, 6), round(gh$r, 6))])
}
m1_test <- function(Xobs, nsim, bw, seed) {
  if (nsim < 19) stop("GET ERL test at alpha 0.05 needs nsim >= 19")
  set.seed(seed)
  f0 <- fast_fit(Xobs); ob <- fast_curves(Xobs, f0, bw)
  lam_im <- lam_im_of(f0$beta)
  sL <- matrix(NA_real_, length(rL), nsim); sG <- matrix(NA_real_, length(rG), nsim)
  for (k in seq_len(nsim)) {
    Xs <- to_poly(rpoispp(lam_im)); fs <- fast_fit(Xs); cs <- fast_curves(Xs, fs, bw)
    sL[, k] <- cs$L; sG[, k] <- cs$g
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  res <- GET::global_envelope_test(curve_sets = list(
    GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE])),
    GET::create_curve_set(list(r = rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))), type = "erl", alpha = TESTS$alpha)
  list(p = as.numeric(attr(res, "p")), usable = sum(okc))
}
null_rep <- function(i, nsim) {
  t0 <- unname(proc.time()["elapsed"])
  set.seed(SEEDS$master + SEEDS$offsets["empirical"] + 1000 + i)
  X <- to_poly(rpoispp(truth$lamW))
  r <- m1_test(X, nsim, bw_for_n(RECOVERY_ENV$expected_n), SEEDS$master + SEEDS$offsets["empirical"] + 2000 + i)
  c(p = r$p, usable = r$usable, n = npoints(X), sec = unname(proc.time()["elapsed"]) - t0)
}
)---"

while (sink.number() > 0) sink()
source(REPO_CONFIG)
source(file.path(REPO_DIR, "R", "recovery_setup.R"))
worker_f <- file.path(REPO_DIR, "R", "recovery_worker.R")
func_f   <- file.path(REPO_DIR, "R", "recovery_functions.R")
cache_f  <- file.path(REPO_DIR, "data_derived", "recovery_cache.rds")
writeLines(WORKER_TEXT, worker_f); writeLines(FUNCTIONS_TEXT, func_f)
saveRDS(list(covs = covs, W_km = W_km, X_obs_km = X_obs_km, formula_E3 = formula_E3), cache_f)
source(func_f)
out <- dir_of("10_simulation_recovery"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
maxdiff <- function(a, b) suppressWarnings(max(abs(a - b), na.rm = TRUE))
M <- list()

sink(file.path(out, "null_calibration_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== 0. MEMORY ====\n")
  g0 <- gc(); cat(sprintf("R memory in use: %.0f MB | cache file: %.1f MB | covariate images: %s | window object: %s\n",
                          sum(g0[, 2]), file.size(cache_f) / 1e6, format(object.size(covs), units = "MB"), format(object.size(W_km), units = "MB")))
  cat("physical cores:", parallel::detectCores(logical = FALSE), "\n")
  
  cat("\n==== A. FAST FIT RE-CHECK ====\n")
  set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 60)
  Xa <- to_poly(rpoispp(truth$lamW)); fa <- fast_fit(Xa)
  pa <- ppm(Xa, trend = formula_E3, covariates = covs, quad.args = list(eps = QUAD_KM))
  cat(sprintf("Poisson pattern n = %d | max |coef diff| fast_fit vs ppm: %.2e\n", npoints(Xa), maxdiff(fa$beta, coef(pa)[names(fa$beta)])))
  tp <- numeric(5); for (i in 1:5) tp[i] <- system.time(fast_fit(Xa))["elapsed"]; M$fast_fit <- median(tp)
  bwn <- bw_for_n(RECOVERY_ENV$expected_n)
  tc <- numeric(5); for (i in 1:5) tc[i] <- system.time(fast_curves(Xa, fa, bwn))["elapsed"]; M$curves <- median(tc)
  M$pipe <- M$fast_fit + M$curves; M$geyer_fit <- 1.5
  cat(sprintf("per-pattern pipeline: fit %.2f s + curves %.2f s = %.2f s\n", M$fast_fit, M$curves, M$pipe))
  
  mk_cluster <- function(n) {
    cl <- parallel::makeCluster(n)
    parallel::clusterExport(cl, "REPO_CONFIG", envir = globalenv())
    msg <- parallel::clusterEvalQ(cl, tryCatch({
      source(REPO_CONFIG); source(file.path(REPO_DIR, "R", "recovery_worker.R")); source(file.path(REPO_DIR, "R", "recovery_functions.R"))
      paste0("ok (", round(sum(gc()[, 2])), " MB)")
    }, error = function(e) paste("SETUP ERROR:", conditionMessage(e), "| R memory", round(sum(gc()[, 2])), "MB")))
    list(cl = cl, msg = unlist(msg))
  }
  
  cat("\n==== B. PARALLEL SETUP DIAGNOSTIC (one worker, ", DIAG_NSIM, " simulations) ====\n", sep = "")
  t_setup <- system.time(c1 <- mk_cluster(1))["elapsed"]
  cat(sprintf("worker setup %.1f s | worker message: %s\n", t_setup, paste(c1$msg, collapse = " | ")))
  d1 <- try(parallel::clusterCall(c1$cl, function(n) tryCatch({ x <- null_rep(1, n); paste(names(x), round(x, 2), sep = "=", collapse = " ") },
                                                              error = function(e) paste("ERROR:", conditionMessage(e), "| call:", paste(deparse(conditionCall(e))[1], collapse = " "))), DIAG_NSIM), silent = TRUE)
  parallel::stopCluster(c1$cl)
  diag_msg <- if (inherits(d1, "try-error")) conditionMessage(attr(d1, "condition")) else unlist(d1)[1]
  cat("one-worker test replicate:", diag_msg, "\n")
  if (grepl("ERROR", diag_msg)) stop("Worker diagnostic failed. See the message above.")
  
  cat("\n==== C. NULL CALIBRATION OF THE M1 TEST (", N_NULL_REPS, " Poisson-truth replicates x ", NSIM_NULL, " simulations) ====\n", sep = "")
  ncl <- min(4, parallel::detectCores(logical = FALSE))
  cc <- mk_cluster(ncl)
  cat("workers:", ncl, "| setup messages:", paste(unique(cc$msg), collapse = " | "), "\n")
  tw <- system.time(pr <- parallel::parLapply(cc$cl, seq_len(N_NULL_REPS), function(i, nsim) tryCatch(null_rep(i, nsim), error = function(e) paste("ERROR:", conditionMessage(e))), nsim = NSIM_NULL))["elapsed"]
  parallel::stopCluster(cc$cl)
  is_err <- vapply(pr, is.character, logical(1))
  if (any(is_err)) cat("worker errors:", paste(unique(unlist(pr[is_err])), collapse = " | "), "\n")
  M$par_speedup <- NA_real_
  if (any(!is_err)) {
    P <- do.call(rbind, pr[!is_err]); pv <- P[, "p"]
    cat(sprintf("replicates completed %d of %d | wall time %.1f min | mean worker time per replicate %.1f s | usable simulations per test: min %d\n",
                nrow(P), N_NULL_REPS, tw / 60, mean(P[, "sec"]), min(P[, "usable"])))
    cat(sprintf("synthetic pattern sizes: mean n %.1f (min %d, max %d)\n", mean(P[, "n"]), min(P[, "n"]), max(P[, "n"])))
    cat(sprintf("p-values: rejections at 0.05: %d of %d (expect about %.1f) | p <= 0.10: %d | median p %.2f | min p %.3f\n",
                sum(pv <= 0.05), length(pv), 0.05 * length(pv), sum(pv <= 0.10), median(pv), min(pv)))
    print(round(sort(pv), 3))
    M$par_speedup <- sum(P[, "sec"]) / tw
    cat(sprintf("parallel efficiency: %.2f workers' worth on %d workers\n", M$par_speedup, ncl))
    write.csv(data.frame(rep = which(!is_err), p = pv, n = P[, "n"], usable = P[, "usable"], sec = P[, "sec"]), file.path(out, "null_calibration.csv"), row.names = FALSE)
  }
  
  cat("\n==== D. COST RE-ESTIMATE (hierarchy fractions ASSUMED; M3 and M4 parameter draws, no Gibbs refit per simulation) ====\n")
  kfit <- 0.34; ksim <- 0.34; rmh_ps <- 0.61; rmh_ov <- 10.8       # measured in script 19
  tgrid <- 8 * M$geyer_fit
  cost <- function(nsim, reps) {
    m1 <- nsim * (0.02 + M$pipe); m3 <- nsim * (ksim + kfit + M$pipe); m4 <- nsim * (rmh_ps + M$pipe) + rmh_ov
    fit_all <- M$pipe + kfit + tgrid
    pois <- fit_all + m1 + 0.1 * m3; lgcp <- fit_all + m1 + 0.9 * m3 + 0.1 * m4; gey <- fit_all + m1 + 0.9 * m3 + 0.7 * m4
    per_set <- 6 * pois + 7 * lgcp + 6 * gey
    c(hours_per_replicate_index = per_set / 3600, study_single_core_h = reps * per_set / 3600)
  }
  sp <- if (!is.null(M$par_speedup) && is.finite(M$par_speedup) && M$par_speedup > 1) M$par_speedup else 2.3
  est <- rbind(data.frame(workflow = "frozen counts (200 reps)", nsim = TESTS$nsim_recovery, reps = TESTS$recovery_reps, t(cost(TESTS$nsim_recovery, TESTS$recovery_reps))),
               data.frame(workflow = "199 sims, 100 reps (after batch 5)", nsim = TESTS$nsim_recovery, reps = 100, t(cost(TESTS$nsim_recovery, 100))),
               data.frame(workflow = "floors (100 reps, 99 sims)", nsim = TESTS$nsim_recovery_floor, reps = TESTS$recovery_reps_floor, t(cost(TESTS$nsim_recovery_floor, TESTS$recovery_reps_floor))))
  est$days_on_workers <- round(est$study_single_core_h / sp / 24, 1); est[, 4:5] <- round(est[, 4:5], 1)
  print(est, row.names = FALSE)
  cat(sprintf("parallel speed-up used: x%.2f (%s)\n", sp, if (sp == 2.3) "default from script 18; calibration did not provide one" else "measured here"))
  bat <- cost(TESTS$nsim_recovery, 20)["study_single_core_h"]
  cat(sprintf("one batch of 20 replicates across all 19 cells: %.1f single-core hours = %.1f hours on %d workers; 10 batches reach 200 replicates\n", bat, bat / sp, ncl))
  cat("Measured components (s):", paste(names(M), round(unlist(M), 3), sep = "=", collapse = " | "), "\n")
  write.csv(est, file.path(out, "cost_estimate_v3.csv"), row.names = FALSE)
  
  dl_f <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv")
  dl <- read.csv(dl_f, stringsAsFactors = FALSE)
  note <- data.frame(date = format(Sys.Date()), item = "NOTE null calibration of the M1 test is reported by script 20 (script 19 section D failed: NSIM_NULL not exported to workers)",
                     reason = "Correction of the log text of the third clarification row.", effect_on_inference = "None.", approved_by = "author-delegated (assistant)", stringsAsFactors = FALSE)
  if (!note$item %in% dl$item) { write.csv(rbind(dl[, names(note)], note), dl_f, row.names = FALSE); cat("\nAppended correcting note to docs/protocol_deviation_log.csv\n") }
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_20.txt"))
cat("\nDONE. Outputs in:", out, "\n")