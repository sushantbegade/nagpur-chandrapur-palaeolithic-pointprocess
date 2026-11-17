# =============================================================================
# 21_window_speed_test.R
# Stage: 10_simulation_recovery | window fallback test, real parallel speed-up, ERL variants
# Purpose: (1) rewrite R/recovery_functions.R and R/recovery_worker.R (named cache, both
#   ERL combination variants, timing probe); (2) simplified-window check and cache;
#   (3) per-pattern timing, exact vs 100 m simplified window; (4) parallel scaling on
#   1/2/4 workers with identical task lists; (5) null calibration, 48 replicates x 99
#   simulations, nstep 2 vs 1; (6) cost re-estimate. SYNTHETIC patterns only.
# Does NOT adopt any window or any ERL variant by itself.
# Requires R/recovery_setup.R (script 17) and data_derived/recovery_cache.rds (script 20 v3).
# Outputs (10_simulation_recovery/outputs/): window_speed_report.txt, window_variant_check.csv,
#   parallel_scaling.csv, null_calibration_v2.csv, cost_estimate_v4.csv, sessionInfo_21.txt
# Not executed by assistant. RESTART R FIRST, then run top to bottom (40-60 min).
# WARNING: first line clears the objects in your R workspace (not files).
# =============================================================================

rm(list = ls(envir = globalenv()), envir = globalenv()); invisible(gc())

REPO_CONFIG   <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
N_NULL_REPS   <- 48L
NSIM_NULL     <- 99L
PAR_TASKS     <- 8L
PAR_NSIM      <- 19L
PAR_WORKERS   <- c(1L, 2L, 4L)
PROBE_N       <- 10L
SEED_OFFSET   <- 100L      # null replicate indices 101..148 (script 20 used 1..24)

WORKER_TEXT <- r"---(
# recovery_worker.R -- written by 21_window_speed_test.R. Light loader for parallel workers (no terra).
# Cache file name comes from environment variable RECOVERY_CACHE (default recovery_cache.rds).
params_f <- file.path(REPO_DIR, "R", "freeze_v1_params.R"); stopifnot(file.exists(params_f)); source(params_f)
if (!grepl("^FROZEN", FREEZE_STATUS)) stop("Freeze v1 is not FROZEN. Aborting.")
suppressPackageStartupMessages({ library(spatstat.geom); library(spatstat.explore); library(spatstat.model)
                                 library(spatstat.random); library(GET) })
CACHE_NAME <- Sys.getenv("RECOVERY_CACHE", "recovery_cache.rds")
cache <- readRDS(file.path(REPO_DIR, "data_derived", CACHE_NAME))
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
# recovery_functions.R -- written by 21_window_speed_test.R; sourced by recovery scripts and workers.
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
erl_p <- function(csL, csG, nstep) tryCatch(
  as.numeric(attr(GET::global_envelope_test(curve_sets = list(csL, csG), type = "erl", alpha = TESTS$alpha, nstep = nstep), "p")),
  error = function(e) NA_real_)
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
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), p_nstep1 = erl_p(csL, csG, 1), usable = sum(okc))
}
null_rep <- function(i, nsim) {
  t0 <- unname(proc.time()["elapsed"])
  set.seed(SEEDS$master + SEEDS$offsets["empirical"] + 1000 + i)
  X <- to_poly(rpoispp(truth$lamW))
  r <- m1_test(X, nsim, bw_for_n(RECOVERY_ENV$expected_n), SEEDS$master + SEEDS$offsets["empirical"] + 2000 + i)
  c(p = r$p, p_nstep1 = r$p_nstep1, usable = r$usable, n = npoints(X), sec = unname(proc.time()["elapsed"]) - t0)
}
time_probe <- function(npat) {
  ts <- tf <- tc <- numeric(npat)
  bw <- bw_for_n(RECOVERY_ENV$expected_n)
  for (i in seq_len(npat)) {
    set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 3000 + i)
    ts[i] <- system.time(X <- to_poly(rpoispp(truth$lamW)))["elapsed"]
    tf[i] <- system.time(f <- fast_fit(X))["elapsed"]
    tc[i] <- system.time(cv <- fast_curves(X, f, bw))["elapsed"]
  }
  c(sim = median(ts), fit = median(tf), curves = median(tc), vertices = sum(vapply(W_km$bdry, function(b) length(b$x), 1L)), area = area_km2)
}
)---"

while (sink.number() > 0) sink()
source(REPO_CONFIG)
source(file.path(REPO_DIR, "R", "recovery_setup.R"))
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
worker_f <- file.path(REPO_DIR, "R", "recovery_worker.R"); func_f <- file.path(REPO_DIR, "R", "recovery_functions.R")
cache_exact <- file.path(REPO_DIR, "data_derived", "recovery_cache.rds")
cache_simpl <- file.path(REPO_DIR, "data_derived", "recovery_cache_simpl.rds")
stopifnot(file.exists(cache_exact))
writeLines(WORKER_TEXT, worker_f); writeLines(FUNCTIONS_TEXT, func_f)
out <- dir_of("10_simulation_recovery"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
M <- list()

mk_cluster <- function(n, cache_name) {
  cl <- parallel::makeCluster(n)
  parallel::clusterExport(cl, "REPO_CONFIG", envir = globalenv())
  parallel::clusterCall(cl, Sys.setenv, RECOVERY_CACHE = cache_name)
  msg <- parallel::clusterEvalQ(cl, tryCatch({
    source(REPO_CONFIG); source(file.path(REPO_DIR, "R", "recovery_worker.R")); source(file.path(REPO_DIR, "R", "recovery_functions.R"))
    paste0("ok (", CACHE_NAME, ", ", round(sum(gc()[, 2])), " MB)")
  }, error = function(e) paste("SETUP ERROR:", conditionMessage(e), "| R memory", round(sum(gc()[, 2])), "MB")))
  list(cl = cl, msg = unlist(msg))
}

sink(file.path(out, "window_speed_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== A. SIMPLIFIED WINDOW (100 m, st_simplify preserveTopology, as in script 08) ====\n")
  bnd <- st_transform(st_make_valid(st_read(file.path(RAW_DIR, man$path[man$layer == "Boundary"]), quiet = TRUE)), CRS_PROJ)
  win_sf <- st_union(st_geometry(bnd))
  simp <- st_union(st_make_valid(st_simplify(win_sf, preserveTopology = TRUE, dTolerance = 100)))
  Ws <- rescale(as.owin(simp), 1000, c("kilometre", "kilometres"))
  stopifnot(Ws$type == "polygonal")
  vx_e <- sum(vapply(W_km$bdry, function(b) length(b$x), 1L)); vx_s <- sum(vapply(Ws$bdry, function(b) length(b$x), 1L))
  a_e <- area.owin(W_km); a_s <- area.owin(Ws)
  sym <- as.numeric(st_area(st_sym_difference(win_sf, simp))) / 1e6
  ins <- inside.owin(X_obs_km$x, X_obs_km$y, Ws)
  cat(sprintf("vertices: exact %d | simplified %d\narea (km2): exact %.2f | simplified %.2f | change %.4f%% | symmetric difference %.1f km2\nobserved sites inside simplified window: %d of %d\n",
              vx_e, vx_s, a_e, a_s, 100 * (a_s - a_e) / a_e, sym, sum(ins), npoints(X_obs_km)))
  valid_simpl <- all(ins) && abs(a_s - a_e) / a_e < 1e-4
  cat("validation:", if (valid_simpl) "PASS (all sites inside, area change below 0.01%)" else "FAIL: simplified window not usable", "\n")
  write.csv(data.frame(item = c("vertices_exact", "vertices_simplified", "area_exact_km2", "area_simplified_km2", "symdiff_km2", "sites_inside", "valid"),
                       value = c(vx_e, vx_s, round(a_e, 2), round(a_s, 2), round(sym, 1), sum(ins), valid_simpl)), file.path(out, "window_variant_check.csv"), row.names = FALSE)
  if (valid_simpl) saveRDS(list(covs = covs, W_km = Ws, X_obs_km = X_obs_km[Ws], formula_E3 = formula_E3), cache_simpl)
  
  cat("\n==== B. PER-PATTERN TIMING: one worker per window variant ====\n")
  probe <- list()
  for (v in c("exact", if (valid_simpl) "simplified")) {
    cname <- if (v == "exact") "recovery_cache.rds" else "recovery_cache_simpl.rds"
    ts0 <- system.time(cl1 <- mk_cluster(1, cname))["elapsed"]
    cat(sprintf("%-10s worker setup %.1f s | message: %s\n", v, ts0, paste(cl1$msg, collapse = " | ")))
    r <- try(parallel::clusterCall(cl1$cl, time_probe, PROBE_N)[[1]], silent = TRUE)
    parallel::stopCluster(cl1$cl)
    if (inherits(r, "try-error")) { cat("  probe FAILED:", conditionMessage(attr(r, "condition")), "\n"); next }
    probe[[v]] <- r
    cat(sprintf("%-10s vertices %d | area %.1f | per pattern: simulate %.3f s + fit %.3f s + curves %.3f s = %.3f s\n", v, r["vertices"], r["area"], r["sim"], r["fit"], r["curves"], r["sim"] + r["fit"] + r["curves"]))
  }
  speed <- if (length(probe) == 2) (sum(probe$exact[c("sim", "fit", "curves")]) / sum(probe$simplified[c("sim", "fit", "curves")])) else NA
  cat(sprintf("speed-up of the simplified window per pattern: x%.2f\n", speed))
  use_variant <- if (valid_simpl && is.finite(speed) && speed >= 1.5) "simplified" else "exact"
  use_cache <- if (use_variant == "simplified") "recovery_cache_simpl.rds" else "recovery_cache.rds"
  cat("window used for the scaling and calibration steps below:", use_variant, "(measurement only; no decision is made here)\n")
  
  cat("\n==== C. REAL PARALLEL SPEED-UP (", PAR_TASKS, " tasks x ", PAR_NSIM, " simulations, same task list for every worker count) ====\n", sep = "")
  scal <- list()
  for (k in PAR_WORKERS) {
    cl <- mk_cluster(k, use_cache)
    tw <- system.time(rr <- parallel::parLapplyLB(cl$cl, seq_len(PAR_TASKS), function(i, nsim) tryCatch(null_rep(i, nsim), error = function(e) paste("ERROR:", conditionMessage(e))), nsim = PAR_NSIM))["elapsed"]
    parallel::stopCluster(cl$cl)
    errs <- vapply(rr, is.character, logical(1))
    scal[[length(scal) + 1]] <- data.frame(workers = k, wall_sec = round(tw, 1), errors = sum(errs))
    cat(sprintf("workers %d: wall %.1f s | errors %d\n", k, tw, sum(errs)))
    if (any(errs)) cat("  first error:", unique(unlist(rr[errs]))[1], "\n")
  }
  sc <- do.call(rbind, scal); sc$speedup_vs_1 <- round(sc$wall_sec[sc$workers == 1] / sc$wall_sec, 2)
  print(sc, row.names = FALSE); write.csv(sc, file.path(out, "parallel_scaling.csv"), row.names = FALSE)
  ok_sc <- sc[sc$errors == 0, ]; best_k <- ok_sc$workers[which.max(ok_sc$speedup_vs_1)]
  M$speedup <- max(ok_sc$speedup_vs_1)
  cat("best worker count:", best_k, "| real speed-up x", M$speedup, "\n")
  
  cat("\n==== D. NULL CALIBRATION, ", N_NULL_REPS, " replicates x ", NSIM_NULL, " simulations (", use_variant, " window, ", best_k, " workers) ====\n", sep = "")
  cc <- mk_cluster(best_k, use_cache)
  cat("setup messages:", paste(unique(cc$msg), collapse = " | "), "\n")
  tw <- system.time(pr <- parallel::parLapplyLB(cc$cl, SEED_OFFSET + seq_len(N_NULL_REPS), function(i, nsim) tryCatch(null_rep(i, nsim), error = function(e) paste("ERROR:", conditionMessage(e))), nsim = NSIM_NULL))["elapsed"]
  parallel::stopCluster(cc$cl)
  is_err <- vapply(pr, is.character, logical(1))
  if (any(is_err)) cat("worker errors:", paste(unique(unlist(pr[is_err])), collapse = " | "), "\n")
  if (any(!is_err)) {
    P <- do.call(rbind, pr[!is_err])
    cat(sprintf("completed %d of %d | wall %.1f min | mean worker time per replicate %.1f s | min usable simulations %d | mean n %.1f (min %d, max %d)\n",
                nrow(P), N_NULL_REPS, tw / 60, mean(P[, "sec"]), min(P[, "usable"]), mean(P[, "n"]), min(P[, "n"]), max(P[, "n"])))
    rep_stats <- function(pv, label) {
      pv <- pv[is.finite(pv)]
      if (!length(pv)) { cat(sprintf("%-8s no p-values\n", label)); return(invisible()) }
      ks <- suppressWarnings(ks.test(pv, "punif")$p.value)
      cat(sprintf("%-8s n %d | rejections at 0.05: %d (expect %.1f) | p <= 0.10: %d | median %.2f | share above 0.5: %.2f | KS vs uniform p = %.3f\n",
                  label, length(pv), sum(pv <= 0.05), 0.05 * length(pv), sum(pv <= 0.10), median(pv), mean(pv > 0.5), ks))
    }
    rep_stats(P[, "p"], "nstep 2"); rep_stats(P[, "p_nstep1"], "nstep 1")
    cat("sorted p-values, nstep 2:\n"); print(round(sort(P[, "p"]), 3))
    cat("sorted p-values, nstep 1:\n"); print(round(sort(P[, "p_nstep1"]), 3))
    write.csv(data.frame(rep = SEED_OFFSET + which(!is_err), P), file.path(out, "null_calibration_v2.csv"), row.names = FALSE)
  }
  
  cat("\n==== E. COST RE-ESTIMATE (hierarchy fractions ASSUMED; measured per-pattern time; real speed-up) ====\n")
  pr_use <- probe[[use_variant]]
  pipe <- if (!is.null(pr_use)) unname(sum(pr_use[c("sim", "fit", "curves")])) else 0.85
  kfit <- 0.34; ksim <- 0.34; rmh_ps <- 0.61; rmh_ov <- 10.8; geyer_fit_s <- 1.5      # kppm, rmh from script 19; Geyer fit assumed
  cost <- function(nsim, reps) {
    m1 <- nsim * pipe; m3 <- nsim * (ksim + kfit + pipe); m4 <- nsim * (rmh_ps + pipe) + rmh_ov
    fit_all <- pipe + kfit + 8 * geyer_fit_s
    pois <- fit_all + m1 + 0.1 * m3; lgcp <- fit_all + m1 + 0.9 * m3 + 0.1 * m4; gey <- fit_all + m1 + 0.9 * m3 + 0.7 * m4
    per_set <- 6 * pois + 7 * lgcp + 6 * gey
    c(hours_per_replicate_index = per_set / 3600, study_single_core_h = reps * per_set / 3600)
  }
  est <- rbind(data.frame(workflow = "frozen counts (200 reps)", nsim = TESTS$nsim_recovery, reps = TESTS$recovery_reps, t(cost(TESTS$nsim_recovery, TESTS$recovery_reps))),
               data.frame(workflow = "199 sims, 100 reps", nsim = TESTS$nsim_recovery, reps = 100, t(cost(TESTS$nsim_recovery, 100))),
               data.frame(workflow = "floors (100 reps, 99 sims)", nsim = TESTS$nsim_recovery_floor, reps = TESTS$recovery_reps_floor, t(cost(TESTS$nsim_recovery_floor, TESTS$recovery_reps_floor))))
  est$days_on_workers <- round(est$study_single_core_h / M$speedup / 24, 1); est[, 4:5] <- round(est[, 4:5], 1)
  print(est, row.names = FALSE)
  cat(sprintf("window variant %s | per-pattern pipeline %.3f s | measured real speed-up x%.2f on %d workers\n", use_variant, pipe, M$speedup, best_k))
  write.csv(est, file.path(out, "cost_estimate_v4.csv"), row.names = FALSE)
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_21.txt"))
cat("\nDONE. Outputs in:", out, "\n")