# =============================================================================
# 19_fast_pipeline.R   (v2: scoping fixed, workspace cleared, fixed-dummy variant, debug prints)
# Stage: 10_simulation_recovery | cost reduction of the per-pattern fit
# Purpose: (A) decompose ppm fit time; (B) hand-built Berman-Turner Poisson fit via
#   glm.fit on spatstat's own quadrature scheme (per-pattern vs prebuilt dummy grid),
#   validated against ppm; (C) full-size M1 test timing; (D) null calibration of the
#   M1 test (parallel); (E) Geyer truth multiplier calibration; (F) LGCP truth mean n;
#   (G) cost re-estimate; (H) log implementation clarifications. SYNTHETIC patterns only.
# Requires R/recovery_setup.R (script 17).
# Outputs (10_simulation_recovery/outputs/): fast_pipeline_report.txt, fast_pipeline.csv,
#   null_calibration.csv, geyer_multipliers.csv, lgcp_mean_n.csv, cost_estimate_v2.csv,
#   sessionInfo_19.txt ; docs/protocol_deviation_log.csv (rows appended)
# Not executed by assistant. Run top to bottom in RStudio (25-45 min).
# WARNING: first line clears the objects in your R workspace (not files).
# =============================================================================

rm(list = ls(envir = globalenv()), envir = globalenv())

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
N_NULL_REPS <- 24L; NSIM_NULL <- 99L; NSIM_FULL <- 199L
RMH_NREP <- 5e5; RMH_BATCH <- 8L

while (sink.number() > 0) sink()
source(REPO_CONFIG)
source(file.path(REPO_DIR, "R", "recovery_setup.R"))
out <- dir_of("10_simulation_recovery"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)

covsW <- lapply(covs, function(im) im[W_km, drop = FALSE])
area_km2 <- area.owin(W_km)
bw_nom <- (0.15 / sqrt(RECOVERY_ENV$expected_n / area_km2)) / sqrt(5)
L_ARGS    <- list(renormalise = TRUE, normpower = 1)
G_ARGS_NS <- list(kernel = "epanechnikov", zerocor = "convolution", divisor = "r", renormalise = TRUE, normpower = 2)
maxdiff <- function(a, b) suppressWarnings(max(abs(a - b), na.rm = TRUE))
truth <- make_truth()
DUMMY_FIXED <- NULL

fast_fit <- function(Xp, debug = FALSE) {
  dbg <- function(label, x) if (debug) cat(sprintf("  [fast_fit] %-12s class %-14s length %d\n", label, class(x)[1], length(x)))
  Q <- if (is.null(DUMMY_FIXED)) quadscheme(Xp, eps = QUAD_KM) else quadscheme(Xp, dummy = DUMMY_FIXED, method = "grid")
  dbg("quadscheme", Q)
  w <- w.quad(Q); isd <- is.data(Q); xq <- x.quad(Q); yq <- y.quad(Q)
  dbg("weights", w); dbg("is.data", isd)
  Z <- vapply(covs, function(im) lookup.im(im, xq, yq, naok = TRUE), numeric(length(xq)))
  dbg("covariates", Z)
  ok <- stats::complete.cases(Z) & w > 0
  X <- cbind(`(Intercept)` = 1, Z)
  gfit <- suppressWarnings(stats::glm.fit(X[ok, , drop = FALSE], (isd / w)[ok], weights = w[ok], family = poisson(), intercept = FALSE))
  dbg("glm.fit", gfit)
  beta <- gfit$coefficients
  Xd <- X[isd & ok, , drop = FALSE]
  if (nrow(Xd) != npoints(Xp)) stop("data points dropped in fast_fit: ", nrow(Xd), " of ", npoints(Xp))
  list(beta = beta, lam_data = as.numeric(exp(Xd %*% beta)), n_quad = sum(ok))
}
lam_im_of <- function(beta) exp(beta[1] + Reduce(`+`, Map(function(a, b) a * b, covsW, beta[-1])))
fast_curves <- function(Xp, ftb, bw = bw_nom) {
  Lh <- do.call(Linhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "isotropic"), L_ARGS))
  gh <- do.call(pcfinhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "Ripley", bw = bw), G_ARGS_NS))
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(rG, 6), round(gh$r, 6))])
}
null_rep <- function(i, nsim) {
  set.seed(SEEDS$master + SEEDS$offsets["empirical"] + 1000 + i)
  t0 <- unname(proc.time()["elapsed"])
  X <- to_poly(rpoispp(truth$lamW)); f0 <- fast_fit(X); ob <- fast_curves(X, f0)
  lam_im <- lam_im_of(f0$beta)
  sL <- matrix(NA_real_, length(rL), nsim); sG <- matrix(NA_real_, length(rG), nsim)
  for (k in seq_len(nsim)) { Xs <- to_poly(rpoispp(lam_im)); fs <- fast_fit(Xs); cs <- fast_curves(Xs, fs); sL[, k] <- cs$L; sG[, k] <- cs$g }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  res <- GET::global_envelope_test(curve_sets = list(
    GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE])),
    GET::create_curve_set(list(r = rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))), type = "erl", alpha = TESTS$alpha)
  c(p = as.numeric(attr(res, "p")), usable = sum(okc), sec = unname(proc.time()["elapsed"]) - t0)
}

tab <- list(); add <- function(step, variant, sec = NA, value = NA) tab[[length(tab) + 1]] <<- data.frame(step = step, variant = variant, sec = sec, value = value)
M <- list()

sink(file.path(out, "fast_pipeline_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== 0. SETUP ====\n")
  cat(sprintf("window area %.1f km2 | nominal bandwidth: half-width h = %.3f km, sd bw = %.3f km | quadrature %g km\n", area_km2, bw_nom * sqrt(5), bw_nom, QUAD_KM))
  set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 1)
  Xp <- to_poly(rpoispp(truth$lamW)); cat("synthetic pattern n =", npoints(Xp), "\n")
  
  sec("A decomposition", {
    cat("\n==== A. WHERE THE ppm FIT TIME GOES ====\n")
    tq <- system.time(Q <- quadscheme(Xp, eps = QUAD_KM))["elapsed"]
    tpQ <- system.time(f1 <- ppm(Q, trend = formula_E3, covariates = covs))["elapsed"]
    tp  <- system.time(f0 <- ppm(Xp, trend = formula_E3, covariates = covs, quad.args = list(eps = QUAD_KM)))["elapsed"]
    cat(sprintf("quadscheme build %.2f s | ppm(prebuilt Q) %.2f s | ppm(pattern) %.2f s | quad points %d\n", tq, tpQ, tp, n.quad(Q)))
    add("ppm decomposition", "quadscheme", tq); add("ppm decomposition", "ppm(prebuilt Q)", tpQ); add("ppm decomposition", "ppm(pattern)", tp)
  })
  
  cat("\n==== B. HAND-BUILT BERMAN-TURNER FIT vs ppm ====\n")
  cat("step-by-step debug run (per-pattern quadscheme):\n")
  dbg_res <- tryCatch(fast_fit(Xp, debug = TRUE), error = function(e) { cat("fast_fit FAILED at:", paste(deparse(conditionCall(e)), collapse = " "), "\n  message:", conditionMessage(e), "\n"); NULL })
  if (is.null(dbg_res)) stop("fast_fit failed; see the message above.")
  
  DUMMY_FIXED <- NULL
  t1 <- numeric(5); for (i in 1:5) t1[i] <- system.time(ftb <- fast_fit(Xp))["elapsed"]
  cat(sprintf("variant 1, quadscheme per pattern: median %.3f s over 5 | quad points used %d\n", median(t1), ftb$n_quad))
  t_dum <- system.time(D_fixed <- default.dummy(Xp, eps = QUAD_KM))["elapsed"]
  cat(sprintf("building the fixed dummy grid once: %.2f s (%d points)\n", t_dum, npoints(D_fixed)))
  DUMMY_FIXED <- D_fixed
  t2 <- numeric(5); v2_ok <- TRUE
  v2 <- tryCatch({ for (i in 1:5) t2[i] <- system.time(ftb2 <- fast_fit(Xp))["elapsed"]; TRUE }, error = function(e) { cat("variant 2 FAILED:", conditionMessage(e), "\n"); FALSE })
  if (v2) cat(sprintf("variant 2, prebuilt dummy grid: median %.3f s over 5 | quad points used %d | max |coef diff| vs variant 1 %.2e\n", median(t2), ftb2$n_quad, maxdiff(ftb$beta, ftb2$beta)))
  DUMMY_FIXED <- NULL
  
  check_variant <- function(label, dummy) {
    DUMMY_FIXED <<- dummy; worst <- 0; worst_l <- 0
    for (i in 1:3) {
      set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 50 + i)
      Xi <- to_poly(rpoispp(truth$lamW)); fi <- fast_fit(Xi)
      pi_ <- ppm(Xi, trend = formula_E3, covariates = covs, quad.args = list(eps = QUAD_KM))
      cf <- coef(pi_)[names(fi$beta)]
      dc <- maxdiff(fi$beta, cf); dl <- max(abs(fi$lam_data / as.numeric(fitted(pi_, dataonly = TRUE)) - 1))
      cat(sprintf("  %s | pattern %d (n=%d): max |coef diff| %.2e | max relative |lambda diff| %.2e\n", label, i, npoints(Xi), dc, dl))
      worst <- max(worst, dc); worst_l <- max(worst_l, dl)
    }
    DUMMY_FIXED <<- NULL
    c(coef = worst, lam = worst_l)
  }
  w1 <- check_variant("variant 1", NULL)
  w2 <- if (v2) check_variant("variant 2", D_fixed) else c(coef = Inf, lam = Inf)
  use_fixed <- v2 && is.finite(w2["coef"]) && w2["coef"] < 1e-4 && median(t2) < median(t1)
  valid1 <- is.finite(w1["coef"]) && w1["coef"] < 1e-4
  cat(sprintf("\nvalidation: variant 1 worst coef diff %.2e (%s) | variant 2 worst coef diff %.2e (%s)\n", w1["coef"], if (valid1) "valid" else "NOT valid", w2["coef"], if (w2["coef"] < 1e-4) "valid" else "NOT valid"))
  if (!valid1 && !use_fixed) stop("Hand-built fit does not reproduce ppm. Not adopted; later sections skipped.")
  DUMMY_FIXED <- if (use_fixed) D_fixed else NULL
  cat("ADOPTED:", if (use_fixed) "variant 2 (prebuilt dummy grid)" else "variant 1 (quadscheme per pattern)", "\n")
  M$fast_fit <- if (use_fixed) median(t2) else median(t1)
  add("fast_fit", "variant 1 median", median(t1)); if (v2) add("fast_fit", "variant 2 median", median(t2)); add("fast_fit", "adopted", M$fast_fit)
  
  ftb <- fast_fit(Xp)
  lam_im <- lam_im_of(ftb$beta)
  tr_v <- numeric(20); for (i in 1:20) tr_v[i] <- system.time(Xs <- to_poly(rpoispp(lam_im)))["elapsed"]
  M$sim_pois <- median(tr_v)
  cat(sprintf("rpoispp + window restriction: %.3f s per pattern\n", M$sim_pois))
  tc_v <- numeric(5); for (i in 1:5) tc_v[i] <- system.time(cc <- fast_curves(Xp, ftb))["elapsed"]
  M$curves <- median(tc_v)
  cat(sprintf("L_inhom + g_inhom (plain intensity, fixed bw): %.2f s | NA: L %d g %d\n", M$curves, sum(is.na(cc$L)), sum(is.na(cc$g))))
  M$pipe <- M$fast_fit + M$curves
  cat(sprintf("full per-pattern pipeline: %.2f s (was 3.39 s plain, 30.18 s leave-one-out, 56.7 s original)\n", M$pipe))
  add("per-pattern pipeline", "fast_fit + curves", M$pipe); add("rpoispp", "per pattern", M$sim_pois)
  
  sec("C full-size M1 test", {
    cat("\n==== C. ONE FULL-SIZE M1 TEST (", NSIM_FULL, " bootstrap patterns, serial) ====\n", sep = "")
    r <- run("M1 test, Poisson truth, 199 simulations", null_rep(0, NSIM_FULL))
    if (r$ok) { cat(sprintf("  p = %.3f | usable simulations %d | worker seconds %.1f\n", r$res["p"], r$res["usable"], r$res["sec"])); M$m1_test <- unname(r$res["sec"]); add("M1 test", "199 sims serial", unname(r$res["sec"]), unname(r$res["p"])) }
  })
  
  sec("D null calibration", {
    cat("\n==== D. NULL CALIBRATION OF THE M1 TEST (", N_NULL_REPS, " Poisson-truth replicates x ", NSIM_NULL, " simulations, parallel) ====\n", sep = "")
    ncl <- min(4, parallel::detectCores(logical = FALSE)); cl <- parallel::makeCluster(ncl)
    exports <- c("covs", "covsW", "W_km", "QUAD_KM", "RGRID", "rL", "rG", "col_of", "L_ARGS", "G_ARGS_NS", "bw_nom", "fast_fit", "lam_im_of",
                 "fast_curves", "to_poly", "truth", "SEEDS", "TESTS", "null_rep", "RECOVERY_ENV", "DUMMY_FIXED", "maxdiff")
    res <- tryCatch({
      parallel::clusterExport(cl, exports, envir = globalenv())
      parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(spatstat.geom); library(spatstat.explore); library(spatstat.model); library(spatstat.random); library(GET) }))
      tw <- system.time(pr <- parallel::parLapply(cl, seq_len(N_NULL_REPS), function(i) tryCatch(null_rep(i, NSIM_NULL), error = function(e) c(p = NA, usable = NA, sec = NA))))["elapsed"]
      list(pr = do.call(rbind, pr), wall = unname(tw))
    }, error = function(e) { cat("parallel FAILED:", conditionMessage(e), "\n"); NULL }, finally = parallel::stopCluster(cl))
    if (!is.null(res)) {
      P <- res$pr; ok <- !is.na(P[, "p"]); pv <- P[ok, "p"]
      cat(sprintf("replicates completed %d of %d | wall time %.1f min | mean worker time per replicate %.1f s\n", sum(ok), N_NULL_REPS, res$wall / 60, mean(P[ok, "sec"])))
      cat(sprintf("p-values: rejections at 0.05: %d of %d (expect about %.1f) | p <= 0.10: %d | median p %.2f | min p %.3f\n",
                  sum(pv <= 0.05), sum(ok), 0.05 * sum(ok), sum(pv <= 0.10), median(pv), min(pv)))
      print(round(sort(pv), 3))
      M$par_speedup <- sum(P[ok, "sec"]) / res$wall; M$rep_sec <- mean(P[ok, "sec"])
      cat(sprintf("parallel efficiency: %.2f workers' worth of speed-up on %d workers\n", M$par_speedup, ncl))
      write.csv(data.frame(rep = which(ok), p = pv, sec = P[ok, "sec"]), file.path(out, "null_calibration.csv"), row.names = FALSE)
      add("null calibration", "rejections at 0.05", NA, sum(pv <= 0.05)); add("null calibration", "wall minutes", res$wall / 60, NA)
    }
  })
  
  sec("E Geyer calibration", {
    cat("\n==== E. GEYER TRUTH MULTIPLIER CALIBRATION (n ~ 197, converged chains) ====\n")
    mk <- function(g, r, m) rmhmodel(cif = "geyer", par = list(beta = 1, gamma = g, r = r, sat = 3), w = W_km, trend = truth$lamW * m)
    set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 700)
    tg1 <- system.time(P1 <- rmh(mk(2.5, 3, 0.165), start = list(n.start = 197), control = list(nrep = RMH_NREP, expand = 1), verbose = FALSE))["elapsed"]
    tg8 <- system.time(P8 <- rmh(mk(2.5, 3, 0.165), start = list(n.start = 197), control = list(nrep = RMH_NREP, expand = 1), nsim = RMH_BATCH, verbose = FALSE))["elapsed"]
    M$rmh_per_sim <- (tg8 - tg1) / (RMH_BATCH - 1); M$rmh_overhead <- tg1 - M$rmh_per_sim
    cat(sprintf("rmh at %g steps: nsim=1 %.1f s | nsim=%d %.1f s | per extra simulation %.2f s | overhead %.1f s\n", RMH_NREP, tg1, RMH_BATCH, tg8, M$rmh_per_sim, M$rmh_overhead))
    add("rmh", "per extra simulation", M$rmh_per_sim); add("rmh", "overhead", M$rmh_overhead)
    cells <- data.frame(cell = c("G1", "G2", "G3", "G4"), gamma = c(1.5, 1.5, 2.5, 2.5), r = c(3, 7, 3, 7), start = c(0.30, 0.28, 0.165, 0.15))
    res_g <- list()
    for (k in seq_len(nrow(cells))) {
      m <- cells$start[k]; hist <- NULL; ok <- FALSE
      for (it in 1:6) {
        set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 800 + k * 10 + it)
        PP <- rmh(mk(cells$gamma[k], cells$r[k], m), start = list(n.start = 197), control = list(nrep = RMH_NREP, expand = 1), nsim = RMH_BATCH, verbose = FALSE)
        ns <- sapply(PP, function(z) npoints(z[W_km])); nb <- mean(ns)
        hist <- rbind(hist, data.frame(cell = cells$cell[k], iter = it, mult = round(m, 4), mean_n = round(nb, 1), sd_n = round(sd(ns), 1)))
        if (abs(nb / 197 - 1) < 0.04) { ok <- TRUE; break }
        m <- m * (197 / nb)^(1 / 1.5)
      }
      print(hist, row.names = FALSE); res_g[[k]] <- cbind(hist[nrow(hist), ], converged = ok)
    }
    gm <- do.call(rbind, res_g); write.csv(gm, file.path(out, "geyer_multipliers.csv"), row.names = FALSE)
    cat("Calibrated multipliers (final iteration):\n"); print(gm, row.names = FALSE)
  })
  
  sec("F LGCP mean n", {
    cat("\n==== F. LGCP TRUTH MEAN n AND kppm / SIMULATION TIMING ====\n")
    lg <- data.frame(cell = c("L1", "L2", "L3", "L4"), var = c(0.5, 0.5, 1.5, 1.5), scale = c(3, 8, 3, 8)); rows <- list()
    for (k in 1:4) {
      set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 900 + k)
      PL <- rLGCP("exp", mu = log(truth$lamW) - lg$var[k] / 2, param = list(var = lg$var[k], scale = lg$scale[k]), saveLambda = FALSE, dimyx = c(504, 362), nsim = 20)
      ns <- sapply(PL, function(z) npoints(z[W_km])); rows[[k]] <- data.frame(lg[k, ], mean_n = round(mean(ns), 1), sd_n = round(sd(ns), 1), min_n = min(ns), max_n = max(ns))
    }
    lm_ <- do.call(rbind, rows); print(lm_, row.names = FALSE); write.csv(lm_, file.path(out, "lgcp_mean_n.csv"), row.names = FALSE)
    XL <- to_poly(PL[[1]])
    tkf <- system.time(kp <- kppm(XL, trend = formula_E3, clusters = "LGCP", covariates = covs, method = "clik2"))["elapsed"]; M$kfit <- unname(tkf)
    tks <- tryCatch(system.time(S <- simulate(kp, nsim = 20, verbose = FALSE, dimyx = c(504, 362)))["elapsed"] / 20,
                    error = function(e) { cat("simulate(kppm, dimyx) failed:", conditionMessage(e), "\n"); system.time(S <- simulate(kp, nsim = 20, verbose = FALSE))["elapsed"] / 20 })
    M$ksim <- unname(tks)
    cat(sprintf("kppm fit %.2f s | simulate from fitted kppm %.2f s per pattern (batch of 20)\n", M$kfit, M$ksim)); add("kppm", "fit", M$kfit); add("kppm", "simulate per pattern", M$ksim)
  })
  
  sec("G cost", {
    cat("\n==== G. COST RE-ESTIMATE (single core; hierarchy fractions are ASSUMED) ====\n")
    g <- function(x, d = NA_real_) if (is.null(x) || !is.finite(x)) d else unname(x)
    tpipe <- g(M$pipe); tsp <- g(M$sim_pois); tk <- g(M$kfit); tks <- g(M$ksim); trm <- g(M$rmh_per_sim, 0.6); tov <- g(M$rmh_overhead, 10); tgrid <- 8 * 3.2
    cost <- function(nsim, reps) {
      m1 <- nsim * (tsp + tpipe); m3 <- nsim * (tks + tk + tpipe); m4 <- nsim * (trm + tpipe) + tov
      fit_all <- tpipe + tk + tgrid
      pois <- fit_all + m1 + 0.1 * m3; lgcp <- fit_all + m1 + 0.9 * m3 + 0.1 * m4; gey <- fit_all + m1 + 0.9 * m3 + 0.7 * m4
      per_set <- 6 * pois + 7 * lgcp + 6 * gey
      c(per_replicate_index_h = per_set / 3600, study_single_core_h = reps * per_set / 3600)
    }
    sp <- g(M$par_speedup, 2.3)
    est <- rbind(data.frame(workflow = "frozen counts", nsim = TESTS$nsim_recovery, reps = TESTS$recovery_reps, t(cost(TESTS$nsim_recovery, TESTS$recovery_reps))),
                 data.frame(workflow = "floors", nsim = TESTS$nsim_recovery_floor, reps = TESTS$recovery_reps_floor, t(cost(TESTS$nsim_recovery_floor, TESTS$recovery_reps_floor))))
    est$days_on_4_cores <- round(est$study_single_core_h / sp / 24, 1); est[, 4:5] <- round(est[, 4:5], 1)
    print(est, row.names = FALSE); cat(sprintf("measured parallel speed-up used: x%.2f\n", sp))
    cat("Measured components (s):", paste(names(M), round(unlist(M), 3), sep = "=", collapse = " | "), "\n")
    write.csv(est, file.path(out, "cost_estimate_v2.csv"), row.names = FALSE)
  })
  
  sec("H log", {
    dl_f <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv")
    dl <- if (file.exists(dl_f)) read.csv(dl_f, stringsAsFactors = FALSE) else data.frame(date = character(), item = character(), reason = character(), effect_on_inference = character(), approved_by = character())
    new <- data.frame(date = format(Sys.Date()), approved_by = "author-delegated (assistant)", stringsAsFactors = FALSE,
                      item = c("CLARIFICATION pcf bandwidth: Stoyan value is kernel half-width; bandwidth held fixed within each test",
                               "CLARIFICATION test statistics use plain (not leave-one-out) fitted intensity",
                               "CLARIFICATION rmh simulation in one batch call; >= 5e5 steps; truth multipliers calibrated at converged chains"),
                      reason = c("Script 18: spatstat default bw = h/sqrt(5) where h = 0.15/sqrt(lambda); freeze values 1.56 km (A) and 2.22 km (C) are half-widths. Bandwidth fixed at the nominal-n value so observed and simulated patterns share one estimator.",
                                 "Script 18: leave-one-out intensity costs 30 s per pattern; plain 0 s; applied identically to observed and simulated patterns inside a parametric bootstrap with refit. Leave-one-out kept as sensitivity on the observed pattern.",
                                 "Script 18: rmh has about 10 s fixed overhead and 0.6 s per extra batched simulation; mean n stable between 5e5 and 2e6 steps, not at 1e3-1e5 for some replicates."),
                      effect_on_inference = c("None on frozen quantities (lower bounds 2.0 km A, 2.5 km C unchanged). Removes bandwidth variability among simulated patterns.",
                                              "None on frozen quantities. Test validity rests on the null calibration reported by script 19.",
                                              "None on frozen quantities. Gibbs simulations are better converged."))
    new <- new[!new$item %in% dl$item, ]
    if (nrow(new)) { write.csv(rbind(dl[, names(new)], new), dl_f, row.names = FALSE); cat("\nAppended", nrow(new), "clarification row(s) to docs/protocol_deviation_log.csv\n") }
  })
  
  write.csv(do.call(rbind, tab), file.path(out, "fast_pipeline.csv"), row.names = FALSE)
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_19.txt"))
cat("\nDONE. Outputs in:", out, "\n")