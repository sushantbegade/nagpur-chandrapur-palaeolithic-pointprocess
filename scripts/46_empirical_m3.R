# =====================================================================
# 46_empirical_m3.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW)
# Purpose : M3 = environmental LGCP (E3 trend, O1) on the archaeological pattern:
#           kppm (exponential, clik2), conditional simulation (rLGCP eps 0.5, saveLambda,
#           in chunks of 100, then rpoint(n | Lambda)), nsim 999, joint ERL test
#           (L_inhom 1-10 km, g_inhom 2-10 km, nstep 2), floored z-max, parametric
#           bootstrap refits (B = 199) at the point estimate for sigma2 and alpha,
#           curves and pair counts vs the M3 simulations. Same functions as recovery.
# Outputs : outputs/46_report.txt, 46_m3_result.csv, 46_m3_curves.csv, 46_m3_boot.csv,
#           46_sessionInfo.txt
# Notes   : RUN ONCE. Stops if 46_m3_result.csv exists. Dataset A, simplified window,
#           nominal coordinates. Reading rules (set before running): p < 0.05 = deficiency
#           -> M4 (script 47, descriptive under the Gate 6 pivot). p in [0.04, 0.06] ->
#           precision rule: rerun at nsim 10000 (script 46b). p >= 0.06 -> hierarchy stops
#           (latent field adequate; NOT distinguishable from interaction, Gate 6).
#           sigma2_hat and alpha_hat are estimator outputs (Stage 1: sigma2 inflated, alpha
#           low-biased, alpha > 10 km in some fits), not physical parameters.
#           Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))
source(file.path(REPO_DIR, "R", "recovery_tests.R"))
source(file.path(REPO_DIR, "R", "recovery_cells.R"))

OUT_DIR <- dir_of("11_generative_gof", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
f_res <- file.path(OUT_DIR, "46_m3_result.csv")
if (file.exists(f_res)) stop("46_m3_result.csv exists: this analysis runs once. Stop; do not rerun.")

NSIM <- TESTS$nsim_empirical; N_A <- N_ANALYTICAL
CHK <- 100; B_BOOT <- 199
SEED_E <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 30000 + 4
bw <- bw_for_n(N_A); rGv <- rG_for(N_A)
PAIR_R <- c(1, 2, 3, 5)

WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
pairs_at <- function(P, r) sum(closepaircounts(P, r)) / 2
fit_k <- function(X) tryCatch(kppm(X, trend = formula_E3, clusters = "LGCP", model = "exponential",
                                   covariates = covs, method = "clik2"), error = function(e) NULL)

REPORT <- file.path(OUT_DIR, "46_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("46_empirical_m3 | ", format(Sys.time()), "\n", sep = "")
    cat(sprintf("frozen settings: E3 trend (O1), LGCP exponential, clik2 | dataset A, simplified window | nsim %d conditional on n | ERL nstep 2 | alpha %.2f | seed %d\n",
                NSIM, TESTS$alpha, SEED_E))
    cat(sprintf("sites %d | NA-covariate sites %d (must be 0)\n", npoints(X_obs_km), sum(na_cov_pts(X_obs_km))))
    stopifnot(npoints(X_obs_km) == N_A, sum(na_cov_pts(X_obs_km)) == 0)
    
    cat("\n== 1. kppm fit ==\n")
    fit <- fit_k(X_obs_km)
    if (is.null(fit)) stop("kppm fit failed on the observed pattern; paste output")
    s2 <- unname(fit$par["sigma2"]); al <- unname(fit$par["alpha"])
    cat(sprintf("sigma2_hat %.4f | alpha_hat %.3f km (decision range upper bound %.0f km)\n", s2, al, R_MAX_DECISION_KM))
    b <- coef(fit); print(round(b, 3))
    cat("Stage 1 reference (estimator, L4 truth n about 190): sigma2_hat median 1.82-2.02, alpha_hat median 5.3-6.4 km\n")
    
    cat("\n== 2. Conditional simulation (chunks of ", CHK, ") and curves ==\n", sep = "")
    eta <- trend_im(coef(fit)); n <- npoints(X_obs_km)
    f0 <- fast_fit(X_obs_km); ob <- fast_curves_r(X_obs_km, f0, bw, rGv)
    sL <- matrix(NA_real_, length(rL), NSIM); sG <- matrix(NA_real_, length(rGv), NSIM)
    pr <- matrix(NA_real_, length(PAIR_R), NSIM); sims_first <- list(); ndrop <- 0
    t0 <- unname(proc.time()["elapsed"]); k0 <- 0
    nch <- ceiling(NSIM / CHK)
    for (ch in seq_len(nch)) {
      m <- min(CHK, NSIM - k0)
      set.seed(SEED_E + 100 * ch)
      S <- rLGCP(model = "exp", mu = eta - s2 / 2, param = list(var = s2, scale = al), win = W_km,
                 eps = EPS_LGCP, saveLambda = TRUE, nsim = m, drop = TRUE)
      Sl <- as_list(S)
      for (j in seq_len(m)) {
        k <- k0 + j
        Xr <- to_poly(rpoint(n, f = attr(Sl[[j]], "Lambda"), win = W_km, forcewin = TRUE))
        Xs <- clean_pat(Xr); ndrop <- ndrop + npoints(Xr) - npoints(Xs)
        res <- tryCatch({ fs <- fast_fit(Xs); cs <- fast_curves_r(Xs, fs, bw, rGv)
        list(cs = cs, pr = vapply(PAIR_R, function(r) pairs_at(Xs, r), numeric(1))) }, error = function(e) NULL)
        if (!is.null(res)) { sL[, k] <- res$cs$L; sG[, k] <- res$cs$g; pr[, k] <- res$pr }
        if (k <= B_BOOT) sims_first[[k]] <- Xs
      }
      k0 <- k0 + m
      cat(sprintf("  chunk %2d of %d done [%.0f s elapsed]\n", ch, nch, unname(proc.time()["elapsed"]) - t0))
    }
    okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
    csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
    csG <- GET::create_curve_set(list(r = rGv, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
    p3 <- erl_p(csL, csG, 2)
    cat(sprintf("\nusable simulations %d of %d | synthetic points dropped by cleaning %d\n", sum(okc), NSIM, ndrop))
    cat(sprintf("ERL combined global envelope p (nstep 2): %.4f | Monte Carlo SE %.4f | deficiency at 0.05: %s\n", p3,
                sqrt(p3 * (1 - p3) / sum(okc)), if (is.na(p3)) "NA" else if (p3 < TESTS$alpha) "YES" else "NO"))
    zL <- zsum(ob$L, sL[, okc, drop = FALSE]); zG <- zsum(ob$g, sG[, okc, drop = FALSE])
    zfl <- max(abs(c(zL$flo, zG$flo)), na.rm = TRUE)
    zr <- abs(c(zL$raw, zG$raw)); zr[!is.finite(zr)] <- NA; rr <- c(rL, rGv); ia <- which.max(zr)
    cat(sprintf("max |z| (floored sd) %.2f | raw max |z| %.2f at r = %.1f km on the %s curve (M1 reference: 6.58 floored)\n",
                zfl, zr[ia], rr[ia], if (ia <= length(rL)) "L_inhom" else "g_inhom"))
    
    mk <- function(curve, r, o, S) data.frame(curve = curve, r_km = r, obs = o, sim_mean = rowMeans(S),
                                              lo = apply(S, 1, quantile, 0.025), hi = apply(S, 1, quantile, 0.975), z_raw = (o - rowMeans(S)) / apply(S, 1, sd))
    CV <- rbind(mk("L_inhom_minus_r", rL, ob$L, sL[, okc, drop = FALSE]), mk("g_inhom", rGv, ob$g, sG[, okc, drop = FALSE]))
    write.csv(CV, file.path(OUT_DIR, "46_m3_curves.csv"), row.names = FALSE)
    cat("\n-- Observed vs M3 simulations (pointwise 95% band, display only) --\n")
    print(CV[round(CV$r_km, 6) %in% c(1, 2, 3, 4, 5, 7, 10), ], digits = 3, row.names = FALSE)
    po <- vapply(PAIR_R, function(r) pairs_at(X_obs_km, r), numeric(1))
    cat("\n-- Pair counts: observed vs M3 simulations --\n")
    print(data.frame(r_km = PAIR_R, observed = po, sim_mean = rowMeans(pr[, okc, drop = FALSE]),
                     sim_lo = apply(pr[, okc, drop = FALSE], 1, quantile, 0.025), sim_hi = apply(pr[, okc, drop = FALSE], 1, quantile, 0.975),
                     share_sim_ge_obs = vapply(seq_along(PAIR_R), function(i) mean(pr[i, okc] >= po[i]), numeric(1))), digits = 3, row.names = FALSE)
    
    cat("\n== 3. Parametric bootstrap refits at the point estimate (B = ", B_BOOT, ") ==\n", sep = "")
    bs <- t(vapply(seq_len(min(B_BOOT, length(sims_first))), function(k) {
      f <- fit_k(sims_first[[k]]); if (is.null(f)) c(NA_real_, NA_real_) else unname(c(f$par["sigma2"], f$par["alpha"])) }, numeric(2)))
    colnames(bs) <- c("sigma2", "alpha"); write.csv(bs, file.path(OUT_DIR, "46_m3_boot.csv"), row.names = FALSE)
    q <- function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)
    cat(sprintf("refits ok %d of %d\n", sum(!is.na(bs[, 1])), nrow(bs)))
    cat(sprintf("sigma2*: median %.3f (95%% %.3f-%.3f) | estimate %.3f | median/estimate %.2f | share < 0.01: %.2f\n",
                q(bs[, 1])[2], q(bs[, 1])[1], q(bs[, 1])[3], s2, q(bs[, 1])[2] / s2, mean(bs[, 1] < 0.01, na.rm = TRUE)))
    cat(sprintf("alpha*:  median %.2f (95%% %.2f-%.2f) | estimate %.2f | median/estimate %.2f | share > 10 km: %.2f\n",
                q(bs[, 2])[2], q(bs[, 2])[1], q(bs[, 2])[3], al, q(bs[, 2])[2] / al, mean(bs[, 2] > R_MAX_DECISION_KM, na.rm = TRUE)))
    
    act <- if (is.na(p3)) "UNDETERMINED (NA p)" else if (p3 >= 0.04 && p3 <= 0.06)
      "PRECISION RULE TRIGGERED: rerun M3 at nsim 10000 (script 46b)" else if (p3 < TESTS$alpha)
        "DEFICIENCY under M3: continue to M4 (script 47, descriptive under the Gate 6 pivot)" else
          "NO deficiency under M3: hierarchy stops (latent field adequate; not distinguishable from interaction)"
    cat("\nDECISION (frozen threshold, precision rule):", act, "\n")
    write.csv(data.frame(test = "M3", n_obs = N_A, nsim = NSIM, usable = sum(okc), p_erl_nstep2 = p3, zmax_floored = zfl,
                         zraw_max = zr[ia], zraw_r_km = rr[ia], sigma2_hat = s2, alpha_hat = al, seed = SEED_E,
                         decision = act, date = format(Sys.time())), f_res, row.names = FALSE)
    cat("\n== Warning tally (main session) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "46_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")