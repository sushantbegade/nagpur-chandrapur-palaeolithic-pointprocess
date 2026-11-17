# =====================================================================
# 44_empirical_m2_offsets.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW)
# Purpose : M2 = environmental Poisson (E3) with a FIXED log-E offset: O2 (kappa 1.1192)
#           and O3 (kappa 2.5366), frozen layers at 500 m. Dataset A, simplified window,
#           nominal coordinates. Per pattern: refit E3 + offset (glm.fit offset), intensity
#           for L_inhom / g_inhom = fitted value INCLUDING the offset; simulate from the
#           fitted intensity incl. offset (unconditional rpoispp); joint ERL test of
#           L_inhom (1-10 km) and g_inhom (2-10 km), nsim 999, nstep 2. Same summary
#           functions as scripts 40-43.
# Outputs : outputs/44_report.txt, 44_m2_result.csv, 44_m2_curves.csv, 44_m2_coef.csv,
#           44_sessionInfo.txt
# Notes   : RUN ONCE. Stops if 44_m2_result.csv exists. Frozen rule (freeze_v1.md):
#           sufficiency under O2 stops the hierarchy; sufficiency only under O3 is labelled
#           observation-sensitive and the hierarchy continues. Offsets are ASSUMED strengths,
#           not estimates; O3 contains a pattern-informed core term (not blind). M2 has not
#           been calibrated in the recovery study (nuisance cells N1, N5 fitted without offset).
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
f_res <- file.path(OUT_DIR, "44_m2_result.csv")
if (file.exists(f_res)) stop("44_m2_result.csv exists: this analysis runs once. Stop; do not rerun.")
f_obs <- file.path(REPO_DIR, "data_derived", "recovery_obs_cache.rds"); stopifnot(file.exists(f_obs))
OBSL <- readRDS(f_obs); stopifnot(all(c("logE_O2", "logE_O3") %in% names(OBSL)))

NSIM <- TESTS$nsim_empirical; N_A <- N_ANALYTICAL
SEED_BASE <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 30000
bw <- bw_for_n(N_A); rGv <- rG_for(N_A)
PAIR_R <- c(1, 2, 3, 5)

WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
pairs_at <- function(P, r) sum(closepaircounts(P, r)) / 2
fast_fit_off <- function(Xp, offim) {
  Q <- make_quad(Xp); w <- w.quad(Q); isd <- is.data(Q); xq <- x.quad(Q); yq <- y.quad(Q)
  Z <- vapply(covs, function(im) lookup.im(im, xq, yq, naok = TRUE), numeric(length(xq)))
  o <- lookup.im(offim, xq, yq, naok = TRUE)
  ok <- stats::complete.cases(Z) & is.finite(o) & w > 0
  X <- cbind(`(Intercept)` = 1, Z)
  gfit <- suppressWarnings(stats::glm.fit(X[ok, , drop = FALSE], (isd / w)[ok], weights = w[ok],
                                          family = poisson(), offset = o[ok], intercept = FALSE))
  beta <- gfit$coefficients
  Xd <- X[isd & ok, , drop = FALSE]; od <- o[isd & ok]
  if (nrow(Xd) != npoints(Xp)) stop("data points dropped in fast_fit_off: ", nrow(Xd), " of ", npoints(Xp))
  list(beta = beta, lam_data = as.numeric(exp(Xd %*% beta + od)))
}
emp_m2 <- function(Xobs, nsim, bw, rGv, seed, offim) {
  set.seed(seed)
  offW <- offim[W_km, drop = FALSE]
  f0 <- fast_fit_off(Xobs, offim); ob <- fast_curves_r(Xobs, f0, bw, rGv)
  lam_im <- lam_im_of(f0$beta) * exp(offW)
  sL <- matrix(NA_real_, length(rL), nsim); sG <- matrix(NA_real_, length(rGv), nsim)
  nn <- rep(NA_real_, nsim); pr <- matrix(NA_real_, length(PAIR_R), nsim)
  for (k in seq_len(nsim)) {
    res <- tryCatch({
      Xs <- to_poly(rpoispp(lam_im)); fs <- fast_fit_off(Xs, offim); cs <- fast_curves_r(Xs, fs, bw, rGv)
      list(cs = cs, n = npoints(Xs), pr = vapply(PAIR_R, function(r) pairs_at(Xs, r), numeric(1)))
    }, error = function(e) NULL)
    if (!is.null(res)) { sL[, k] <- res$cs$L; sG[, k] <- res$cs$g; nn[k] <- res$n; pr[, k] <- res$pr }
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rGv, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc), ob = ob, f0 = f0, int = integral(lam_im),
       sL = sL[, okc, drop = FALSE], sG = sG[, okc, drop = FALSE], n = nn[okc], pairs = pr[, okc, drop = FALSE])
}
env_df <- function(scn, curve, r, obs, S) data.frame(
  scenario = scn, curve = curve, r_km = r, obs = obs, sim_mean = rowMeans(S),
  lo = apply(S, 1, quantile, 0.025), hi = apply(S, 1, quantile, 0.975),
  z_raw = (obs - rowMeans(S)) / apply(S, 1, sd))

REPORT <- file.path(OUT_DIR, "44_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("44_empirical_m2_offsets | ", format(Sys.time()), "\n", sep = "")
    cat(sprintf("frozen settings: E3 + fixed offset, dataset A, simplified window | nsim %d | ERL nstep 2 | alpha %.2f\n", NSIM, TESTS$alpha))
    SC <- list(O2 = list(name = "logE_O2", kappa = OBS$kappa[["O2"]], seed = SEED_BASE + 2),
               O3 = list(name = "logE_O3", kappa = OBS$kappa[["O3"]], seed = SEED_BASE + 3))
    OUTR <- data.frame(); CVR <- data.frame(); CFR <- data.frame(); PV <- c()
    for (sn in names(SC)) {
      offim <- OBSL[[SC[[sn]]$name]]
      cat(sprintf("\n================ %s (kappa %.4f) ================\n", sn, SC[[sn]]$kappa))
      cat(sprintf("offset range %.4f-%.4f (mean %.3f)\n", min(offim$v, na.rm = TRUE), max(offim$v, na.rm = TRUE), mean(offim$v, na.rm = TRUE)))
      f0 <- fast_fit_off(X_obs_km, offim)
      fitO <- tryCatch(ppm(make_quad(X_obs_km), trend = update(formula_E3, . ~ . + offset(logEoff)),
                           covariates = c(covs, list(logEoff = offim))), error = function(e) NULL)
      if (is.null(fitO)) cat("ppm cross-check FAILED (continuing)\n") else {
        b <- coef(fitO); V <- tryCatch(vcov(fitO), error = function(e) NULL)
        se <- if (is.null(V)) rep(NA_real_, length(b)) else sqrt(diag(V))
        cat(sprintf("max |coef| difference fast_fit_off vs ppm(offset): %.2e\n", max(abs(f0$beta - b[names(f0$beta)]))))
        tab <- data.frame(scenario = sn, term = names(b), coef = unname(b), se = unname(se), z = unname(b / se))
        print(tab[, -1], digits = 3, row.names = FALSE); CFR <- rbind(CFR, tab)
      }
      t0 <- unname(proc.time()["elapsed"])
      R <- emp_m2(X_obs_km, NSIM, bw, rGv, SC[[sn]]$seed, offim)
      cat(sprintf("\nsimulations done in %.0f s | usable %d of %d | integral of fitted intensity %.2f (n = %d) | simulated n mean %.1f sd %.1f\n",
                  unname(proc.time()["elapsed"]) - t0, R$usable, NSIM, R$int, N_A, mean(R$n), sd(R$n)))
      cat(sprintf("ERL combined global envelope p (nstep 2): %.4f | deficiency at 0.05: %s\n", R$p,
                  if (is.na(R$p)) "NA" else if (R$p < TESTS$alpha) "YES" else "NO"))
      zL <- zsum(R$ob$L, R$sL); zG <- zsum(R$ob$g, R$sG)
      zr <- abs(c(zL$raw, zG$raw)); zr[!is.finite(zr)] <- NA; rr <- c(rL, rGv); ia <- which.max(zr)
      cat(sprintf("raw max |z| %.2f at r = %.1f km on the %s curve\n", zr[ia], rr[ia], if (ia <= length(rL)) "L_inhom" else "g_inhom"))
      CV <- rbind(env_df(sn, "L_inhom_minus_r", rL, R$ob$L, R$sL), env_df(sn, "g_inhom", rGv, R$ob$g, R$sG))
      print(CV[round(CV$r_km, 6) %in% c(2, 3, 5, 7, 10), -1], digits = 3, row.names = FALSE)
      po <- vapply(PAIR_R, function(r) pairs_at(X_obs_km, r), numeric(1))
      cat("pairs within 1, 2, 3, 5 km: observed", paste(po, collapse = " "), "| simulated mean", paste(round(rowMeans(R$pairs), 1), collapse = " "), "\n")
      CVR <- rbind(CVR, CV); PV[sn] <- R$p
      OUTR <- rbind(OUTR, data.frame(scenario = sn, kappa = SC[[sn]]$kappa, n_obs = N_A, nsim = NSIM, usable = R$usable,
                                     p_erl_nstep2 = R$p, zraw_max = zr[ia], zraw_r_km = rr[ia], seed = SC[[sn]]$seed))
    }
    cat("\n== Frozen decision rule ==\n")
    act <- if (is.na(PV["O2"]) || is.na(PV["O3"])) "UNDETERMINED (NA p)" else if (PV["O2"] >= TESTS$alpha)
      "STOP: sufficiency under O2 (assumed strength; not proof of sufficiency)" else if (PV["O3"] >= TESTS$alpha)
        "CONTINUE to M3: sufficiency only under O3 = observation-sensitive" else "CONTINUE to M3: deficient under both O2 and O3"
    cat(sprintf("M1 (O1) p 0.0310 | M2-O2 p %.4f | M2-O3 p %.4f -> %s\n", PV["O2"], PV["O3"], act))
    OUTR$decision <- act; OUTR$date <- format(Sys.time())
    write.csv(OUTR, f_res, row.names = FALSE)
    write.csv(CVR, file.path(OUT_DIR, "44_m2_curves.csv"), row.names = FALSE)
    write.csv(CFR, file.path(OUT_DIR, "44_m2_coef.csv"), row.names = FALSE)
    cat("\n== Warning tally (main session) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "44_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")