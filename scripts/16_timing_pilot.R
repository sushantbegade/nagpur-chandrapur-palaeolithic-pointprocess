# =============================================================================
# 16_timing_pilot.R   (v2: polygon-window enforcement, rLGCP param, section isolation)
# Stage: 10_simulation_recovery | pre-recovery feasibility
# Purpose: measure run times of the frozen workflow on SYNTHETIC patterns and
#   estimate recovery-study cost. NO archaeological model is fitted.
# Requires: freeze v1 FROZEN, data_derived/ppp_NA_window.rds, env rasters (500 m),
#   packages spatstat.geom/explore/model/random, GET, terra, sf.
# Outputs (10_simulation_recovery/outputs/): timing_pilot_report.txt,
#   timing_pilot.csv, cost_estimate.csv, sessionInfo_16.txt
# Not executed by assistant. Run top to bottom in RStudio. Long run (30-90 min).
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
NSIM_PILOT   <- 19L
NREP_RMH     <- 100000L
N_GEYER_SIMS <- 3L
GEYER_CAL    <- 0.5            # trend multiplier for the Geyer truth (pilot only)
QUAD_KM      <- 1
RGRID        <- seq(0, 10, by = 0.5)

while (sink.number() > 0) sink()
source(REPO_CONFIG)
params_f <- file.path(REPO_DIR, "R", "freeze_v1_params.R"); stopifnot(file.exists(params_f)); source(params_f)
if (!grepl("^FROZEN", FREEZE_STATUS)) stop("Freeze v1 is not FROZEN. Aborting.")
need <- c("spatstat.geom", "spatstat.explore", "spatstat.model", "spatstat.random", "GET", "terra", "sf")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: install.packages(c(", paste(sprintf('"%s"', miss), collapse = ", "), "))")
suppressPackageStartupMessages({ library(sf); library(spatstat.geom); library(spatstat.explore)
  library(spatstat.model); library(spatstat.random); library(GET) })

out <- dir_of("10_simulation_recovery"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
env_dir <- file.path(REPO_DIR, "data_derived", "env", "500m")
obj <- readRDS(file.path(REPO_DIR, "data_derived", "ppp_NA_window.rds"))
X_obs_km <- obj$X_km; W_km <- Window(X_obs_km)
stopifnot(W_km$type == "polygonal")
rl <- function(n) terra::rast(file.path(env_dir, paste0(n, ".tif")))
for (n in c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5", "lith5", "geomorph3"))
  stopifnot(file.exists(file.path(env_dir, paste0(n, ".tif"))))

to_im <- function(r) {
  m <- terra::as.matrix(r, wide = TRUE); m <- m[nrow(m):1, , drop = FALSE]
  im(m, xcol = terra::xFromCol(r, seq_len(ncol(r))) / 1000, yrow = rev(terra::yFromRow(r, seq_len(nrow(r)))) / 1000)
}
cont <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5")
lay <- lapply(cont, rl); names(lay) <- cont
l5 <- rl("lith5"); g3 <- rl("geomorph3")
lay$lith_basalt <- terra::ifel(l5 == 1, 1, 0); lay$lith_sandstone <- terra::ifel(l5 == 3, 1, 0)
lay$lith_carbonate <- terra::ifel(l5 == 4, 1, 0); lay$lith_other <- terra::ifel(l5 == 5, 1, 0)
lay$geom_dissected <- terra::ifel(g3 == 2, 1, 0); lay$geom_alluvial <- terra::ifel(g3 == 3, 1, 0)
covs <- lapply(lay, to_im)
formula_E3 <- as.formula(paste("~", paste(names(covs), collapse = " + ")))

# every pattern must live in the exact polygon window (isotropic correction needs it)
to_poly <- function(P, label = "") {
  w0 <- Window(P)$type
  Q <- P[W_km]
  w1 <- Window(Q)$type
  cat(sprintf("  window %s: %s -> %s | n %d -> %d\n", label, w0, w1, npoints(P), npoints(Q)))
  if (w1 != "polygonal") stop("pattern window is not polygonal after restriction")
  Q
}
T <- list()
run <- function(label, expr) {
  cat(sprintf("%-58s", label)); flush.console()
  t <- system.time(r <- try(expr, silent = TRUE))["elapsed"]; ok <- !inherits(r, "try-error")
  cat(sprintf("%8.1f s  %s\n", t, if (ok) "ok" else paste("FAILED:", conditionMessage(attr(r, "condition")))))
  list(res = if (ok) r else NULL, sec = unname(t), ok = ok)
}
sec <- function(label, expr) tryCatch(expr, error = function(e) cat(sprintf("\n[SECTION FAILED: %s] %s\n", label, conditionMessage(e))))
probe <- function(f) paste(tryCatch(names(formals(f)), error = function(e) "n/a"), collapse = ", ")
col_of <- function(f, pref) { for (nm in pref) if (nm %in% names(f)) return(f[[nm]]); stop("no column among: ", paste(names(f), collapse = ",")) }
rL <- RGRID[RGRID >= R_MIN_KM["A"] & RGRID <= R_MAX_DECISION_KM]
rG <- RGRID[RGRID >= GINHOM_LOWER_KM["A"] & RGRID <= R_MAX_DECISION_KM]
curves <- function(Xp) {
  if (Window(Xp)$type != "polygonal") stop("pattern window is not polygonal")
  fit <- ppm(Xp, trend = formula_E3, covariates = covs, quad.args = list(eps = QUAD_KM))
  Lh <- Linhom(Xp, lambda = fit, r = RGRID, correction = "isotropic", update = FALSE)
  gh <- pcfinhom(Xp, lambda = fit, r = RGRID, correction = "Ripley", update = FALSE)
  Lv <- col_of(Lh, c("iso", "trans", "border")) - Lh$r; gv <- col_of(gh, c("iso", "trans"))
  list(L = Lv[match(round(rL, 6), round(Lh$r, 6))], g = gv[match(round(rG, 6), round(gh$r, 6))],
       bw = attr(gh, "bw"), fit = fit)
}

sink(file.path(out, "timing_pilot_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== A. VERSIONS ====\n")
  for (p in c("spatstat.geom", "spatstat.explore", "spatstat.model", "spatstat.random", "GET", "terra", "sf"))
    cat(sprintf("%-18s %s\n", p, as.character(packageVersion(p))))
  cat("cores (physical):", parallel::detectCores(logical = FALSE), "\n")
  
  cat("\n==== B. SYNTHETIC TRUTH (freeze recovery environment) ====\n")
  set.seed(SEEDS$master + RECOVERY_ENV$seed_offset)
  beta <- setNames(rnorm(length(covs), 0, RECOVERY_ENV$coef_sd), names(covs))
  eta <- Reduce(`+`, Map(function(a, b) a * b, covs, beta))
  lam0 <- exp(eta)[W_km, drop = FALSE]
  b0 <- log(RECOVERY_ENV$expected_n / integral(lam0))
  lamW <- lam0 * exp(b0)
  cat(sprintf("intercept %.3f | integral of intensity over window: %.1f (target %d) | truth image window type: %s\n",
              b0, integral(lamW), RECOVERY_ENV$expected_n, Window(lamW)$type))
  
  sec("C Poisson pipeline", {
    cat("\n==== C. POISSON PIPELINE ====\n")
    set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 1)
    Xp <- to_poly(rpoispp(lamW), "Poisson truth")
    r1 <- run("trend-only fit + L_inhom + g_inhom (one pattern)", curves(Xp)); T$stats <- r1$sec
    if (r1$ok) {
      cat(sprintf("  pcf bandwidth used (km): %.3f | NA in L: %d | NA in g: %d\n", r1$res$bw, sum(is.na(r1$res$L)), sum(is.na(r1$res$g))))
      cf <- coef(r1$res$fit)[names(beta)]; cat("  coefficient recovery, correlation with truth:", round(cor(cf, beta), 3), "\n")
      cat("  fitted-model simulation window check:\n")
      r2 <- run("simulate one pattern from fitted Poisson model", to_poly(simulate(r1$res$fit, nsim = 1, drop = TRUE, verbose = FALSE), "sim from fit")); T$sim_pois <- r2$sec
      cat("\nParametric bootstrap with refit,", NSIM_PILOT, "patterns:\n")
      simL <- matrix(NA_real_, length(rL), NSIM_PILOT); simG <- matrix(NA_real_, length(rG), NSIM_PILOT)
      t_sim <- numeric(NSIM_PILOT); t_cur <- numeric(NSIM_PILOT)
      for (i in seq_len(NSIM_PILOT)) {
        set.seed(SEEDS$master + SEEDS$offsets["bootstrap"] + i)
        t_sim[i] <- system.time(Xs <- try({ Z <- simulate(r1$res$fit, nsim = 1, drop = TRUE, verbose = FALSE); Z[W_km] }, silent = TRUE))["elapsed"]
        if (inherits(Xs, "try-error")) next
        t_cur[i] <- system.time(cv <- try(curves(Xs), silent = TRUE))["elapsed"]
        if (!inherits(cv, "try-error")) { simL[, i] <- cv$L; simG[, i] <- cv$g }
        if (i %% 5 == 0) cat("  done", i, "of", NSIM_PILOT, "\n")
      }
      T$boot_sim <- median(t_sim); T$boot_curves <- median(t_cur)
      ok <- colSums(is.na(simL)) == 0 & colSums(is.na(simG)) == 0
      cat(sprintf("median per pattern: simulate %.2f s | refit + curves %.2f s | usable patterns: %d of %d\n", T$boot_sim, T$boot_curves, sum(ok), NSIM_PILOT))
      if (sum(ok) >= 9) {
        r3 <- run("combined ERL test via GET (L_inhom + g_inhom)", {
          csL <- create_curve_set(list(r = rL, obs = r1$res$L, sim_m = simL[, ok, drop = FALSE]))
          csG <- create_curve_set(list(r = rG, obs = r1$res$g, sim_m = simG[, ok, drop = FALSE]))
          global_envelope_test(curve_sets = list(csL, csG), type = "erl", alpha = TESTS$alpha) })
        T$erl <- r3$sec
        if (r3$ok) cat(sprintf("  combined ERL p-value (pilot, %d sims): %s  [Poisson truth: expect not small]\n", sum(ok), format(attr(r3$res, "p"))))
      } else cat("Too few usable bootstrap patterns for the ERL test.\n")
    }
  })
  
  sec("D LGCP", {
    cat("\n==== D. LGCP (cell L4: sigma2 1.5, scale 8 km) ====\n")
    var_t <- 1.5; scale_t <- 8
    set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 2)
    mu_im <- log(lamW) - var_t / 2
    rg <- run("simulate LGCP truth (param route, dimyx)", rLGCP("exp", mu = mu_im, param = list(var = var_t, scale = scale_t), saveLambda = TRUE, dimyx = c(504, 362)))
    if (!rg$ok) rg <- run("simulate LGCP truth (var/scale via ..., dimyx)", rLGCP("exp", mu = mu_im, var = var_t, scale = scale_t, saveLambda = TRUE, dimyx = c(504, 362)))
    T$sim_lgcp_truth <- rg$sec
    if (rg$ok) {
      Lam <- attr(rg$res, "Lambda")
      if (!is.null(Lam)) { Z <- try(log(Lam) - mu_im, silent = TRUE)
      if (!inherits(Z, "try-error")) cat(sprintf("  realised variance of log-intensity field minus mean: %.3f (nominal %.1f; finite window reduces it)\n", var(as.vector(as.matrix(Z)), na.rm = TRUE), var_t)) }
      XL <- to_poly(rg$res, "LGCP truth")
    } else { cat("Using the Poisson pattern for kppm timing.\n"); XL <- Xp }
    kp <- run("kppm LGCP fit, method clik2", kppm(XL, trend = formula_E3, clusters = "LGCP", covariates = covs, method = "clik2"))
    if (!kp$ok) kp <- run("kppm LGCP fit, method mincon (fallback)", kppm(XL, trend = formula_E3, clusters = "LGCP", covariates = covs, method = "mincon"))
    T$kppm_fit <- kp$sec
    if (kp$ok) {
      cat("  fitted cluster parameters (truth: sigma2 1.5, scale 8 km; single-pattern sanity only):\n")
      print(tryCatch(parameters(kp$res), error = function(e) "n/a"))
      rs <- run("simulate one pattern from fitted kppm", to_poly(simulate(kp$res, nsim = 1, drop = TRUE, verbose = FALSE), "sim from kppm")); T$sim_kppm <- rs$sec
    }
  })
  
  sec("E Geyer", {
    cat("\n==== E. GEYER (cell G3: gamma 2.5, r 3 km, sat 3) ====\n")
    mod <- try(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = 2.5, r = 3, sat = 3), w = W_km, trend = lamW * GEYER_CAL), silent = TRUE)
    if (inherits(mod, "try-error")) { cat("rmhmodel FAILED:", conditionMessage(attr(mod, "condition")), "\n") } else {
      XG <- list(); tg <- numeric(0)
      for (i in seq_len(N_GEYER_SIMS)) {
        set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 10 + i)
        r <- run(sprintf("rmh Geyer truth %d (nrep %d)", i, NREP_RMH), to_poly(rmh(mod, start = list(n.start = 100), control = list(nrep = NREP_RMH, expand = 1), verbose = FALSE), "Geyer truth"))
        if (r$ok) { XG[[length(XG) + 1]] <- r$res; tg <- c(tg, r$sec) }
      }
      if (length(XG)) {
        ns <- sapply(XG, npoints); T$rmh <- median(tg)
        cat(sprintf("Geyer truth n: %s | mean %.1f | trend multiplier giving n = 197 (rough, pilot only): %.2f\n", paste(ns, collapse = ", "), mean(ns), GEYER_CAL * 197 / mean(ns)))
        Xg <- XG[[1]]
        combos <- expand.grid(r = GEYER_R_KM, sat = GEYER_SAT)
        fits <- list(); tab <- list()
        for (rb in c(10, 4)) for (k in seq_len(nrow(combos))) {
          f <- run(sprintf("Geyer fit r=%g sat=%g rbord=%g", combos$r[k], combos$sat[k], rb),
                   ppm(Xg, trend = formula_E3, interaction = Geyer(r = combos$r[k], sat = combos$sat[k]), covariates = covs, rbord = rb, quad.args = list(eps = QUAD_KM)))
          gh <- NA_real_
          if (f$ok) { cf <- coef(f$res); ii <- grep("Interaction", names(cf)); if (length(ii)) gh <- as.numeric(exp(cf[ii[1]])) }
          tab[[length(tab) + 1]] <- data.frame(rbord = rb, r = combos$r[k], sat = combos$sat[k], sec = round(f$sec, 1), ok = f$ok,
                                               logPL = if (f$ok) as.numeric(logLik(f$res)) else NA, gamma_hat = round(gh, 3))
          if (f$ok) fits[[paste(rb, combos$r[k], combos$sat[k])]] <- f$res
        }
        tb <- do.call(rbind, tab); cat("\nProfile table (synthetic Geyer pattern; truth gamma 2.5, r 3, sat 3):\n"); print(tb, row.names = FALSE)
        T$geyer_fit <- median(tb$sec[tb$ok]); T$geyer_grid <- sum(tb$sec[tb$ok & tb$rbord == 10])
        best <- which.max(ifelse(tb$ok & tb$rbord == 10, tb$logPL, -Inf))
        if (length(best) && is.finite(tb$logPL[best])) {
          key <- paste(10, tb$r[best], tb$sat[best]); cat("Best (rbord 10):", key, "\n")
          rsim <- run("simulate one pattern from the best Geyer fit", to_poly(simulate(fits[[key]], nsim = 1, drop = TRUE, verbose = FALSE, control = list(nrep = NREP_RMH, expand = 1)), "sim from Geyer fit"))
          T$sim_geyer_fit <- rsim$sec
        }
      }
    }
  })
  
  sec("F border cost", {
    cat("\n==== F. BORDER COST: observed sites beyond distance d from the window edge (geometry only) ====\n")
    bd <- bdist.points(X_obs_km)
    print(data.frame(rbord_km = c(2, 4, 5, 7, 10), sites_retained = sapply(c(2, 4, 5, 7, 10), function(d) sum(bd > d)),
                     share_of_197 = round(sapply(c(2, 4, 5, 7, 10), function(d) mean(bd > d)), 3)), row.names = FALSE)
    cat("Eroded window area share (script 07): 5 km 0.763 | 10 km 0.570.\n")
  })
  
  sec("G cost estimate", {
    cat("\n==== G. COST ESTIMATE (single core, hours) ====\n")
    g <- function(x, d = NA_real_) if (is.null(x) || !is.finite(x)) d else x
    ts <- g(T$stats); tsp <- g(T$sim_pois); tk <- g(T$kppm_fit); tsk <- g(T$sim_kppm); tgrid <- g(T$geyer_grid); trmh <- g(T$rmh)
    cost <- function(nsim, reps, refit) {
      m1 <- nsim * (tsp + ts)
      m3 <- nsim * (tsk + ts + if (refit) tk else 0)
      m4 <- nsim * (trmh + ts + if (refit) tgrid else 0)
      per_rep <- ts + tk + tgrid + m1 + m3 + m4
      c(per_replicate_min = per_rep / 60, per_cell_h = reps * per_rep / 3600, study_19_cells_h = 19 * reps * per_rep / 3600)
    }
    est <- rbind(
      data.frame(workflow = "frozen counts, model refit in each simulation", nsim = TESTS$nsim_recovery, reps = TESTS$recovery_reps, t(cost(TESTS$nsim_recovery, TESTS$recovery_reps, TRUE))),
      data.frame(workflow = "floors, model refit in each simulation", nsim = TESTS$nsim_recovery_floor, reps = TESTS$recovery_reps_floor, t(cost(TESTS$nsim_recovery_floor, TESTS$recovery_reps_floor, TRUE))),
      data.frame(workflow = "frozen counts, no model refit (trend-only refit in curves)", nsim = TESTS$nsim_recovery, reps = TESTS$recovery_reps, t(cost(TESTS$nsim_recovery, TESTS$recovery_reps, FALSE))),
      data.frame(workflow = "floors, no model refit", nsim = TESTS$nsim_recovery_floor, reps = TESTS$recovery_reps_floor, t(cost(TESTS$nsim_recovery_floor, TESTS$recovery_reps_floor, FALSE))))
    est[, -(1:3)] <- round(est[, -(1:3)], 2)
    options(width = 220); print(est, row.names = FALSE)
    cat("Measured components (s):", paste(names(T), round(unlist(T), 2), sep = "=", collapse = " | "), "\n")
    cat("Divide hours by usable cores for a parallel estimate. No decision is made by this script.\n")
    write.csv(est, file.path(out, "cost_estimate.csv"), row.names = FALSE)
    write.csv(data.frame(component = names(T), seconds = round(unlist(T), 3)), file.path(out, "timing_pilot.csv"), row.names = FALSE)
  })
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_16.txt"))
cat("\nDONE. Outputs in:", out, "\n")