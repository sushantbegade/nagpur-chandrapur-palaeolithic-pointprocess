# =====================================================================
# 35_na_cov_fix_m3_rerun_and_log.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : (A) diagnose NA-covariate points in synthetic patterns (script 34,
#               M3 rep 4: "data points dropped in fast_fit: 214 of 215"):
#               NA-covariate pixels inside the window; NA-covariate points in
#               the observed pattern, rLGCP, rpoint(n | Lambda), rpoispp, rmh;
#           (B) rerun the script-34 M3 pilot (same 4 L4 truth patterns, same
#               seeds) with NA-covariate points removed from every synthetic
#               pattern before fitting (clean_pat), counts logged;
#           (C) append clarification rows to docs/protocol_deviation_log.csv.
# Outputs : outputs/35_report.txt, 35_na_cov_diag.csv, 35_m3_pilot_rerun.csv,
#           35_sessionInfo.txt, protocol_deviation_log_before_35.csv (backup);
#           docs/protocol_deviation_log.csv (appended)
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           The observed archaeological pattern is only checked for NA covariates.
#           Hypothesis "NA-covariate pixels inside the polygon" is an INFERENCE.
#           Whole run under withCallingHandlers (warnings tallied).
#           Restart R before running. Main session only (no workers).
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))

OUT_DIR <- dir_of("10_simulation_recovery", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEED7      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 23000   # as script 34 (same 4 truth patterns)
SEED8      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 25000
N_PILOT    <- 4
NSIM_PILOT <- 99
NREP_RMH   <- 5e5
LGCP_VAR   <- 1.5
LGCP_SCALE <- 8
EPS_LGCP   <- 0.5
N_DIAG     <- 40
N_DIAG_GY  <- 8
P34        <- c(0.390, 0.590, 0.500, NA)   # script-34 ERL p (VERIFIED paste); rep 4 failed

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULT_G3 <- gm$mult[gm$cell == "G3"]; stopifnot(length(MULT_G3) == 1)
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")

# ---- helpers ----
WARN <- new.env(); WARN$msgs <- character(0)
note_w  <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
fill_im <- function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
trend_im <- function(b) fill_im(Reduce(`+`, Map(function(a, nmc) a * unname(b[nmc]), covs, names(covs))) + unname(b["(Intercept)"]))
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
na_report <- function(label, L) {
  tot  <- sum(vapply(L, npoints, numeric(1)))
  bads <- lapply(L, na_cov_pts)
  nb   <- vapply(bads, sum, numeric(1))
  bd   <- unlist(Map(function(P, b) if (any(b)) bdist.points(P)[b] else numeric(0), L, bads))
  cat(sprintf("  %-36s patterns %2d | points %5d | NA-cov points %3d (%.2f%%) | patterns affected %2d | boundary dist km of NA points: %s\n",
              label, length(L), tot, sum(nb), 100 * sum(nb) / tot, sum(nb > 0),
              if (length(bd)) sprintf("min %.2f median %.2f max %.2f", min(bd), median(bd), max(bd)) else "none"))
  data.frame(source = label, patterns = length(L), points = tot, na_points = sum(nb),
             patterns_affected = sum(nb > 0), bd_median_km = if (length(bd)) median(bd) else NA_real_)
}

# generic ERL test with NA-covariate cleaning of observed and every simulated pattern
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
  dobs <- npoints(Xobs) - npoints(Xc)
  fit <- tryCatch(kppm(Xc, trend = formula_E3, clusters = "LGCP", model = "exponential",
                       covariates = covs, method = "clik2"), error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA_real_, usable = 0L, drop_sim = NA_integer_, sims_dropped = NA_integer_,
                                sigma2 = NA_real_, alpha = NA_real_, drop_obs = dobs))
  s2 <- unname(fit$par["sigma2"]); al <- unname(fit$par["alpha"])
  eta <- trend_im(coef(fit)); n <- npoints(Xc)
  S <- rLGCP(model = "exp", mu = eta - s2 / 2, param = list(var = s2, scale = al), win = W_km,
             eps = EPS_LGCP, saveLambda = TRUE, nsim = nsim, drop = TRUE)
  sims <- lapply(as_list(S), function(P) rpoint(n, f = attr(P, "Lambda"), win = W_km, forcewin = TRUE))
  r <- curves_test(Xc, sims, bw)
  c(r, list(sigma2 = s2, alpha = al, drop_obs = dobs))
}

REPORT <- file.path(OUT_DIR, "35_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("35_na_cov_fix_m3_rerun_and_log | ", format(Sys.time()), "\n", sep = "")
    bw <- bw_for_n(RECOVERY_ENV$expected_n)
    mu0 <- Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta)) + truth$b0
    mu_fill <- fill_im(mu0)
    
    ## ---------------- A. diagnosis ----------------
    cat("\n== A1. NA-covariate pixels inside the simplified window ==\n")
    same_grid <- length(unique(vapply(covs, function(im) paste(dim(im$v), collapse = "x"), ""))) == 1
    cat("All covariate images on one grid:", same_grid, "\n")
    if (same_grid) {
      pix <- covs[[1]]; nr <- nrow(pix$v); nc <- ncol(pix$v)
      V   <- Reduce(`|`, lapply(covs, function(im) is.na(as.vector(im$v))))
      gx  <- rep(pix$xcol, each = nr); gy <- rep(pix$yrow, times = nc)
      ins <- inside.owin(gx, gy, W_km); sel <- ins & V
      cat(sprintf("Pixels with centre inside W: %d | of these NA in at least one covariate: %d (area %.1f km2, %.3f%% of W)\n",
                  sum(ins), sum(sel), sum(sel) * pix$xstep * pix$ystep, 100 * sum(sel) * pix$xstep * pix$ystep / area.owin(W_km)))
      if (sum(sel) > 0) {
        bdp <- bdist.points(ppp(gx[sel], gy[sel], window = W_km, check = FALSE))
        cat(sprintf("Boundary distance (km) of those pixel centres: min %.2f median %.2f max %.2f\n", min(bdp), median(bdp), max(bdp)))
      }
    }
    cat("\n== A2. NA-covariate points by source ==\n")
    n_obs_na <- sum(na_cov_pts(X_obs_km))
    cat(sprintf("  OBSERVED archaeological pattern: %d of %d sites have NA covariates (must be 0)\n", n_obs_na, npoints(X_obs_km)))
    DG <- list()
    set.seed(SEED8 + 1)
    LG <- rLGCP(model = "exp", mu = mu_fill - LGCP_VAR / 2, param = list(var = LGCP_VAR, scale = LGCP_SCALE),
                win = W_km, eps = EPS_LGCP, saveLambda = TRUE, nsim = N_DIAG, drop = TRUE)
    LGl <- as_list(LG)
    DG[[1]] <- na_report("rLGCP L4 (eps 0.5), to_poly", lapply(LGl, to_poly))
    set.seed(SEED8 + 2)
    RPl <- lapply(LGl, function(P) rpoint(197, f = attr(P, "Lambda"), win = W_km, forcewin = TRUE))
    DG[[2]] <- na_report("rpoint(197 | Lambda), forcewin", RPl)
    set.seed(SEED8 + 3)
    POl <- lapply(seq_len(N_DIAG), function(i) to_poly(rpoispp(truth$lamW)))
    DG[[3]] <- na_report("rpoispp(truth$lamW) control", POl)
    GYl <- as_list(sim_geyer(MULT_G3, 2.5, 3, 3, 197, N_DIAG_GY, SEED8 + 4))
    DG[[4]] <- na_report("rmh Geyer G3 truth", lapply(GYl, to_poly))
    write.csv(do.call(rbind, DG), file.path(OUT_DIR, "35_na_cov_diag.csv"), row.names = FALSE)
    
    ## ---------------- B. M3 pilot rerun ----------------
    cat("\n== B. M3 pilot rerun (same 4 L4 truth patterns and seeds as script 34; NA-covariate cleaning on) ==\n")
    set.seed(SEED7 + 1)
    TL <- rLGCP(model = "exp", mu = mu_fill - LGCP_VAR / 2, param = list(var = LGCP_VAR, scale = LGCP_SCALE),
                win = W_km, eps = EPS_LGCP, nsim = N_PILOT, drop = TRUE)
    TLl <- lapply(as_list(TL), to_poly)
    PIL <- data.frame()
    for (k in seq_along(TLl)) {
      Xk <- TLl[[k]]; t0 <- unname(proc.time()["elapsed"])
      r <- tryCatch(m3_test(Xk, NSIM_PILOT, bw, SEED7 + 100 + k),
                    error = function(e) { cat("  M3 FAILED:", conditionMessage(e), "\n"); NULL })
      sec <- unname(proc.time()["elapsed"]) - t0
      if (is.null(r)) next
      cat(sprintf("  rep %d: n %d (observed points dropped %d) | sigma2_hat %.3f alpha_hat %.2f | ERL p %.3f (script 34: %s) | usable %d | sim points dropped %d in %d sims | %.0f s\n",
                  k, npoints(Xk), r$drop_obs, r$sigma2, r$alpha, r$p, if (is.na(P34[k])) "failed" else sprintf("%.3f", P34[k]),
                  r$usable, r$drop_sim, r$sims_dropped, sec))
      PIL <- rbind(PIL, data.frame(rep = k, n = npoints(Xk), drop_obs = r$drop_obs, sigma2 = r$sigma2, alpha = r$alpha,
                                   p = r$p, p_script34 = P34[k], usable = r$usable, drop_sim = r$drop_sim,
                                   sims_dropped = r$sims_dropped, sec = sec))
    }
    if (nrow(PIL)) write.csv(PIL, file.path(OUT_DIR, "35_m3_pilot_rerun.csv"), row.names = FALSE)
    completed <- nrow(PIL) == N_PILOT
    cat(sprintf("\nM3 pilot completed %d of %d | usable sims all >= 90: %s | reproduces script-34 p for reps 1-3: %s\n",
                nrow(PIL), N_PILOT, if (nrow(PIL)) all(PIL$usable >= 90) else NA,
                if (nrow(PIL) >= 3) all(abs(PIL$p[1:3] - P34[1:3]) < 0.0051) else NA))
    
    ## ---------------- C. clarification log ----------------
    cat("\n== C. protocol_deviation_log.csv ==\n")
    D <- "2026-11-10"; AB <- "author-delegated (assistant decision; clarification, no frozen item changed)"
    ROWS <- data.frame(
      date = D,
      item = c("Clarification: Geyer fits use rbord = 0 (register value was 10 km)",
               "Clarification: M3 and M4 tests simulate conditional on n; M1 test unchanged",
               "Clarification: parametric bootstrap at point estimates in recovery and confirmatory tests",
               "Clarification: LGCP truth generator grid eps = 0.5 km accepted outside pre-set tolerance"),
      reason = c("rbord was a register value, not in freeze_v1. At 10 km complete separation occurred in 4 of 8 synthetic G3 fits (default vcov failed), and the active sites fell to about half of n. rbord 0 gave no separation in 8 of 8 (scripts 26-27).",
                 "Unconditional simulation from fitted Geyer/LGCP models carried an n spread of about +/-28% (scripts 26-27), which would widen envelopes spuriously. M3: rpoint on saved Lambda. M4: rmhcontrol(p = 1, expand = 1, n.start = n). Both are approximate conditional simulations (scripts 28-31).",
                 "Research_Design asks for parameter-uncertainty propagation where feasible. vcov of Geyer fits was singular in 4 of 8 fits at rbord 10 and vcov.kppm covers trend terms only; refit-based draws are not feasible at 19 cells x 100-200 reps x 99-199 sims (scripts 25-26). Bootstrap refits for gamma, r, sigma2, alpha are planned for the empirical stage.",
                 "Default rLGCP grid (1.41 km) did not match the 500 m covariates and produced excess pair correlation. eps = 0.5 km removed most of it. Strict rule (ratio 0.85-1.15 at r >= 2 km for both cells) still failed: L1 1.03-1.15 (unrounded above 1.15), L4 1.10-1.19, 10 sims at x4 intensity (script 33)."),
      effect_on_inference = c("M4 uses all sites and treats the window as closed (edge truncation of neighbour counts). rbord 5 is a required exploratory sensitivity. No archaeological model has been fitted.",
                              "M1 and M3/M4 tests differ in conditioning; size of each is checked by null calibration before the recovery run. Methods must state the asymmetry. No archaeological model has been fitted.",
                              "Plug-in envelopes ignore parameter uncertainty; direction of the effect on size is unknown and is measured by null calibration. Interaction claims need bootstrap-based gamma intervals.",
                              "Residual pair-correlation excess of about 3-19% at x4 intensity remains. sigma2_hat is about twice the nominal variance in L1; Gate 6 bracketing for M3 uses estimator output, nominal (var, scale) are cell labels only."),
      approved_by = AB, stringsAsFactors = FALSE)
    if (completed && sum(PIL$usable >= 90) == N_PILOT) {
      ROWS <- rbind(ROWS, data.frame(
        date = D, item = "Clarification: NA-covariate points removed from synthetic patterns before fitting",
        reason = "Script 34 M3 rep 4 failed: fast_fit dropped a data point (214 of 215). Synthetic patterns from rLGCP/rpoint are generated on the NA-filled field and can place points in pixels with NA covariates (script 35 diagnosis).",
        effect_on_inference = "Each synthetic pattern loses the few points in NA-covariate pixels (counts logged per test). Observed archaeological sites are checked for NA covariates (must be 0). No archaeological model has been fitted.",
        approved_by = AB, stringsAsFactors = FALSE))
    }
    if (!identical(names(LOG0), LOG_COLS)) {
      cat("Existing log columns differ from expected; NOT appending. Rows written to 35_log_rows_for_manual_append.csv\n")
      write.csv(ROWS, file.path(OUT_DIR, "35_log_rows_for_manual_append.csv"), row.names = FALSE)
    } else {
      NEW <- ROWS[!(ROWS$item %in% LOG0$item), , drop = FALSE]
      file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_35.csv"), overwrite = TRUE)
      write.csv(rbind(LOG0, NEW), f_log, row.names = FALSE)
      cat(sprintf("Appended %d row(s) to docs/protocol_deviation_log.csv (existing rows %d; backup in outputs).\n", nrow(NEW), nrow(LOG0)))
    }
    
    cat("\n== Warning tally (whole run) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "35_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")