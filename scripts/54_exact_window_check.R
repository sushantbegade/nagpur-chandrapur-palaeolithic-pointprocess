# =====================================================================
# 54_exact_window_check.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW)
# Purpose : reduced exact-window sensitivity of the frozen M1 test (exact boundary polygon vs
#           the 100 m simplified window used for all simulations):
#           (A) window geometry (vertices, area, symmetric-difference area, sites inside);
#           (B) observed pattern: E3 fit and L_inhom / g_inhom curves on both windows, difference
#               relative to the null sd of the 999 saved simplified-window simulations (43_m1_sims.rds);
#           (C) ERL p with the exact-window observed curves against the saved simulations;
#           (D) up to 12 regenerated simulated patterns: curves on both windows, same comparison.
# Outputs : outputs/54_report.txt, 54_window_check.csv, 54_sessionInfo.txt
# Notes   : Exact window = vertex-rich polygon from data_derived/recovery_cache.rds (W_km) or
#           ppp_NA_window.rds (W_m rescaled to km); the one with more than 5000 vertices.
#           Reading rules (set before running): window difference immaterial if max over radii of
#           |curve difference| / null sd < 0.05 for the observed pattern and all simulated patterns,
#           and |p(exact observed vs saved sims) - 0.0310| <= 0.005; otherwise escalate to a full
#           exact-window M1 test at nsim 199 on workers. Bandwidth is held at the simplified-window
#           value to isolate the window effect. Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))
GENERIC_FILE <- file.path(REPO_DIR, "R", "robust_generic.R"); stopifnot(file.exists(GENERIC_FILE))
source(GENERIC_FILE)

OUT_DIR <- dir_of("11_generative_gof", "outputs")
f_sims <- file.path(OUT_DIR, "43_m1_sims.rds"); f_r42 <- file.path(OUT_DIR, "42_m1_result.csv")
stopifnot(file.exists(f_sims), file.exists(f_r42))
f_c0 <- file.path(REPO_DIR, "data_derived", "recovery_cache.rds")
f_pp <- file.path(REPO_DIR, "data_derived", "ppp_NA_window.rds"); stopifnot(file.exists(f_pp))
SEED_W <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 110000
NS_MAX <- 12; BUDGET_S <- 1500
SM <- readRDS(f_sims); R42 <- read.csv(f_r42, stringsAsFactors = FALSE)
nv <- function(W) sum(vapply(W$bdry, function(b) length(b$x), 1L))
CAND <- list()
if (file.exists(f_c0)) CAND$recovery_cache <- tryCatch(readRDS(f_c0)$W_km, error = function(e) NULL)
OBJ <- readRDS(f_pp); CAND$ppp_W_m <- tryCatch(rescale(OBJ$W_m, 1000), error = function(e) NULL)
CAND <- CAND[!vapply(CAND, is.null, logical(1))]

REPORT <- file.path(OUT_DIR, "54_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("54_exact_window_check | ", format(Sys.time()), "\n", sep = "")
  cat("window route note in ppp object:", if (!is.null(OBJ$window_route)) OBJ$window_route else "none", "\n")
  
  cat("\n== A. Windows ==\n")
  cat(sprintf("simplified window (W_km): %d vertices, area %.2f km2\n", nv(W_km), area.owin(W_km)))
  for (nm in names(CAND)) cat(sprintf("candidate %-15s: type %s, %d vertices, area %.2f km2\n", nm, CAND[[nm]]$type, nv(CAND[[nm]]), area.owin(CAND[[nm]])))
  vn <- vapply(CAND, nv, 1L); stopifnot(length(CAND) > 0, max(vn) > 5000)
  W_ex <- CAND[[which.max(vn)]]; cat("exact window used:", names(CAND)[which.max(vn)], "\n")
  stopifnot(W_ex$type == "polygonal")
  X_ex <- X_obs_km[W_ex]
  cat(sprintf("sites inside the exact window: %d of %d\n", npoints(X_ex), npoints(X_obs_km))); stopifnot(npoints(X_ex) == 197)
  sd_area <- tryCatch(area.owin(setminus.owin(W_ex, W_km)) + area.owin(setminus.owin(W_km, W_ex)), error = function(e) NA_real_)
  cat(sprintf("symmetric-difference area %s km2\n", if (is.na(sd_area)) "not computed" else sprintf("%.3f (%.4f%% of the window)", sd_area, 100 * sd_area / area.owin(W_km))))
  
  cat("\n== B. Observed pattern on both windows ==\n")
  S_s <- mk("SIMPL", X_obs_km, W_km, covs, DUMMY_FIXED)
  t0 <- unname(proc.time()["elapsed"])
  S_x <- mk("EXACT", X_ex, W_ex, covs, default.dummy(X_ex, eps = QUAD_KM)); S_x$bw <- S_s$bw
  t_build <- unname(proc.time()["elapsed"]) - t0
  t1 <- unname(proc.time()["elapsed"])
  ft_x <- rob_fit(X_ex, S_x); ob_x <- rob_curves(X_ex, ft_x, S_x)
  t_obs <- unname(proc.time()["elapsed"]) - t1
  ft_s <- rob_fit(X_obs_km, S_s); ob_s <- rob_curves(X_obs_km, ft_s, S_s)
  cat(sprintf("exact-window setting built in %.0f s; fit + curves in %.0f s\n", t_build, t_obs))
  cat(sprintf("simplified-window curves equal the saved script-43 observed curves: max |L diff| %.2e, max |g diff| %.2e\n",
              max(abs(ob_s$L - SM$ob$L)), max(abs(ob_s$g - SM$ob$g))))
  sdL <- apply(SM$sL, 1, sd); sdG <- apply(SM$sG, 1, sd)
  dL <- abs(ob_x$L - ob_s$L) / sdL; dG <- abs(ob_x$g - ob_s$g) / sdG
  cat(sprintf("E3 coefficients: max |difference| %.2e | fitted intensity at sites: max relative difference %.2e\n",
              max(abs(ft_x$beta - ft_s$beta)), max(abs(ft_x$lam_data / ft_s$lam_data - 1))))
  cat(sprintf("observed curves exact vs simplified: max |difference| / null sd: L %.4f (r = %.1f km), g %.4f (r = %.1f km)\n",
              max(dL), S_s$rL[which.max(dL)], max(dG), S_s$rG[which.max(dG)]))
  
  cat("\n== C. ERL p with exact-window observed curves vs the 999 saved simulations ==\n")
  csL <- GET::create_curve_set(list(r = S_s$rL, obs = ob_x$L, sim_m = SM$sL)); csG <- GET::create_curve_set(list(r = S_s$rG, obs = ob_x$g, sim_m = SM$sG))
  p_x <- erl_p(csL, csG, 2)
  cat(sprintf("p (exact-window observed curves) %.4f | script 42 (simplified) %.4f | difference %.4f\n", p_x, R42$p_erl_nstep2, p_x - R42$p_erl_nstep2))
  
  cat("\n== D. Simulated patterns on both windows ==\n")
  NS <- max(0, min(NS_MAX, floor(BUDGET_S / max(t_obs, 1)) - 1))
  cat(sprintf("patterns to compare: %d (budget %d s, %.0f s per exact-window evaluation)\n", NS, BUDGET_S, t_obs))
  f0 <- fast_fit(X_obs_km); lam <- lam_im_of(f0$beta); set.seed(SEED_W)
  ROWS <- data.frame(id = 0, type = "observed", sec = t_obs, pts_lost = 0L, max_dz_L = max(dL), max_dz_g = max(dG))
  for (k in seq_len(NS)) {
    P <- to_poly(rpoispp(lam)); Px <- P[W_ex]; t2 <- unname(proc.time()["elapsed"])
    cs <- rob_curves(P, rob_fit(P, S_s), S_s); cx <- rob_curves(Px, rob_fit(Px, S_x), S_x)
    ROWS <- rbind(ROWS, data.frame(id = k, type = "simulated", sec = unname(proc.time()["elapsed"]) - t2, pts_lost = npoints(P) - npoints(Px),
                                   max_dz_L = max(abs(cx$L - cs$L) / sdL), max_dz_g = max(abs(cx$g - cs$g) / sdG)))
    cat(sprintf("  pattern %2d: n %d, points lost %d | max |diff|/sd: L %.4f g %.4f | %.0f s\n", k, npoints(P), npoints(P) - npoints(Px),
                ROWS$max_dz_L[nrow(ROWS)], ROWS$max_dz_g[nrow(ROWS)], ROWS$sec[nrow(ROWS)]))
  }
  write.csv(ROWS, file.path(OUT_DIR, "54_window_check.csv"), row.names = FALSE)
  
  cat("\n== Verdict (pre-set rule) ==\n")
  mx <- max(ROWS$max_dz_L, ROWS$max_dz_g)
  ok1 <- mx < 0.05; ok2 <- abs(p_x - R42$p_erl_nstep2) <= 0.005
  cat(sprintf("max |difference| / null sd over observed and %d simulated patterns: %.4f (< 0.05: %s) | |p difference| %.4f (<= 0.005: %s)\n", NS, mx, ok1, abs(p_x - R42$p_erl_nstep2), ok2))
  cat("verdict:", if (ok1 && ok2) "WINDOW DIFFERENCE IMMATERIAL (reduced design)" else "ESCALATE: full exact-window M1 test at nsim 199 on workers", "\n")
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "54_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")