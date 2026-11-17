
# recovery_setup.R -- written by 17_stats_profile.R; sourced by recovery scripts.
# Requires REPO_DIR and MASTER_SEED from R/config.R.
params_f <- file.path(REPO_DIR, "R", "freeze_v1_params.R"); stopifnot(file.exists(params_f)); source(params_f)
if (!grepl("^FROZEN", FREEZE_STATUS)) stop("Freeze v1 is not FROZEN. Aborting.")
need <- c("spatstat.geom", "spatstat.explore", "spatstat.model", "spatstat.random", "GET", "terra", "sf")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: install.packages(c(", paste(sprintf('"%s"', miss), collapse = ", "), "))")
suppressPackageStartupMessages({ library(sf); library(spatstat.geom); library(spatstat.explore)
                                 library(spatstat.model); library(spatstat.random); library(GET) })
QUAD_KM <- QUADRATURE_KM$primary
RGRID <- seq(0, R_MAX_DECISION_KM, by = R_STEP_KM)
rL <- RGRID[RGRID >= R_MIN_KM["A"] & RGRID <= R_MAX_DECISION_KM]
rG <- RGRID[RGRID >= GINHOM_LOWER_KM["A"] & RGRID <= R_MAX_DECISION_KM]
env_dir <- file.path(REPO_DIR, "data_derived", "env", "500m")
obj <- readRDS(file.path(REPO_DIR, "data_derived", "ppp_NA_window.rds"))
X_obs_km <- obj$X_km; W_km <- Window(X_obs_km); stopifnot(W_km$type == "polygonal")
rl <- function(n) terra::rast(file.path(env_dir, paste0(n, ".tif")))
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
to_poly <- function(P, label = "") {
  Q <- P[W_km]
  if (Window(Q)$type != "polygonal") stop("pattern window is not polygonal after restriction")
  Q
}
make_truth <- function(env_mult = 1, n_target = RECOVERY_ENV$expected_n) {
  set.seed(SEEDS$master + RECOVERY_ENV$seed_offset)
  beta <- setNames(rnorm(length(covs), 0, RECOVERY_ENV$coef_sd), names(covs)) * env_mult
  eta <- Reduce(`+`, Map(function(a, b) a * b, covs, beta))
  lam0 <- exp(eta)[W_km, drop = FALSE]
  b0 <- log(n_target / integral(lam0))
  list(beta = beta, b0 = b0, lamW = lam0 * exp(b0))
}
col_of <- function(f, pref) { for (nm in pref) if (nm %in% names(f)) return(f[[nm]]); stop("no column among: ", paste(names(f), collapse = ",")) }
run <- function(label, expr) {
  cat(sprintf("%-60s", label)); flush.console()
  t <- system.time(r <- try(expr, silent = TRUE))["elapsed"]; ok <- !inherits(r, "try-error")
  cat(sprintf("%8.2f s  %s\n", t, if (ok) "ok" else paste("FAILED:", conditionMessage(attr(r, "condition")))))
  list(res = if (ok) r else NULL, sec = unname(t), ok = ok)
}
sec <- function(label, expr) tryCatch(expr, error = function(e) cat(sprintf("\n[SECTION FAILED: %s] %s\n", label, conditionMessage(e))))

