
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

