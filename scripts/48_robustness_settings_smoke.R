# =====================================================================
# 48_robustness_settings_smoke.R   (v2: alignment tolerance 1e-3 km, rounding check)
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (OVERWRITES v1; writes data_derived/robust_grid_covs.rds)
# Purpose : build and smoke-test all robustness settings for the frozen M1 test before
#           any run: alignment of sites_NA.csv with X_km; covariate-builder validation
#           (500 m rebuild vs cached covs); 250 m and 1000 m covariate sets; subsets
#           (redeposited-excluded, Dataset C, LP, MP, UP); eroded window; E1, E2, E4;
#           quadrature 0.5 km; edge corrections; free z_S2 / z_S3; unwinsorised TPI (if found).
#           Generic M1 test code (own fit, curves, simulation) is checked against the
#           verified recovery code on the baseline pattern; each setting runs nsim = 29.
# Outputs : outputs/48_report.txt, 48_settings_smoke.csv, 48_sessionInfo.txt;
#           data_derived/robust_grid_covs.rds
# Notes   : Archaeological pattern used only for geometry, NA-covariate checks and timing;
#           smoke p-values (nsim 29) are NOT results and are printed only as completion checks.
#           No site is dropped to make a setting work: NA-covariate sites => BLOCKED.
#           Alignment rule (fixed before rerun): max |difference| between sites_NA.csv
#           coordinates and X_km < 1e-3 km (1 m); v1 used 1e-6 km, which was below the
#           0.1 m rounding of the CSV (observed maximum 5.00e-05 km).
#           Restart R before running. Main session only (terra is used here).
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))
suppressPackageStartupMessages(library(terra))

OUT_DIR <- dir_of("11_generative_gof", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
NSIM_SMOKE <- 29
ALIGN_TOL_KM <- 1e-3
SEED_S <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 40000
CONT <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5")
LITH <- c("lith_basalt", "lith_sandstone", "lith_carbonate", "lith_other")
EXPECT_N <- c(R00_baseline = 197, R01_redep = 195, R02_datasetC = 97, R03_LP = 75, R04_MP = 114, R05_UP = 79, R06_eroded5 = 169)

# ---- open all inputs BEFORE sink ----
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
SITES <- read.csv(f_sites, stringsAsFactors = FALSE)
LAY <- MAN[MAN$group == "derived_layer", , drop = FALSE]
LAYP <- gsub("\\\\", "/", LAY$path)
LAYNAMES <- basename(LAYP[grepl("/env/500m/", LAYP)])
ENV <- file.path(REPO_DIR, "data_derived", "env")
to_im2 <- function(r) {
  m <- terra::as.matrix(r, wide = TRUE); m <- m[nrow(m):1, , drop = FALSE]
  im(m, xcol = terra::xFromCol(r, seq_len(ncol(r))) / 1000, yrow = rev(terra::yFromRow(r, seq_len(nrow(r)))) / 1000)
}
build_covs <- function(dirp) {
  rl <- function(n) terra::rast(file.path(dirp, paste0(n, ".tif")))
  lay <- lapply(CONT, rl); names(lay) <- CONT
  l5 <- rl("lith5"); g3 <- rl("geomorph3")
  lay$lith_basalt <- terra::ifel(l5 == 1, 1, 0); lay$lith_sandstone <- terra::ifel(l5 == 3, 1, 0)
  lay$lith_carbonate <- terra::ifel(l5 == 4, 1, 0); lay$lith_other <- terra::ifel(l5 == 5, 1, 0)
  lay$geom_dissected <- terra::ifel(g3 == 2, 1, 0); lay$geom_alluvial <- terra::ifel(g3 == 3, 1, 0)
  lapply(lay, to_im2)
}
C500 <- build_covs(file.path(ENV, "500m")); C250 <- build_covs(file.path(ENV, "250m")); C1000 <- build_covs(file.path(ENV, "1000m"))
ZS2 <- to_im2(terra::rast(file.path(ENV, "500m", "z_S2.tif"))); ZS3 <- to_im2(terra::rast(file.path(ENV, "500m", "z_S3.tif")))
TPI_CAND <- LAYNAMES[grepl("tpi", LAYNAMES, ignore.case = TRUE) & !grepl("^z_tpi2km\\.tif$", LAYNAMES)]
TPIU <- if (length(TPI_CAND) == 1) to_im2(terra::rast(file.path(ENV, "500m", TPI_CAND))) else NULL
saveRDS(list(C250 = C250, C1000 = C1000, ZS2 = ZS2, ZS3 = ZS3, TPIU = TPIU, tpi_name = TPI_CAND),
        file.path(REPO_DIR, "data_derived", "robust_grid_covs.rds"))

# ---- generic M1 test code ----
WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
mk <- function(name, X, W, cv, dum, Lmin = 1, corrL = "isotropic", corrG = "Ripley") {
  n <- npoints(X); a <- area.owin(W); half <- PCF_STOYAN_C / sqrt(n / a); gmin <- max(Lmin, ceiling(half * 2) / 2)
  list(name = name, X = X, W = W, cov = cv, n = n, area = a, bw = half / sqrt(5),
       rL = RGRID[RGRID >= Lmin & RGRID <= R_MAX_DECISION_KM], rG = RGRID[RGRID >= gmin & RGRID <= R_MAX_DECISION_KM],
       corrL = corrL, corrG = corrG, dum = dum, covW = lapply(cv, function(im) im[W, drop = FALSE]))
}
cov_na <- function(P, S) {
  if (npoints(P) == 0) return(logical(0))
  Z <- do.call(cbind, lapply(S$cov, function(im) lookup.im(im, P$x, P$y, naok = TRUE)))
  rowSums(is.na(Z)) > 0
}
rob_clean <- function(P, S) { b <- cov_na(P, S); if (any(b)) P[!b] else P }
rob_poly <- function(P, S) { Q <- P[S$W]; if (Window(Q)$type != "polygonal") stop("not polygonal"); Q }
rob_fit <- function(Xp, S) {
  Q <- quadscheme(Xp, dummy = S$dum, method = "grid")
  w <- w.quad(Q); isd <- is.data(Q); xq <- x.quad(Q); yq <- y.quad(Q)
  Z <- vapply(S$cov, function(im) lookup.im(im, xq, yq, naok = TRUE), numeric(length(xq)))
  ok <- stats::complete.cases(Z) & w > 0
  X1 <- cbind(`(Intercept)` = 1, Z)
  g <- suppressWarnings(stats::glm.fit(X1[ok, , drop = FALSE], (isd / w)[ok], weights = w[ok], family = poisson(), intercept = FALSE))
  Xd <- X1[isd & ok, , drop = FALSE]
  if (nrow(Xd) != npoints(Xp)) stop("data points dropped in rob_fit: ", nrow(Xd), " of ", npoints(Xp))
  list(beta = g$coefficients, lam_data = as.numeric(exp(Xd %*% g$coefficients)))
}
rob_curves <- function(Xp, ft, S) {
  Lh <- do.call(Linhom, c(list(Xp, lambda = ft$lam_data, r = RGRID, correction = S$corrL), L_ARGS))
  gh <- do.call(pcfinhom, c(list(Xp, lambda = ft$lam_data, r = RGRID, correction = S$corrG, bw = S$bw), G_ARGS_NS))
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(S$rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(S$rG, 6), round(gh$r, 6))])
}
rob_test <- function(S, nsim, seed) {
  set.seed(seed)
  ft <- rob_fit(S$X, S); ob <- rob_curves(S$X, ft, S)
  lam <- exp(Reduce(`+`, Map(function(a, b) a * b, S$covW, ft$beta[-1])) + ft$beta[1])
  sL <- matrix(NA_real_, length(S$rL), nsim); sG <- matrix(NA_real_, length(S$rG), nsim); nn <- rep(NA_real_, nsim); nd <- 0
  for (k in seq_len(nsim)) {
    res <- tryCatch({ Xs <- rob_poly(rpoispp(lam), S); Xc <- rob_clean(Xs, S); fs <- rob_fit(Xc, S)
    list(cs = rob_curves(Xc, fs, S), n = npoints(Xc), d = npoints(Xs) - npoints(Xc)) }, error = function(e) NULL)
    if (!is.null(res)) { sL[, k] <- res$cs$L; sG[, k] <- res$cs$g; nn[k] <- res$n; nd <- nd + res$d }
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = S$rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = S$rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc), n_mean = mean(nn[okc]), dropped = nd)
}

REPORT <- file.path(OUT_DIR, "48_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("48_robustness_settings_smoke v2 | ", format(Sys.time()), "\n", sep = "")
    
    cat("\n== 0. Alignment and builder checks ==\n")
    dxv <- pmax(abs(X_obs_km$x - SITES$x_m / 1000), abs(X_obs_km$y - SITES$y_m / 1000))
    cat(sprintf("sites_NA.csv rows %d | coordinate difference to X_km (km): max %.2e (%.1f cm), 99th percentile %.2e, median %.2e | rule: max < %.0e km\n",
                nrow(SITES), max(dxv), max(dxv) * 1e5, quantile(dxv, 0.99), median(dxv), ALIGN_TOL_KM))
    r01 <- all(abs(SITES$x_m * 10 - round(SITES$x_m * 10)) < 1e-6) && all(abs(SITES$y_m * 10 - round(SITES$y_m * 10)) < 1e-6)
    cat("x_m and y_m stored at 0.1 m resolution:", r01, "\n")
    stopifnot(nrow(SITES) == 197, max(dxv) < ALIGN_TOL_KM)
    cat(sprintf("subset sizes: redeposited %d | Dataset C %d | LP %d | MP %d | UP %d\n", sum(SITES$sens_redeposited),
                sum(SITES$dataset_C), sum(SITES$LP), sum(SITES$MP), sum(SITES$UP)))
    same <- vapply(names(covs), function(nm) isTRUE(all(is.na(C500[[nm]]$v) == is.na(covs[[nm]]$v))) &&
                     isTRUE(max(abs(C500[[nm]]$v - covs[[nm]]$v), na.rm = TRUE) < 1e-9), logical(1))
    cat("500 m builder reproduces cached covariates for all", length(same), "layers:", all(same), "\n")
    if (!all(same)) cat("  differing layers:", paste(names(same)[!same], collapse = ", "), "\n")
    cat("env/500m layer names (", length(LAYNAMES), "): ", paste(LAYNAMES, collapse = ", "), "\n", sep = "")
    cat("unwinsorised-TPI candidate(s):", if (length(TPI_CAND)) paste(TPI_CAND, collapse = ", ") else "none",
        "| loaded:", !is.null(TPIU), "\n")
    if (!is.null(TPIU)) cat(sprintf("  candidate values: mean %.3f sd %.3f range %.2f to %.2f (cached z_tpi2km: mean %.3f sd %.3f)\n",
                                    mean(TPIU$v, na.rm = TRUE), sd(TPIU$v, na.rm = TRUE), min(TPIU$v, na.rm = TRUE), max(TPIU$v, na.rm = TRUE),
                                    mean(covs$z_tpi2km$v, na.rm = TRUE), sd(covs$z_tpi2km$v, na.rm = TRUE)))
    cat("grid sizes (pixels): 250 m", paste(dim(C250[[1]]$v), collapse = "x"), "| 500 m", paste(dim(C500[[1]]$v), collapse = "x"),
        "| 1000 m", paste(dim(C1000[[1]]$v), collapse = "x"), "\n")
    cat("pcfinhom correction choices:", tryCatch(paste(eval(formals(getS3method("pcfinhom", "ppp"))$correction), collapse = ", "), error = function(e) "not readable"), "\n")
    cat("Linhom formals:", paste(names(formals(Linhom)), collapse = ", "), "\n")
    
    cat("\n== 1. Settings ==\n")
    Wer <- erosion(W_km, 5)
    cat(sprintf("eroded window: type %s, area %.1f km2\n", Wer$type, area.owin(Wer)))
    Xe <- X_obs_km[Wer]
    DUMq <- default.dummy(X_obs_km, eps = 0.5); DUMe <- default.dummy(Xe, eps = QUAD_KM)
    ST <- list()
    ST$R00_baseline <- mk("R00_baseline", X_obs_km, W_km, covs, DUMMY_FIXED)
    ST$R01_redep    <- mk("R01_redep", X_obs_km[!SITES$sens_redeposited], W_km, covs, DUMMY_FIXED)
    ST$R02_datasetC <- mk("R02_datasetC", X_obs_km[SITES$dataset_C], W_km, covs, DUMMY_FIXED, Lmin = 2)
    ST$R03_LP       <- mk("R03_LP", X_obs_km[SITES$LP], W_km, covs, DUMMY_FIXED)
    ST$R04_MP       <- mk("R04_MP", X_obs_km[SITES$MP], W_km, covs, DUMMY_FIXED)
    ST$R05_UP       <- mk("R05_UP", X_obs_km[SITES$UP], W_km, covs, DUMMY_FIXED)
    ST$R06_eroded5  <- mk("R06_eroded5", Xe, Wer, covs, DUMe)
    ST$R07_E1       <- mk("R07_E1", X_obs_km, W_km, covs[CONT], DUMMY_FIXED)
    ST$R08_E2       <- mk("R08_E2", X_obs_km, W_km, covs[c(CONT, LITH)], DUMMY_FIXED)
    ST$R09_E4       <- mk("R09_E4", X_obs_km, W_km, covs[c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", LITH)], DUMMY_FIXED)
    ST$R10_grid250  <- mk("R10_grid250", X_obs_km, W_km, C250, DUMMY_FIXED)
    ST$R11_grid1000 <- mk("R11_grid1000", X_obs_km, W_km, C1000, DUMMY_FIXED)
    ST$R12_quad05   <- mk("R12_quad05", X_obs_km, W_km, covs, DUMq)
    ST$R13_edge_translate <- mk("R13_edge_translate", X_obs_km, W_km, covs, DUMMY_FIXED, corrL = "translate", corrG = "translate")
    ST$R14_edge_border_L  <- mk("R14_edge_border_L", X_obs_km, W_km, covs, DUMMY_FIXED, corrL = "border", corrG = "translate")
    ST$R15_zS2      <- mk("R15_zS2", X_obs_km, W_km, c(covs, list(z_S2 = ZS2)), DUMMY_FIXED)
    ST$R16_zS3      <- mk("R16_zS3", X_obs_km, W_km, c(covs, list(z_S3 = ZS3)), DUMMY_FIXED)
    if (!is.null(TPIU)) { cvt <- covs; cvt$z_tpi2km <- TPIU; ST$R17_tpi_unwins <- mk("R17_tpi_unwins", X_obs_km, W_km, cvt, DUMMY_FIXED) }
    TB <- do.call(rbind, lapply(ST, function(S) data.frame(setting = S$name, n = S$n, expected = if (S$name %in% names(EXPECT_N)) EXPECT_N[[S$name]] else NA,
                                                           area_km2 = S$area, df = length(S$cov), L_range = sprintf("%.1f-%.1f", min(S$rL), max(S$rL)), g_range = sprintf("%.1f-%.1f", min(S$rG), max(S$rG)),
                                                           bw_km = S$bw, corrL = S$corrL, corrG = S$corrG, na_sites = sum(cov_na(S$X, S)), stringsAsFactors = FALSE)))
    print(TB, digits = 4, row.names = FALSE)
    bad <- TB$setting[!is.na(TB$expected) & TB$n != TB$expected]
    cat("site counts differing from the frozen counts:", if (length(bad)) paste(bad, collapse = ", ") else "none", "\n")
    
    cat("\n== 2. Generic code vs verified recovery code (baseline observed pattern) ==\n")
    S0 <- ST$R00_baseline; ft <- rob_fit(S0$X, S0); ob <- rob_curves(S0$X, ft, S0)
    f_ref <- fast_fit(X_obs_km); ref <- fast_curves(X_obs_km, f_ref, bw_for_n(N_ANALYTICAL))
    cat(sprintf("max |coef difference| %.2e | max |L difference| %.2e | max |g difference| %.2e (expect ~0)\n",
                max(abs(ft$beta - f_ref$beta)), max(abs(ob$L - ref$L)), max(abs(ob$g - ref$g))))
    
    cat("\n== 3. Smoke test (nsim ", NSIM_SMOKE, "; completion check only, p-values are NOT results) ==\n", sep = "")
    SM <- data.frame()
    for (nm in names(ST)) {
      S <- ST[[nm]]; t0 <- unname(proc.time()["elapsed"])
      if (TB$na_sites[TB$setting == nm] > 0) {
        cat(sprintf("  %-20s BLOCKED: %d observed site(s) with NA covariates\n", nm, TB$na_sites[TB$setting == nm]))
        SM <- rbind(SM, data.frame(setting = nm, status = "BLOCKED_NA_SITES", sec = NA, usable = NA, n_mean = NA, dropped = NA)); next
      }
      r <- tryCatch(rob_test(S, NSIM_SMOKE, SEED_S + match(nm, names(ST))), error = function(e) { cat(sprintf("  %-20s FAILED: %s\n", nm, conditionMessage(e))); NULL })
      sec <- unname(proc.time()["elapsed"]) - t0
      if (is.null(r)) { SM <- rbind(SM, data.frame(setting = nm, status = "FAILED", sec = sec, usable = NA, n_mean = NA, dropped = NA)); next }
      cat(sprintf("  %-20s ok | %5.1f s | usable %d of %d | simulated n mean %.1f | cleaned points %d | p (smoke) %.3f\n", nm, sec, r$usable, NSIM_SMOKE, r$n_mean, r$dropped, r$p))
      SM <- rbind(SM, data.frame(setting = nm, status = "OK", sec = sec, usable = r$usable, n_mean = r$n_mean, dropped = r$dropped))
    }
    write.csv(merge(TB, SM, by = "setting"), file.path(OUT_DIR, "48_settings_smoke.csv"), row.names = FALSE)
    cat(sprintf("\nsettings OK %d | BLOCKED %d | FAILED %d\n", sum(SM$status == "OK"), sum(SM$status == "BLOCKED_NA_SITES"), sum(SM$status == "FAILED")))
    cat(sprintf("projected nsim-999 time per setting (min): %s\n", paste(sprintf("%s %.0f", SM$setting[SM$status == "OK"], SM$sec[SM$status == "OK"] / NSIM_SMOKE * 999 / 60), collapse = " | ")))
    cat("\n== Warning tally (main session) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "48_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")