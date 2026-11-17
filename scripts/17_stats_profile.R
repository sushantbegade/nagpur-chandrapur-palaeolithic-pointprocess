# =============================================================================
# 17_stats_profile.R
# Stage: 10_simulation_recovery | cost reduction for the summary-statistic step
# Purpose: (1) write R/recovery_setup.R shared by recovery scripts; (2) find which part
#   of fit + L_inhom + g_inhom is slow; (3) compare explicit pcfinhom settings;
#   (4) compare fast vs slow pipeline; (5) parallel speed-up; (6) Geyer truth
#   calibration (n ~ 197) and rmh convergence. SYNTHETIC patterns only.
# Outputs: R/recovery_setup.R, 10_simulation_recovery/outputs/stats_profile_report.txt,
#   stats_profile.csv, geyer_calibration.csv, sessionInfo_17.txt
# Not executed by assistant. Run top to bottom in RStudio (about 10-20 minutes).
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"

SETUP_TEXT <- r"---(
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
)---"

while (sink.number() > 0) sink()
source(REPO_CONFIG)
setup_f <- file.path(REPO_DIR, "R", "recovery_setup.R")
writeLines(SETUP_TEXT, setup_f)
source(setup_f)
out <- dir_of("10_simulation_recovery"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
maxdiff <- function(a, b) suppressWarnings(max(abs(a - b), na.rm = TRUE))
tab <- list()

sink(file.path(out, "stats_profile_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== A. INSTALLED DEFAULTS ====\n")
  for (p in c("spatstat.geom", "spatstat.explore", "spatstat.model", "spatstat.random", "GET")) cat(sprintf("%-18s %s\n", p, as.character(packageVersion(p))))
  fm <- formals(pcfinhom)
  cat("pcfinhom defaults:\n"); print(sapply(fm[c("kernel", "bw", "stoyan", "correction", "divisor", "zerocor", "renormalise", "normpower", "update", "leaveoneout", "fast")],
                                            function(x) paste(deparse(x), collapse = " ")))
  cat("fitted.ppm args:", paste(names(formals(getS3method("fitted", "ppm"))), collapse = ", "), "\n")
  cat("Kinhom args    :", paste(names(formals(Kinhom)), collapse = ", "), "\n")
  
  cat("\n==== B. SYNTHETIC PATTERN AND FIT ====\n")
  truth <- make_truth()
  set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 1)
  Xp <- to_poly(rpoispp(truth$lamW)); cat("n =", npoints(Xp), "\n")
  tt <- replicate(3, system.time(fit <- ppm(Xp, trend = formula_E3, covariates = covs, quad.args = list(eps = QUAD_KM)))["elapsed"])
  cat(sprintf("trend-only ppm fit (quad %g km): median %.2f s over 3\n", QUAD_KM, median(tt)))
  tab[[1]] <- data.frame(step = "ppm trend-only fit", variant = "quad 1 km", sec = median(tt), maxdiff = NA)
  lam_loo <- tryCatch(as.numeric(fitted(fit, dataonly = TRUE, leaveoneout = TRUE)), error = function(e) { cat("fitted(leaveoneout) FAILED:", conditionMessage(e), "\n"); NULL })
  t_loo <- system.time(try(fitted(fit, dataonly = TRUE, leaveoneout = TRUE), silent = TRUE))["elapsed"]
  lam_plain <- tryCatch(as.numeric(fitted(fit, dataonly = TRUE)), error = function(e) { cat("fitted(plain) FAILED:", conditionMessage(e), "\n"); NULL })
  t_plain <- system.time(try(fitted(fit, dataonly = TRUE), silent = TRUE))["elapsed"]
  cat(sprintf("intensity at data points: leave-one-out %.2f s | plain %.2f s", t_loo, t_plain))
  if (!is.null(lam_loo) && !is.null(lam_plain)) cat(sprintf(" | max relative difference %.4f\n", max(abs(lam_loo / lam_plain - 1)))) else cat("\n")
  if (is.null(lam_loo)) lam_loo <- lam_plain
  tab[[2]] <- data.frame(step = "fitted() at data points", variant = c("leaveoneout", "plain"), sec = c(t_loo, t_plain), maxdiff = NA)
  
  cat("\n==== C. L_inhom VARIANTS ====\n")
  Lcalc <- function(Xp, lam, ...) { Lh <- Linhom(Xp, lambda = lam, r = RGRID, correction = "isotropic", ...)
  (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(rL, 6), round(Lh$r, 6))] }
  vL <- list("model lambda, update=FALSE" = function() Lcalc(Xp, fit, update = FALSE),
             "numeric leave-one-out"      = function() Lcalc(Xp, lam_loo),
             "numeric plain"              = function() Lcalc(Xp, lam_plain))
  Lres <- list(); rows <- list()
  for (nm in names(vL)) { r <- run(paste("Linhom:", nm), vL[[nm]]()); Lres[[nm]] <- r$res
  rows[[length(rows) + 1]] <- data.frame(step = "L_inhom", variant = nm, sec = r$sec, maxdiff = if (r$ok && !is.null(Lres[[1]])) maxdiff(r$res, Lres[[1]]) else NA) }
  print(do.call(rbind, rows), row.names = FALSE); tab[[3]] <- do.call(rbind, rows)
  
  cat("\n==== D. g_inhom VARIANTS (numeric leave-one-out intensity) ====\n")
  gcalc <- function(Xp, lam, args) { gh <- do.call(pcfinhom, c(list(Xp, lambda = lam, r = RGRID, correction = "Ripley"), args))
  list(g = col_of(gh, c("iso", "trans"))[match(round(rG, 6), round(gh$r, 6))], bw = as.numeric(attr(gh, "bw"))[1]) }
  gv <- list(default                = list(),
             weighted_r_np1         = list(zerocor = "weighted", divisor = "r", normpower = 1),
             weighted_r_np2         = list(zerocor = "weighted", divisor = "r", normpower = 2),
             none_r_np2             = list(zerocor = "none", divisor = "r", normpower = 2),
             convolution_r_np2      = list(zerocor = "convolution", divisor = "r", normpower = 2),
             reflection_r_np2       = list(zerocor = "reflection", divisor = "r", normpower = 2))
  Gres <- list(); rows <- list()
  for (nm in names(gv)) { r <- run(paste("pcfinhom:", nm), gcalc(Xp, lam_loo, gv[[nm]])); Gres[[nm]] <- r$res
  rows[[length(rows) + 1]] <- data.frame(step = "g_inhom", variant = nm, sec = r$sec, bw_km = if (r$ok) round(r$res$bw, 3) else NA,
                                         NA_in_g = if (r$ok) sum(is.na(r$res$g)) else NA, maxdiff = if (r$ok && !is.null(Gres$default)) maxdiff(r$res$g, Gres$default$g) else NA) }
  gt <- do.call(rbind, rows); print(gt, row.names = FALSE)
  tab[[4]] <- gt[, c("step", "variant", "sec", "maxdiff")]
  okg <- gt$variant != "default" & !is.na(gt$sec) & gt$NA_in_g == 0 & !is.na(gt$NA_in_g)
  g_choice <- if (any(okg)) gt$variant[okg][which.min(gt$sec[okg])] else "weighted_r_np1"
  cat("fastest explicit g_inhom setting without NA:", g_choice, "\n")
  
  cat("\n==== E. FAST vs SLOW PIPELINE (one pattern) ====\n")
  make_fast <- function(g_args, L_args = list(renormalise = TRUE, normpower = 1)) {
    force(g_args); force(L_args)
    function(Xp) {
      fit <- ppm(Xp, trend = formula_E3, covariates = covs, quad.args = list(eps = QUAD_KM))
      lam <- tryCatch(as.numeric(fitted(fit, dataonly = TRUE, leaveoneout = TRUE)), error = function(e) as.numeric(fitted(fit, dataonly = TRUE)))
      Lh <- do.call(Linhom, c(list(Xp, lambda = lam, r = RGRID, correction = "isotropic"), L_args))
      gh <- do.call(pcfinhom, c(list(Xp, lambda = lam, r = RGRID, correction = "Ripley"), g_args))
      Lv <- col_of(Lh, c("iso", "trans", "border")) - Lh$r; gv <- col_of(gh, c("iso", "trans"))
      list(L = Lv[match(round(rL, 6), round(Lh$r, 6))], g = gv[match(round(rG, 6), round(gh$r, 6))], bw = as.numeric(attr(gh, "bw"))[1])
    }
  }
  curves_slow <- function(Xp) { fit <- ppm(Xp, trend = formula_E3, covariates = covs, quad.args = list(eps = QUAD_KM))
  Lh <- Linhom(Xp, lambda = fit, r = RGRID, correction = "isotropic", update = FALSE)
  gh <- pcfinhom(Xp, lambda = fit, r = RGRID, correction = "Ripley", update = FALSE)
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(rG, 6), round(gh$r, 6))]) }
  fast <- make_fast(gv[[g_choice]])
  rs <- run("slow pipeline (as in script 16)", curves_slow(Xp))
  rf <- run(paste("fast pipeline, g setting:", g_choice), fast(Xp))
  if (rs$ok && rf$ok) cat(sprintf("speed-up x%.1f | max |dL| %.4f | max |dg| %.4f | g range %.2f-%.2f\n", rs$sec / rf$sec,
                                  maxdiff(rs$res$L, rf$res$L), maxdiff(rs$res$g, rf$res$g), min(rf$res$g, na.rm = TRUE), max(rf$res$g, na.rm = TRUE)))
  tab[[5]] <- data.frame(step = "full pipeline", variant = c("slow", paste("fast", g_choice)), sec = c(rs$sec, rf$sec),
                         maxdiff = c(NA, if (rs$ok && rf$ok) maxdiff(rs$res$g, rf$res$g) else NA))
  
  cat("\n==== F. PARALLEL SPEED-UP (fast pipeline, 8 patterns) ====\n")
  sec("F parallel", {
    patterns <- lapply(1:8, function(i) { set.seed(SEEDS$master + SEEDS$offsets["bootstrap"] + 100 + i); to_poly(simulate(fit, nsim = 1, drop = TRUE, verbose = FALSE)) })
    ncl <- min(4, parallel::detectCores(logical = FALSE))
    ts <- system.time(ser <- lapply(patterns, fast))["elapsed"]
    cl <- parallel::makeCluster(ncl)
    tp <- tryCatch({
      parallel::clusterExport(cl, c("covs", "W_km", "formula_E3", "QUAD_KM", "RGRID", "rL", "rG", "col_of"), envir = environment())
      parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(spatstat.geom); library(spatstat.explore); library(spatstat.model); library(spatstat.random) }))
      t0 <- system.time(par_res <- parallel::parLapply(cl, patterns, fast))["elapsed"]
      list(t = t0, ok = all(sapply(par_res, function(x) is.list(x) && !anyNA(x$L))), d = max(sapply(seq_along(ser), function(i) maxdiff(ser[[i]]$g, par_res[[i]]$g))))
    }, error = function(e) { cat("parallel FAILED:", conditionMessage(e), "\n"); NULL }, finally = parallel::stopCluster(cl))
    cat(sprintf("serial %.1f s | %d workers", ts, ncl))
    if (!is.null(tp)) cat(sprintf(" %.1f s | speed-up x%.2f | identical results: %s (max diff %.2g)\n", tp$t, ts / tp$t, tp$ok && tp$d < 1e-8, tp$d)) else cat("\n")
    tab[[6]] <<- data.frame(step = "8 patterns", variant = c("serial", paste0("parallel x", ncl)), sec = c(ts, if (!is.null(tp)) tp$t else NA), maxdiff = NA)
  })
  
  cat("\n==== G. GEYER TRUTH: TREND CALIBRATION (n ~ 197) AND rmh CONVERGENCE ====\n")
  sec("G Geyer", {
    geyer_run <- function(mult, nrep, seed) {
      set.seed(seed)
      mod <- rmhmodel(cif = "geyer", par = list(beta = 1, gamma = 2.5, r = 3, sat = 3), w = W_km, trend = truth$lamW * mult)
      t <- system.time(P <- rmh(mod, start = list(n.start = 197), control = list(nrep = nrep, expand = 1), verbose = FALSE))["elapsed"]
      c(n = npoints(P[W_km]), sec = unname(t))
    }
    mults <- c(0.10, 0.14, 0.18); cal <- list()
    for (m in mults) for (s in 1:3) { x <- geyer_run(m, 5e4, SEEDS$master + SEEDS$offsets["recovery"] + 200 + s); cal[[length(cal) + 1]] <- data.frame(mult = m, nrep = 5e4, rep = s, n = x["n"], sec = x["sec"]) }
    cal <- do.call(rbind, cal); agg <- aggregate(cbind(n, sec) ~ mult, data = cal, FUN = mean); print(round(agg, 2), row.names = FALSE)
    fitlm <- lm(log(n) ~ log(mult), data = agg); m_star <- exp((log(197) - coef(fitlm)[1]) / coef(fitlm)[2])
    cat(sprintf("interpolated trend multiplier for n = 197: %.3f (rough)\n", m_star))
    conv <- list()
    for (nr in c(2e4, 5e4, 1e5)) for (s in 1:3) { x <- geyer_run(m_star, nr, SEEDS$master + SEEDS$offsets["recovery"] + 300 + s); conv[[length(conv) + 1]] <- data.frame(mult = m_star, nrep = nr, rep = s, n = x["n"], sec = x["sec"]) }
    conv <- do.call(rbind, conv); ca <- aggregate(cbind(n, sec) ~ nrep, data = conv, FUN = mean)
    cat("convergence check at the interpolated multiplier (n should not drift with nrep):\n"); print(round(ca, 2), row.names = FALSE)
    write.csv(rbind(cal, conv), file.path(out, "geyer_calibration.csv"), row.names = FALSE)
  })
  
  write.csv(do.call(rbind, lapply(tab, function(d) d[, intersect(names(d), c("step", "variant", "sec", "maxdiff"))])), file.path(out, "stats_profile.csv"), row.names = FALSE)
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_17.txt"))
cat("\nDONE. Outputs in:", out, "\n")