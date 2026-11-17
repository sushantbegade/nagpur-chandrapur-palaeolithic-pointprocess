# =====================================================================
# 49_robustness_m1_runs.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW)
# Purpose : frozen M1 test (E3 Poisson, joint ERL L_inhom + g_inhom, nstep 2, unconditional
#           rpoispp, nsim 999) on the archaeological pattern under the robustness settings:
#           R00 control (baseline, new seed), R01 redeposited-excluded (195), R02 Dataset C (97),
#           R03-R05 chronology LP/MP/UP, R06 eroded window 5 km (169), R07 E1, R08 E2, R09 E4,
#           R10 grid 250 m, R11 grid 1000 m, R12 quadrature 0.5 km, R13 translate edge correction,
#           R14 border (L) + translate (g), R15 free z_S2, R16 free z_S3.
#           Per setting: p, Monte Carlo SE, usable simulations, pair counts (observed vs
#           simulated, 1/2/3/5 km), raw max |z| location, warnings, curves.
# Outputs : outputs/49_report.txt, 49_res_<setting>.csv, 49_curves_<setting>.csv,
#           49_robustness_summary.csv, 49_curves_all.csv, 49_sessionInfo.txt,
#           protocol_deviation_log_before_49.csv; docs/protocol_deviation_log.csv (1 row)
# Notes   : g_inhom lower bound = max(frozen rule, 2.0 km) in every setting (changes only the
#           eroded window). No terra: grid covariates come from data_derived/robust_grid_covs.rds.
#           Unwinsorised-TPI setting deferred (layer not saved). Resumable.
#           Reading rules (set before running): p < 0.04 deficiency retained; p in [0.04, 0.06]
#           precision rule (script 49b, nsim 10000); p > 0.06 no deficiency detected, and
#           INCONCLUSIVE if n <= 114 (Stage 1 M1 power 0.13-0.23 at n 75-114).
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

OUT_DIR <- dir_of("11_generative_gof", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
NSIM <- TESTS$nsim_empirical
SEED_R <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 50000
GMIN_FLOOR <- GINHOM_LOWER_KM[["A"]]
BAND <- c(0.04, 0.06)
PAIR_R <- c(1, 2, 3, 5)
CONT <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5")
LITH <- c("lith_basalt", "lith_sandstone", "lith_carbonate", "lith_other")

# ---- open all inputs BEFORE sink ----
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
SITES <- read.csv(f_sites, stringsAsFactors = FALSE)
stopifnot(nrow(SITES) == 197, max(abs(X_obs_km$x - SITES$x_m / 1000), abs(X_obs_km$y - SITES$y_m / 1000)) < 1e-3)
f_grid <- file.path(REPO_DIR, "data_derived", "robust_grid_covs.rds"); stopifnot(file.exists(f_grid))
GR <- readRDS(f_grid); C250 <- GR$C250; C1000 <- GR$C1000; ZS2 <- GR$ZS2; ZS3 <- GR$ZS3
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")

ROW <- data.frame(
  date = "2026-11-12",
  item = "Clarification: robustness M1 runs - g_inhom lower bound floored at the primary 2.0 km; unwinsorised-TPI setting deferred; retention-profile reporting",
  reason = "Script 48 printed a g_inhom lower bound of 1.5 km for the eroded window (frozen bandwidth rule alone). The lower scale boundary must respect coordinate precision, so every robustness setting uses max(rule, 2.0 km), which changes only the eroded window. The frozen unwinsorised-TPI sensitivity needs a layer that was not saved (env directories hold only the winsorised z_tpi2km) and is deferred until the derivation script is available. Settings use the generic code that reproduced the recovery code exactly on the baseline pattern (script 48).",
  effect_on_inference = "The eroded-window test covers the same scales (2-10 km) as the primary test. Robustness results are reported as a retention profile (p-values, Monte Carlo SE, pair-count ratios) without multiplicity correction; non-rejections at n <= 114 are inconclusive given Stage 1 power. No interaction claim is affected (Gate 6 FAIL stands).",
  approved_by = "author-delegated (assistant decision; clarification, no frozen item changed)", stringsAsFactors = FALSE)
if (identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_49.csv"), overwrite = TRUE)
  if (!(ROW$item %in% LOG0$item)) write.csv(rbind(LOG0, ROW), f_log, row.names = FALSE)
  LOG_MSG <- "Log row appended (backup in outputs)."
} else { write.csv(ROW, file.path(OUT_DIR, "49_log_row_for_manual_append.csv"), row.names = FALSE); LOG_MSG <- "Log columns differ; row saved for manual append." }

# ---- generic M1 test code (identical logic to script 48; reproduced recovery code exactly) ----
WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
pairs_at <- function(P, r) sum(closepaircounts(P, r)) / 2
mk <- function(name, X, W, cv, dum, Lmin = 1, corrL = "isotropic", corrG = "Ripley") {
  n <- npoints(X); a <- area.owin(W); half <- PCF_STOYAN_C / sqrt(n / a); gmin <- max(Lmin, GMIN_FLOOR, ceiling(half * 2) / 2)
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
rob_run <- function(S, nsim, seed) {
  set.seed(seed)
  ft <- rob_fit(S$X, S); ob <- rob_curves(S$X, ft, S)
  lam <- exp(Reduce(`+`, Map(function(a, b) a * b, S$covW, ft$beta[-1])) + ft$beta[1])
  sL <- matrix(NA_real_, length(S$rL), nsim); sG <- matrix(NA_real_, length(S$rG), nsim)
  pr <- matrix(NA_real_, length(PAIR_R), nsim); nn <- rep(NA_real_, nsim); nd <- 0
  for (k in seq_len(nsim)) {
    res <- tryCatch({
      Xs <- rob_poly(rpoispp(lam), S); Xc <- rob_clean(Xs, S); fs <- rob_fit(Xc, S)
      list(cs = rob_curves(Xc, fs, S), n = npoints(Xc), d = npoints(Xs) - npoints(Xc), pr = vapply(PAIR_R, function(r) pairs_at(Xc, r), numeric(1)))
    }, error = function(e) NULL)
    if (!is.null(res)) { sL[, k] <- res$cs$L; sG[, k] <- res$cs$g; nn[k] <- res$n; nd <- nd + res$d; pr[, k] <- res$pr }
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  SLo <- sL[, okc, drop = FALSE]; SGo <- sG[, okc, drop = FALSE]
  csL <- GET::create_curve_set(list(r = S$rL, obs = ob$L, sim_m = SLo))
  csG <- GET::create_curve_set(list(r = S$rG, obs = ob$g, sim_m = SGo))
  p <- erl_p(csL, csG, 2)
  zv <- function(o, M) (o - rowMeans(M)) / apply(M, 1, sd)
  zr <- abs(c(zv(ob$L, SLo), zv(ob$g, SGo))); zr[!is.finite(zr)] <- NA
  rr <- c(S$rL, S$rG); ia <- which.max(zr)
  list(p = p, usable = sum(okc), n_mean = mean(nn[okc]), dropped = nd, zraw = zr[ia], zr = rr[ia],
       zcv = if (ia <= length(S$rL)) "L" else "g", pr = pr[, okc, drop = FALSE], ob = ob, SLo = SLo, SGo = SGo)
}
cv_df <- function(nm, S, r) rbind(
  data.frame(setting = nm, curve = "L_inhom_minus_r", r_km = S$rL, obs = r$ob$L, sim_mean = rowMeans(r$SLo),
             lo = apply(r$SLo, 1, quantile, 0.025), hi = apply(r$SLo, 1, quantile, 0.975)),
  data.frame(setting = nm, curve = "g_inhom", r_km = S$rG, obs = r$ob$g, sim_mean = rowMeans(r$SGo),
             lo = apply(r$SGo, 1, quantile, 0.025), hi = apply(r$SGo, 1, quantile, 0.975)))

REPORT <- file.path(OUT_DIR, "49_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("49_robustness_m1_runs | ", format(Sys.time()), "\n", sep = "")
    cat(LOG_MSG, "\n")
    cat(sprintf("nsim %d | ERL nstep 2 | seed base %d | g lower bound floor %.1f km | precision band %.2f-%.2f\n", NSIM, SEED_R, GMIN_FLOOR, BAND[1], BAND[2]))
    
    ## settings
    Wer <- erosion(W_km, 5); Xe <- X_obs_km[Wer]
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
    cat(sprintf("settings %d | sites: %s\n", length(ST), paste(sprintf("%s %d", names(ST), vapply(ST, function(S) S$n, 1L)), collapse = ", ")))
    cat("g ranges (km, lower bound): ", paste(sprintf("%s %.1f", names(ST), vapply(ST, function(S) min(S$rG), 1)), collapse = ", "), "\n", sep = "")
    
    cat("\n== Border vs isotropic L on the observed pattern (cosmetic-warning check) ==\n")
    S00 <- ST$R00_baseline; S14 <- ST$R14_edge_border_L
    o00 <- rob_curves(S00$X, rob_fit(S00$X, S00), S00); o14 <- rob_curves(S14$X, rob_fit(S14$X, S14), S14)
    ii <- match(c(2, 5, 10), S00$rL)
    print(data.frame(r_km = S00$rL[ii], L_isotropic = o00$L[ii], L_border = o14$L[ii], finite = is.finite(o14$L[ii])), digits = 3, row.names = FALSE)
    
    cat("\n== Runs ==\n")
    for (nm in names(ST)) {
      f_r <- file.path(OUT_DIR, sprintf("49_res_%s.csv", nm))
      if (file.exists(f_r)) { cat(sprintf("  %-20s found, skipped\n", nm)); next }
      S <- ST[[nm]]; k <- match(nm, names(ST)) - 1; w0 <- length(WARN$msgs); t0 <- unname(proc.time()["elapsed"])
      r <- tryCatch(rob_run(S, NSIM, SEED_R + k), error = function(e) { cat(sprintf("  %-20s FAILED: %s\n", nm, conditionMessage(e))); NULL })
      if (is.null(r)) next
      sec <- unname(proc.time()["elapsed"]) - t0; nw <- length(WARN$msgs) - w0
      po <- vapply(PAIR_R, function(x) pairs_at(S$X, x), numeric(1)); ps <- rowMeans(r$pr)
      flag <- if (is.na(r$p)) "NA" else if (r$p >= BAND[1] && r$p <= BAND[2]) "PRECISION" else if (r$p < BAND[1]) "DEFICIENCY" else if (S$n <= 114) "NO_DEFICIENCY_INCONCLUSIVE_LOW_N" else "NO_DEFICIENCY"
      row <- data.frame(setting = nm, n = S$n, df = length(S$cov), L_range = sprintf("%.1f-%.1f", min(S$rL), max(S$rL)),
                        g_range = sprintf("%.1f-%.1f", min(S$rG), max(S$rG)), nsim = NSIM, usable = r$usable, p = r$p,
                        mc_se = sqrt(r$p * (1 - r$p) / r$usable), flag = flag, sim_n_mean = r$n_mean, cleaned = r$dropped,
                        pairs_obs_2 = po[2], pairs_obs_3 = po[3], pairs_obs_5 = po[4], pairs_sim_2 = ps[2], pairs_sim_3 = ps[3], pairs_sim_5 = ps[4],
                        ratio_2 = po[2] / ps[2], ratio_3 = po[3] / ps[3], ratio_5 = po[4] / ps[4],
                        zraw_max = r$zraw, zraw_r_km = r$zr, zraw_curve = r$zcv, warnings = nw, sec = sec, seed = SEED_R + k, stringsAsFactors = FALSE)
      write.csv(cv_df(nm, S, r), file.path(OUT_DIR, sprintf("49_curves_%s.csv", nm)), row.names = FALSE)
      write.csv(row, f_r, row.names = FALSE)
      cat(sprintf("  %-20s n %3d | p %.4f (SE %.4f) | %s | pair ratio 2/3/5 km %.2f/%.2f/%.2f | usable %d | warnings %d | %.0f s [%s]\n",
                  nm, S$n, r$p, row$mc_se, flag, row$ratio_2, row$ratio_3, row$ratio_5, r$usable, nw, sec, format(Sys.time(), "%H:%M")))
    }
    
    ## summary
    fs <- file.path(OUT_DIR, sprintf("49_res_%s.csv", names(ST)))
    RES <- do.call(rbind, lapply(fs[file.exists(fs)], read.csv, stringsAsFactors = FALSE))
    write.csv(RES, file.path(OUT_DIR, "49_robustness_summary.csv"), row.names = FALSE)
    fc <- file.path(OUT_DIR, sprintf("49_curves_%s.csv", names(ST)))
    write.csv(do.call(rbind, lapply(fc[file.exists(fc)], read.csv, stringsAsFactors = FALSE)), file.path(OUT_DIR, "49_curves_all.csv"), row.names = FALSE)
    cat("\n== Summary ==\n")
    print(RES[, c("setting", "n", "p", "mc_se", "flag", "ratio_2", "ratio_3", "ratio_5", "usable", "warnings", "sec")], digits = 3, row.names = FALSE)
    cat("\nsettings flagged for the precision rule (script 49b):", paste(RES$setting[RES$flag == "PRECISION"], collapse = " "), "\n")
    cat("R00 control p", RES$p[RES$setting == "R00_baseline"], "vs script 42 p 0.0310 (different seed)\n")
    cat("\n== Warning tally (main session) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "49_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")