# =====================================================================
# 58_anisotropy_diagnostic.R
# Stage   : 13_sensitivity
# Folder  : 13_sensitivity/scripts/   (NEW; creates 13_sensitivity/outputs)
# Purpose : directional second-order diagnostic (Research_Design: establish whether residual
#           dependence is directional before any anisotropic model; none is fitted here).
#           Statistic per distance band (2-5, 5-10, 2-10 km): Rayleigh statistic on doubled pair
#           angles, S = |sum exp(2 i theta)|^2 / N_pairs (axial: theta and theta + 180 deg equal),
#           dominant orientation phi (degrees from east, 0-180); plus six 30-degree sector shares
#           (2-10 km). Observed vs two nulls: M3 (E3 + LGCP, conditional rLGCP eps 0.5 + rpoint,
#           499 patterns, as in script 46) and M1 (E3 Poisson, 999 patterns, as in script 42).
#           One-sided Monte Carlo p = (1 + #(S_sim >= S_obs)) / (N + 1).
# Outputs : outputs/58_report.txt, 58_anisotropy_bands.csv, 58_anisotropy_sectors.csv,
#           58_sessionInfo.txt, protocol_deviation_log_before_58.csv;
#           docs/protocol_deviation_log.csv (1 row)
# Notes   : Reading rules (set before running): Bonferroni alpha 0.05 / 6 = 0.0083 over 3 bands x
#           2 nulls. ANISOTROPY INDICATED if p(M3) < 0.0083 in any band; SUGGESTED if only p(M1) is;
#           otherwise isotropic primary retained. No anisotropic model is fitted; a direction is not
#           read as a corridor or behavioural axis. Restart R before running. Main session only.
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

OUT_DIR <- dir_of("13_sensitivity", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
NSIM_M1 <- 999; NSIM_M3 <- 499; CHK <- 100
SEED_A <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 130000
BANDS <- list(B1 = c(2, 5), B2 = c(5, 10), B3 = c(2, 10)); NSEC <- 6; ALPHA_B <- 0.05 / 6

f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")
ROW <- data.frame(
  date = "2026-11-14",
  item = "Clarification: anisotropy diagnostic design (Rayleigh statistic on doubled pair angles against M1 and M3 nulls)",
  reason = "Research_Design requires establishing whether residual dependence is directional before any anisotropic model and forbids fitting an anisotropic Gibbs model automatically. Statistic chosen before any result: Rayleigh statistic S = |sum exp(2 i theta)|^2 / N on pair orientations in distance bands 2-5, 5-10 and 2-10 km, Monte Carlo p against M3 (499 conditional LGCP patterns) and M1 (999 Poisson patterns), Bonferroni alpha 0.0083 over 6 tests; six 30-degree sector shares for 2-10 km.",
  effect_on_inference = "Diagnostic only. The isotropic primary analysis is retained unless anisotropy is indicated; no anisotropic model is fitted and no direction is interpreted as a corridor or behavioural axis. Gate 6 unaffected.",
  approved_by = "author-delegated (assistant decision; clarification, no frozen item changed)", stringsAsFactors = FALSE)
if (identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_58.csv"), overwrite = TRUE)
  if (!(ROW$item %in% LOG0$item)) write.csv(rbind(LOG0, ROW), f_log, row.names = FALSE)
  LOG_MSG <- "Log row appended (backup in outputs)."
} else { write.csv(ROW, file.path(OUT_DIR, "58_log_row_for_manual_append.csv"), row.names = FALSE); LOG_MSG <- "Log columns differ; row saved for manual append." }

NMS <- c(unlist(lapply(names(BANDS), function(b) paste0(c("n_", "S_", "phi_"), b))), paste0("sec", seq_len(NSEC)))
dir_stat <- function(P) {
  cp <- closepairs(P, 10, twice = FALSE)
  d <- cp$d; th <- atan2(cp$dy, cp$dx); res <- numeric(0)
  for (b in BANDS) {
    sel <- d >= b[1] & d < b[2]; n <- sum(sel)
    z <- if (n > 0) sum(exp(2i * th[sel])) else 0 + 0i
    res <- c(res, n, if (n > 0) Mod(z)^2 / n else NA_real_, if (n > 0) (Arg(z) / 2 * 180 / pi) %% 180 else NA_real_)
  }
  ang <- (th %% pi) * 180 / pi; sel <- d >= 2 & d < 10
  c(res, tabulate(pmin(floor(ang[sel] / (180 / NSEC)) + 1, NSEC), nbins = NSEC))
}
mc_p <- function(o, s) (1 + sum(s >= o, na.rm = TRUE)) / (sum(!is.na(s)) + 1)

REPORT <- file.path(OUT_DIR, "58_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("58_anisotropy_diagnostic | ", format(Sys.time()), "\n", sep = "")
  cat(LOG_MSG, "\n")
  cat(sprintf("bands %s | nsim M1 %d, M3 %d | seed base %d | Bonferroni alpha %.4f\n", paste(sprintf("%g-%g km", sapply(BANDS, `[`, 1), sapply(BANDS, `[`, 2)), collapse = ", "), NSIM_M1, NSIM_M3, SEED_A, ALPHA_B))
  
  cat("\n== 0. Observed pattern ==\n")
  cpo <- closepairs(X_obs_km, 10, twice = FALSE)
  cat("observed pairs within 1, 2, 3, 5 km:", paste(c(sum(cpo$d < 1), sum(cpo$d < 2), sum(cpo$d < 3), sum(cpo$d < 5)), collapse = " "), "(expected 7 46 106 227)\n")
  stopifnot(identical(as.numeric(c(sum(cpo$d < 1), sum(cpo$d < 2), sum(cpo$d < 3), sum(cpo$d < 5))), c(7, 46, 106, 227)))
  obs <- dir_stat(X_obs_km); names(obs) <- NMS
  
  cat("\n== 1. M1 null (", NSIM_M1, " patterns) ==\n", sep = "")
  f0 <- fast_fit(X_obs_km); lam <- lam_im_of(f0$beta); set.seed(SEED_A + 1); t0 <- unname(proc.time()["elapsed"])
  M1 <- t(vapply(seq_len(NSIM_M1), function(k) dir_stat(to_poly(rpoispp(lam))), numeric(length(NMS)))); colnames(M1) <- NMS
  cat(sprintf("done in %.0f s\n", unname(proc.time()["elapsed"]) - t0))
  
  cat("\n== 2. M3 null (", NSIM_M3, " conditional patterns) ==\n", sep = "")
  fit <- kppm(X_obs_km, trend = formula_E3, clusters = "LGCP", model = "exponential", covariates = covs, method = "clik2")
  s2 <- unname(fit$par["sigma2"]); al <- unname(fit$par["alpha"]); eta <- trend_im(coef(fit)); n <- npoints(X_obs_km)
  cat(sprintf("kppm sigma2_hat %.4f alpha_hat %.3f km (script 46: 1.2298, 6.048)\n", s2, al))
  M3 <- matrix(NA_real_, NSIM_M3, length(NMS)); colnames(M3) <- NMS; k0 <- 0; t0 <- unname(proc.time()["elapsed"])
  while (k0 < NSIM_M3) {
    m <- min(CHK, NSIM_M3 - k0); set.seed(SEED_A + 100 + k0)
    Sl <- as_list(rLGCP(model = "exp", mu = eta - s2 / 2, param = list(var = s2, scale = al), win = W_km, eps = EPS_LGCP,
                        saveLambda = TRUE, nsim = m, drop = TRUE))
    for (j in seq_len(m)) M3[k0 + j, ] <- dir_stat(to_poly(rpoint(n, f = attr(Sl[[j]], "Lambda"), win = W_km, forcewin = TRUE)))
    k0 <- k0 + m; cat(sprintf("  %d of %d [%.0f s]\n", k0, NSIM_M3, unname(proc.time()["elapsed"]) - t0))
  }
  
  cat("\n== 3. Bands ==\n")
  BD <- do.call(rbind, lapply(names(BANDS), function(b) {
    nb <- paste0("n_", b); sb <- paste0("S_", b); pb <- paste0("phi_", b)
    data.frame(band = sprintf("%g-%g km", BANDS[[b]][1], BANDS[[b]][2]), n_obs = obs[[nb]], n_M1 = mean(M1[, nb]), n_M3 = mean(M3[, nb]),
               S_obs = obs[[sb]], S_M1_median = median(M1[, sb], na.rm = TRUE), S_M1_q95 = quantile(M1[, sb], .95, na.rm = TRUE), p_M1 = mc_p(obs[[sb]], M1[, sb]),
               S_M3_median = median(M3[, sb], na.rm = TRUE), S_M3_q95 = quantile(M3[, sb], .95, na.rm = TRUE), p_M3 = mc_p(obs[[sb]], M3[, sb]), phi_obs_deg = obs[[pb]])
  }))
  print(BD, digits = 3, row.names = FALSE); write.csv(BD, file.path(OUT_DIR, "58_anisotropy_bands.csv"), row.names = FALSE)
  
  cat("\n== 4. Sector shares, 2-10 km pairs (orientation mod 180 deg; 0 = east-west) ==\n")
  sc <- paste0("sec", seq_len(NSEC)); lab <- sprintf("%d-%d", seq(0, 150, by = 30), seq(30, 180, by = 30))
  sh <- function(M) M[, sc, drop = FALSE] / rowSums(M[, sc, drop = FALSE])
  s1 <- sh(M1); s3 <- sh(M3); so <- obs[sc] / sum(obs[sc])
  SEC <- data.frame(sector_deg = lab, obs_share = as.numeric(so),
                    M3_mean = colMeans(s3), M3_lo = apply(s3, 2, quantile, .025), M3_hi = apply(s3, 2, quantile, .975),
                    M1_mean = colMeans(s1), M1_lo = apply(s1, 2, quantile, .025), M1_hi = apply(s1, 2, quantile, .975))
  print(SEC, digits = 3, row.names = FALSE); write.csv(SEC, file.path(OUT_DIR, "58_anisotropy_sectors.csv"), row.names = FALSE)
  
  cat("\n== Verdict (pre-set rule) ==\n")
  i3 <- BD$band[BD$p_M3 < ALPHA_B]; i1 <- BD$band[BD$p_M1 < ALPHA_B]
  cat("bands with p(M3) < ", sprintf("%.4f", ALPHA_B), ": ", if (length(i3)) paste(i3, collapse = ", ") else "none", "\n", sep = "")
  cat("bands with p(M1) < ", sprintf("%.4f", ALPHA_B), ": ", if (length(i1)) paste(i1, collapse = ", ") else "none", "\n", sep = "")
  cat("verdict:", if (length(i3)) "ANISOTROPY INDICATED (no anisotropic model fitted; direction not interpreted)" else
    if (length(i1)) "ANISOTROPY SUGGESTED (only against M1); isotropic primary retained" else "no directional structure detected; isotropic primary retained", "\n")
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "58_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")