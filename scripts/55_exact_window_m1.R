# =====================================================================
# 55_exact_window_m1.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW)
# Purpose : full exact-window M1 test (frozen primary window: exact boundary polygon, 18,926
#           vertices) after the reduced check (script 54) triggered the pre-set escalation rule.
#           E3 Poisson, O1, dataset A, nominal coordinates, joint ERL L_inhom + g_inhom, nstep 2,
#           nsim 999 (3 chunks of 333, distinct seeds), main session, generic code from
#           R/robust_generic.R (identical to scripts 49 / 49b). Bandwidth held at the simplified-
#           window value to isolate the window effect.
# Outputs : outputs/55_report.txt, 55_exact_window_m1_result.csv, 55_chunk_NN.rds,
#           55_sessionInfo.txt, protocol_deviation_log_before_55.csv;
#           docs/protocol_deviation_log.csv (1 row appended)
# Notes   : RUN ONCE (resumable by chunk). Reading rules (set before running): p < 0.04 deficiency
#           retained under the exact window; p in [0.04, 0.06] precision rule (nsim 10000, script 55b);
#           p > 0.06 not retained. Compared with script 42 (simplified window, same nsim 999: 0.0310).
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
GENERIC_FILE <- file.path(REPO_DIR, "R", "robust_generic.R"); stopifnot(file.exists(GENERIC_FILE))
source(GENERIC_FILE)

OUT_DIR <- dir_of("11_generative_gof", "outputs")
f_res <- file.path(OUT_DIR, "55_exact_window_m1_result.csv")
if (file.exists(f_res)) stop("55_exact_window_m1_result.csv exists: this analysis runs once. Stop; do not rerun.")
f_c0 <- file.path(REPO_DIR, "data_derived", "recovery_cache.rds"); stopifnot(file.exists(f_c0))
f_r42 <- file.path(OUT_DIR, "42_m1_result.csv"); stopifnot(file.exists(f_r42))
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
R42 <- read.csv(f_r42, stringsAsFactors = FALSE); LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")
W_ex <- readRDS(f_c0)$W_km
nv <- function(W) sum(vapply(W$bdry, function(b) length(b$x), 1L))
stopifnot(W_ex$type == "polygonal", nv(W_ex) > 5000)
SEED_X <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 120000
NCH <- 3; CH <- 333

ROW <- data.frame(
  date = "2026-11-13",
  item = "Clarification: exact-window sensitivity run as a full M1 test at nsim 999 after the reduced check triggered the pre-set escalation rule",
  reason = "Script 54 (reduced design): observed curves on the exact window differ from the simplified-window curves by at most 0.54% of the null sd and the ERL p with exact-window observed curves is identical (0.0310); simulated-pattern comparisons exceeded the pre-set 0.05 sd threshold only for patterns that lost one point in the 28.9 km2 symmetric-difference sliver (2 of 12 above threshold, 3 of 12 with a lost point), which is a point-loss effect, not a window effect on the statistic. The pre-set escalation rule (full exact-window M1) was followed. Measured cost was about 1 s per pattern, not the 57 s of script 21, so the full test runs at the frozen nsim 999 in the main session instead of nsim 199 on workers (nsim 199 p-values are not comparable with nsim 999 p-values for this pattern, scripts 51-52).",
  effect_on_inference = "Exact-window result is reported with the same classification rule as the other robustness settings. The simplified window remains the window of all simulation-based tests.",
  approved_by = "author-delegated (assistant decision; clarification, no frozen item changed)", stringsAsFactors = FALSE)
if (identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_55.csv"), overwrite = TRUE)
  if (!(ROW$item %in% LOG0$item)) write.csv(rbind(LOG0, ROW), f_log, row.names = FALSE)
  LOG_MSG <- "Log row appended (backup in outputs)."
} else { write.csv(ROW, file.path(OUT_DIR, "55_log_row_for_manual_append.csv"), row.names = FALSE); LOG_MSG <- "Log columns differ; row saved for manual append." }

REPORT <- file.path(OUT_DIR, "55_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("55_exact_window_m1 | ", format(Sys.time()), "\n", sep = "")
  cat(LOG_MSG, "\n")
  cat(sprintf("exact window: %d vertices, area %.2f km2 | nsim %d (%d chunks of %d) | seed base %d\n", nv(W_ex), area.owin(W_ex), NCH * CH, NCH, CH, SEED_X))
  X_ex <- X_obs_km[W_ex]; cat("sites inside the exact window:", npoints(X_ex), "of 197\n"); stopifnot(npoints(X_ex) == 197)
  S_s <- mk("SIMPL", X_obs_km, W_km, covs, DUMMY_FIXED)
  t0 <- unname(proc.time()["elapsed"])
  S_x <- mk("EXACT", X_ex, W_ex, covs, default.dummy(X_ex, eps = QUAD_KM)); S_x$bw <- S_s$bw
  cat(sprintf("exact-window setting built in %.0f s | L range %.1f-%.1f, g range %.1f-%.1f (simplified: %.1f-%.1f, %.1f-%.1f)\n", unname(proc.time()["elapsed"]) - t0,
              min(S_x$rL), max(S_x$rL), min(S_x$rG), max(S_x$rG), min(S_s$rL), max(S_s$rL), min(S_s$rG), max(S_s$rG)))
  stopifnot(identical(S_x$rL, S_s$rL), identical(S_x$rG, S_s$rG))
  ft <- rob_fit(X_ex, S_x); ob <- rob_curves(X_ex, ft, S_x)
  
  cat("\n== Chunks ==\n")
  for (ch in seq_len(NCH)) {
    f_c <- file.path(OUT_DIR, sprintf("55_chunk_%02d.rds", ch))
    if (file.exists(f_c)) { cat(sprintf("  chunk %d: found, skipped\n", ch)); next }
    t1 <- unname(proc.time()["elapsed"])
    saveRDS(chunk_run(S_x, CH, SEED_X + ch), f_c)
    cat(sprintf("  chunk %d done in %.0f s [%s]\n", ch, unname(proc.time()["elapsed"]) - t1, format(Sys.time(), "%H:%M")))
  }
  CK <- lapply(file.path(OUT_DIR, sprintf("55_chunk_%02d.rds", seq_len(NCH))), readRDS)
  sL <- do.call(cbind, lapply(CK, function(x) x$sL)); sG <- do.call(cbind, lapply(CK, function(x) x$sG))
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = S_x$rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = S_x$rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  p <- erl_p(csL, csG, 2); N <- sum(okc); se <- sqrt(p * (1 - p) / N)
  flag <- if (is.na(p)) "NA" else if (p >= 0.04 && p <= 0.06) "PRECISION" else if (p < 0.04) "DEFICIENCY retained (exact window)" else "NOT retained (exact window)"
  cat(sprintf("\nusable simulations %d of %d\n", N, NCH * CH))
  cat(sprintf("ERL combined global envelope p (nstep 2, nsim %d): %.4f (SE %.4f) | script 42 simplified window nsim 999: %.4f | difference %.4f\n", N, p, se, R42$p_erl_nstep2, p - R42$p_erl_nstep2))
  zv <- function(o, M) (o - rowMeans(M)) / apply(M, 1, sd)
  zr <- abs(c(zv(ob$L, sL[, okc, drop = FALSE]), zv(ob$g, sG[, okc, drop = FALSE]))); zr[!is.finite(zr)] <- NA
  rr <- c(S_x$rL, S_x$rG); ia <- which.max(zr)
  cat(sprintf("raw max |z| %.2f at r = %.1f km on the %s curve (simplified window: 6.58 at 2.5 km on g_inhom, floored)\n", zr[ia], rr[ia], if (ia <= length(S_x$rL)) "L_inhom" else "g_inhom"))
  cat("classification:", flag, "\n")
  write.csv(data.frame(setting = "R18_exact_window", n = 197, nsim = N, p = p, mc_se = se, flag = flag, p_script42 = R42$p_erl_nstep2,
                       zraw_max = zr[ia], zraw_r_km = rr[ia], seed_base = SEED_X, date = format(Sys.time())), f_res, row.names = FALSE)
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "55_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")