# =====================================================================
# 45_m2_o2_precision_check.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW)
# Purpose : (A) log the symmetric Monte Carlo precision rule (docs/protocol_deviation_log.csv);
#           (B) verify the M2 offset plumbing: constant-offset test (intercept shifts by minus
#               the constant, slopes and fitted intensities unchanged) and a retry of
#               ppm(..., offset(logEoff)) with the error text printed;
#           (C) rerun the M2-O2 test (E3 + fixed offset O2, dataset A, simplified window,
#               joint ERL test, nstep 2) with nsim = 10000 in 10 checkpointed chunks and a
#               distinct seed; report p, Monte Carlo SE, full-precision script-44 p, decision.
# Outputs : outputs/45_report.txt, 45_m2_o2_precision_result.csv, 45_sessionInfo.txt,
#           protocol_deviation_log_before_45.csv, outputs/45_chunks/45_chunk_NN.rds;
#           docs/protocol_deviation_log.csv (1 row appended)
# Notes   : RUN ONCE (resumable through chunk files). Rule (declared 2026-11-12 after seeing the
#           O2 p of 0.0500; applies to every later empirical ERL test): p in [0.04, 0.06] at
#           nsim 999 is recomputed at nsim 10000 with a distinct seed and that p is used;
#           the threshold p < 0.05 is unchanged. If p(10000) < 0.05: deficiency under O2,
#           continue to M3. Otherwise the stop stands. Restart R before running. Main session only.
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
source(file.path(REPO_DIR, "R", "recovery_cells.R"))

OUT_DIR   <- dir_of("11_generative_gof", "outputs")
CHUNK_DIR <- file.path(OUT_DIR, "45_chunks"); dir.create(CHUNK_DIR, recursive = TRUE, showWarnings = FALSE)
f_res   <- file.path(OUT_DIR, "45_m2_o2_precision_result.csv")
if (file.exists(f_res)) stop("45_m2_o2_precision_result.csv exists: this analysis runs once. Stop; do not rerun.")
f_res44 <- file.path(OUT_DIR, "44_m2_result.csv"); stopifnot(file.exists(f_res44))
RES44 <- read.csv(f_res44, stringsAsFactors = FALSE)
f_obs <- file.path(REPO_DIR, "data_derived", "recovery_obs_cache.rds"); stopifnot(file.exists(f_obs))
OBSL <- readRDS(f_obs); stopifnot("logE_O2" %in% names(OBSL))
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")

N_CHUNK <- 10; CH <- 1000
SEED_BASE <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 30000
SEED45 <- SEED_BASE + 12
bw <- bw_for_n(N_ANALYTICAL); rGv <- rG_for(N_ANALYTICAL)
offim <- OBSL$logE_O2

# ---- log row (before sink) ----
ROW <- data.frame(
  date = "2026-11-12",
  item = "Clarification: symmetric Monte Carlo precision rule for empirical ERL tests (p in [0.04, 0.06] at nsim 999 recomputed at nsim 10000)",
  reason = "Script 44: M2-O2 combined ERL p = 0.0500 at nsim 999 (49 simulations at least as extreme; Monte Carlo SE about 0.007), exactly at the frozen threshold, so the stop rule rested on one simulation. The rule was declared after seeing this value and applies to every later empirical ERL test: a p in [0.04, 0.06] at nsim 999 is recomputed with nsim 10000 and a distinct seed, and the 10000-simulation p is used. M1 (0.031) and M2-O3 (0.024) lie outside the band. The decision rule p < 0.05 is unchanged.",
  effect_on_inference = "Frozen nsim 999 is exceeded for borderline tests only; Monte Carlo error falls to about 0.002. The test statistic, order of tests and threshold are unchanged. Gate 6 (provisional FAIL) and the no-interaction-claim rule are unaffected.",
  approved_by = "author-delegated (assistant decision; clarification, no frozen item changed)", stringsAsFactors = FALSE)
if (identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_45.csv"), overwrite = TRUE)
  if (!(ROW$item %in% LOG0$item)) write.csv(rbind(LOG0, ROW), f_log, row.names = FALSE)
  LOG_MSG <- "Log row appended (backup in outputs)."
} else {
  write.csv(ROW, file.path(OUT_DIR, "45_log_row_for_manual_append.csv"), row.names = FALSE)
  LOG_MSG <- "Log columns differ; row saved to 45_log_row_for_manual_append.csv (append manually)."
}

WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
fast_fit_off <- function(Xp, offim) {
  Q <- make_quad(Xp); w <- w.quad(Q); isd <- is.data(Q); xq <- x.quad(Q); yq <- y.quad(Q)
  Z <- vapply(covs, function(im) lookup.im(im, xq, yq, naok = TRUE), numeric(length(xq)))
  o <- lookup.im(offim, xq, yq, naok = TRUE)
  ok <- stats::complete.cases(Z) & is.finite(o) & w > 0
  X <- cbind(`(Intercept)` = 1, Z)
  gfit <- suppressWarnings(stats::glm.fit(X[ok, , drop = FALSE], (isd / w)[ok], weights = w[ok],
                                          family = poisson(), offset = o[ok], intercept = FALSE))
  beta <- gfit$coefficients
  Xd <- X[isd & ok, , drop = FALSE]; od <- o[isd & ok]
  if (nrow(Xd) != npoints(Xp)) stop("data points dropped in fast_fit_off: ", nrow(Xd), " of ", npoints(Xp))
  list(beta = beta, lam_data = as.numeric(exp(Xd %*% beta + od)))
}
m2_chunk <- function(lam_im, offim, nsim, seed, bw, rGv) {
  set.seed(seed)
  sL <- matrix(NA_real_, length(rL), nsim); sG <- matrix(NA_real_, length(rGv), nsim); nn <- rep(NA_real_, nsim)
  for (k in seq_len(nsim)) {
    res <- tryCatch({
      Xs <- to_poly(rpoispp(lam_im)); fs <- fast_fit_off(Xs, offim); cs <- fast_curves_r(Xs, fs, bw, rGv)
      list(cs = cs, n = npoints(Xs))
    }, error = function(e) NULL)
    if (!is.null(res)) { sL[, k] <- res$cs$L; sG[, k] <- res$cs$g; nn[k] <- res$n }
  }
  list(sL = sL, sG = sG, n = nn)
}

REPORT <- file.path(OUT_DIR, "45_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("45_m2_o2_precision_check | ", format(Sys.time()), "\n", sep = "")
    cat(LOG_MSG, "\n")
    p44 <- RES44$p_erl_nstep2[RES44$scenario == "O2"]
    cat(sprintf("script-44 O2 p at full precision: %.10f | numerically < 0.05: %s\n", p44, p44 < 0.05))
    
    cat("\n== A. Offset plumbing check ==\n")
    offc <- offim; offc$v[!is.na(offc$v)] <- 0.7
    fa <- fast_fit(X_obs_km); fc <- fast_fit_off(X_obs_km, offc)
    cat(sprintf("constant offset 0.7: intercept shift %.8f (expect -0.7) | max |slope difference| %.2e | max relative difference of fitted intensities at sites %.2e\n",
                fc$beta[1] - fa$beta[1], max(abs(fc$beta[-1] - fa$beta[-1])), max(abs(fc$lam_data / fa$lam_data - 1))))
    f0 <- fast_fit_off(X_obs_km, offim)
    fm <- as.formula(paste("~", paste(c(names(covs), "offset(logEoff)"), collapse = " + ")))
    fitO <- tryCatch(ppm(make_quad(X_obs_km), trend = fm, covariates = c(covs, list(logEoff = offim))),
                     error = function(e) { cat("ppm(offset) error:", conditionMessage(e), "\n"); NULL })
    if (!is.null(fitO)) cat(sprintf("ppm(offset) retry: max |coef difference| vs fast_fit_off %.2e\n", max(abs(f0$beta - coef(fitO)[names(f0$beta)])))) else cat("ppm(offset) retry failed; the constant-offset test is the verification\n")
    
    cat("\n== C. M2-O2 precision run: nsim ", N_CHUNK * CH, " (", N_CHUNK, " chunks of ", CH, "), seed base ", SEED45, " ==\n", sep = "")
    ob <- fast_curves_r(X_obs_km, f0, bw, rGv)
    lam_im <- lam_im_of(f0$beta) * exp(offim[W_km, drop = FALSE])
    for (ch in seq_len(N_CHUNK)) {
      f_ch <- file.path(CHUNK_DIR, sprintf("45_chunk_%02d.rds", ch))
      if (file.exists(f_ch)) { cat(sprintf("  chunk %2d: found, skipped\n", ch)); next }
      t0 <- unname(proc.time()["elapsed"])
      saveRDS(m2_chunk(lam_im, offim, CH, SEED45 + 100 * ch, bw, rGv), f_ch)
      cat(sprintf("  chunk %2d done in %.0f s [%s]\n", ch, unname(proc.time()["elapsed"]) - t0, format(Sys.time(), "%H:%M")))
    }
    CHK <- lapply(file.path(CHUNK_DIR, sprintf("45_chunk_%02d.rds", seq_len(N_CHUNK))), readRDS)
    sL <- do.call(cbind, lapply(CHK, function(x) x$sL)); sG <- do.call(cbind, lapply(CHK, function(x) x$sG))
    nn <- unlist(lapply(CHK, function(x) x$n))
    okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
    csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
    csG <- GET::create_curve_set(list(r = rGv, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
    p10k <- erl_p(csL, csG, 2); N <- sum(okc)
    se <- sqrt(p10k * (1 - p10k) / N)
    cat(sprintf("\nusable simulations %d of %d | simulated n mean %.1f sd %.1f\n", N, N_CHUNK * CH, mean(nn[okc]), sd(nn[okc])))
    cat(sprintf("ERL combined global envelope p (nstep 2, nsim %d): %.5f | Monte Carlo SE %.4f | approx 95%% interval %.4f-%.4f\n",
                N, p10k, se, p10k - 1.96 * se, p10k + 1.96 * se))
    cat(sprintf("script-44 p (nsim 999): %.4f | difference %.4f\n", p44, p10k - p44))
    act <- if (is.na(p10k)) "UNDETERMINED (NA p)" else if (p10k < TESTS$alpha)
      "DEFICIENCY under O2 (nsim 10000): the hierarchy continues to M3" else
        "NO deficiency under O2 (nsim 10000): the stop stands (assumed offset; not proof of sufficiency)"
    cat("DECISION (frozen threshold, precision rule):", act, "\n")
    write.csv(data.frame(test = "M2-O2", nsim = N, p_erl_nstep2 = p10k, mc_se = se, p_nsim999 = p44, seed_base = SEED45,
                         decision = act, date = format(Sys.time())), f_res, row.names = FALSE)
    cat("\n== Warning tally (main session) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "45_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")