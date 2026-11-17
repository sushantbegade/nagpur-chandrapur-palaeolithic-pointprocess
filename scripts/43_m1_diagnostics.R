# =====================================================================
# 43_m1_diagnostics.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW; does not touch script-42 outputs)
# Purpose : EXPLORATORY diagnostics of the frozen M1 test on the archaeological pattern.
#           (A) deterministic rerun of the script-42 simulations (same seed): p must
#               reproduce the 42 value; (B) pointwise tails: observed vs null quantiles,
#               one-sided Monte Carlo p; (C) single-curve-set global tests (L only, g only)
#               with erl / rank / qdir / st, and the combined nstep = 1 error text;
#               (D) null-tail shape of L_inhom and g_inhom (skewness, maximum vs observed).
# Outputs : outputs/43_report.txt, 43_m1_diag.csv, 43_m1_sims.rds, 43_sessionInfo.txt
# Notes   : Diagnostic of test SENSITIVITY only. The primary result is the frozen script-42
#           p (nstep 2); nothing here replaces it or is used to choose a statistic.
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
source(file.path(REPO_DIR, "R", "recovery_tests.R"))
source(file.path(REPO_DIR, "R", "recovery_cells.R"))

OUT_DIR <- dir_of("11_generative_gof", "outputs")
f_res <- file.path(OUT_DIR, "42_m1_result.csv"); stopifnot(file.exists(f_res))
RES42 <- read.csv(f_res, stringsAsFactors = FALSE)

NSIM   <- TESTS$nsim_empirical
SEED_E <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 30000 + 1   # identical to script 42
bw  <- bw_for_n(N_ANALYTICAL); rGv <- rG_for(N_ANALYTICAL)
R_L <- c(2, 3, 5, 10); R_G <- c(2, 3, 5, 7)
PAIR_R <- c(1, 2, 3, 5)

WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
pairs_at <- function(P, r) sum(closepaircounts(P, r)) / 2
emp_m1 <- function(Xobs, nsim, bw, rGv, seed) {
  set.seed(seed)
  f0 <- fast_fit(Xobs); ob <- fast_curves_r(Xobs, f0, bw, rGv)
  lam_im <- lam_im_of(f0$beta)
  sL <- matrix(NA_real_, length(rL), nsim); sG <- matrix(NA_real_, length(rGv), nsim)
  nn <- rep(NA_real_, nsim); pr <- matrix(NA_real_, length(PAIR_R), nsim)
  for (k in seq_len(nsim)) {
    res <- tryCatch({
      Xs <- to_poly(rpoispp(lam_im)); fs <- fast_fit(Xs); cs <- fast_curves_r(Xs, fs, bw, rGv)
      list(cs = cs, n = npoints(Xs), pr = vapply(PAIR_R, function(r) pairs_at(Xs, r), numeric(1)))
    }, error = function(e) NULL)
    if (!is.null(res)) { sL[, k] <- res$cs$L; sG[, k] <- res$cs$g; nn[k] <- res$n; pr[, k] <- res$pr }
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  list(ob = ob, sL = sL[, okc, drop = FALSE], sG = sG[, okc, drop = FALSE], usable = sum(okc))
}
skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4
tail_row <- function(curve, r, rgrid, obs, S) {
  i <- match(round(r, 6), round(rgrid, 6)); x <- S[i, ]; o <- obs[i]
  data.frame(curve = curve, r_km = r, obs = o, sim_median = median(x), q95 = quantile(x, .95), q99 = quantile(x, .99),
             q999 = quantile(x, .999), sim_max = max(x), n_sims_ge_obs = sum(x >= o), mc_p_upper = (1 + sum(x >= o)) / (length(x) + 1),
             skew = skew(x), kurt = kurt(x), z_obs = (o - mean(x)) / sd(x), z_max_sim = (max(x) - mean(x)) / sd(x))
}
gtest <- function(cs, type, nstep = NULL) {
  a <- list(curve_sets = cs, type = type, alpha = 0.05)
  if (!is.null(nstep)) a$nstep <- nstep
  tryCatch(as.numeric(attr(do.call(GET::global_envelope_test, a), "p")), error = function(e) paste("ERROR:", conditionMessage(e)))
}

REPORT <- file.path(OUT_DIR, "43_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("43_m1_diagnostics (EXPLORATORY) | ", format(Sys.time()), "\n", sep = "")
    cat("primary result (script 42, frozen): p =", RES42$p_erl_nstep2, "| nsim", RES42$nsim, "| seed", RES42$seed, "\n")
    
    cat("\n== A. Deterministic rerun of the M1 simulations ==\n")
    t0 <- unname(proc.time()["elapsed"])
    R <- emp_m1(X_obs_km, NSIM, bw, rGv, SEED_E)
    cat(sprintf("done in %.0f s | usable %d of %d\n", unname(proc.time()["elapsed"]) - t0, R$usable, NSIM))
    csL <- GET::create_curve_set(list(r = rL, obs = R$ob$L, sim_m = R$sL))
    csG <- GET::create_curve_set(list(r = rGv, obs = R$ob$g, sim_m = R$sG))
    p_re <- erl_p(csL, csG, 2)
    cat(sprintf("combined ERL p (nstep 2) rerun %.4f | script 42: %.4f | reproduces: %s\n", p_re, RES42$p_erl_nstep2,
                isTRUE(abs(p_re - RES42$p_erl_nstep2) < 1e-9)))
    saveRDS(list(rL = rL, rG = rGv, ob = R$ob, sL = R$sL, sG = R$sG), file.path(OUT_DIR, "43_m1_sims.rds"))
    
    cat("\n== B/D. Pointwise tails and null shape (one-sided Monte Carlo p; skew/kurt of the null at that r) ==\n")
    TR <- rbind(do.call(rbind, lapply(R_L, function(r) tail_row("L_inhom_minus_r", r, rL, R$ob$L, R$sL))),
                do.call(rbind, lapply(R_G, function(r) tail_row("g_inhom", r, rGv, R$ob$g, R$sG))))
    print(TR, digits = 3, row.names = FALSE)
    write.csv(TR, file.path(OUT_DIR, "43_m1_diag.csv"), row.names = FALSE)
    
    cat("\n== C. Global tests: sensitivity of the construction (diagnostic only; primary stays the frozen p) ==\n")
    cat("combined (L + g), erl, nstep 2:", gtest(list(csL, csG), "erl", 2), "\n")
    cat("combined (L + g), erl, nstep 1:", gtest(list(csL, csG), "erl", 1), "\n")
    for (ty in c("erl", "rank", "qdir", "st")) {
      cat(sprintf("type %-5s | L only p = %s | g only p = %s\n", ty, format(gtest(csL, ty)), format(gtest(csG, ty))))
    }
    cat("(p resolution with 999 simulations is 0.001)\n")
    
    cat("\n== Warning tally (main session) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "43_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")