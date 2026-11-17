# =====================================================================
# 42_empirical_m1.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW; creates 11_generative_gof/outputs)
# Purpose : FIRST FIT TO THE ARCHAEOLOGICAL DATA. Frozen M1 test: environmental
#           Poisson (E3, O1), dataset A (197 sites), 100 m simplified window,
#           nominal coordinates, trend-only E3 refit per pattern, joint ERL test
#           of L_inhom (1-10 km) and g_inhom (2-10 km), isotropic edge correction,
#           nsim = 999 (unconditional rpoispp), nstep 2. Same functions as recovery.
#           Also: E3 coefficient table (ppm, with SE), observed vs simulated curves
#           and pair counts (descriptive).
# Outputs : outputs/42_report.txt, 42_m1_result.csv, 42_m1_curves.csv, 42_m1_coef.csv,
#           42_sessionInfo.txt
# Notes   : RUN ONCE. Stops if 42_m1_result.csv exists. Every setting is frozen;
#           nothing is chosen after seeing output. Reading rules (set before running):
#           deficiency if ERL p < 0.05 (frozen). p >= 0.05 stops the hierarchy but is NOT
#           evidence of environmental sufficiency (Stage 1: M1 rejected 38% of strongly
#           clustered LGCP patterns and about 2% of near-null patterns). p < 0.05 means
#           continue to M2 (script 43). Curve bands are pointwise, display only.
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
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
f_res <- file.path(OUT_DIR, "42_m1_result.csv")
if (file.exists(f_res)) stop("42_m1_result.csv exists: this analysis runs once. Stop; do not rerun.")

NSIM  <- TESTS$nsim_empirical
N_A   <- N_ANALYTICAL
SEED_E <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 30000 + 1
PAIR_R <- c(1, 2, 3, 5)
bw  <- bw_for_n(N_A)
rGv <- rG_for(N_A)

# ---- helpers ----
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
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rGv, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), p_nstep1 = erl_p(csL, csG, 1), usable = sum(okc), ob = ob, f0 = f0,
       sL = sL[, okc, drop = FALSE], sG = sG[, okc, drop = FALSE], n = nn[okc], pairs = pr[, okc, drop = FALSE])
}
env_df <- function(curve, r, obs, S) data.frame(
  curve = curve, r_km = r, obs = obs, sim_mean = rowMeans(S),
  lo = apply(S, 1, quantile, 0.025), hi = apply(S, 1, quantile, 0.975),
  z_raw = (obs - rowMeans(S)) / apply(S, 1, sd))

REPORT <- file.path(OUT_DIR, "42_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("42_empirical_m1 | ", format(Sys.time()), "\n", sep = "")
    cat(sprintf("frozen settings: E3 (12 df), O1, dataset A, simplified window, nominal coordinates | nsim %d | ERL nstep 2 | alpha %.2f | seed %d\n",
                NSIM, TESTS$alpha, SEED_E))
    cat(sprintf("L_inhom r %.1f-%.1f km, g_inhom r %.1f-%.1f km, step %.1f | pcf kernel sd %.3f km\n",
                min(rL), max(rL), min(rGv), max(rGv), R_STEP_KM, bw))
    
    cat("\n== 0. Data checks ==\n")
    cat(sprintf("sites %d (expect %d) | window %s, area %.1f km2 | NA-covariate sites %d (must be 0)\n",
                npoints(X_obs_km), N_A, Window(X_obs_km)$type, area.owin(W_km), sum(na_cov_pts(X_obs_km))))
    stopifnot(npoints(X_obs_km) == N_A, sum(na_cov_pts(X_obs_km)) == 0, Window(X_obs_km)$type == "polygonal")
    cat("observed pairs within 1, 2, 3, 5 km:", paste(vapply(PAIR_R, function(r) pairs_at(X_obs_km, r), numeric(1)), collapse = " "),
        "(project_state: 7 and 46 at 1 and 2 km)\n")
    
    cat("\n== 1. M1 environmental Poisson (E3): coefficients (MPL, Berman-Turner) ==\n")
    fitP <- tryCatch(ppm(make_quad(X_obs_km), trend = formula_E3, covariates = covs), error = function(e) NULL)
    if (is.null(fitP)) cat("ppm fit FAILED\n") else {
      b <- coef(fitP); V <- tryCatch(vcov(fitP), error = function(e) NULL)
      se <- if (is.null(V)) rep(NA_real_, length(b)) else sqrt(diag(V))
      tab <- data.frame(term = names(b), coef = unname(b), se = unname(se), z = unname(b / se))
      print(tab, digits = 3, row.names = FALSE)
      write.csv(tab, file.path(OUT_DIR, "42_m1_coef.csv"), row.names = FALSE)
      f0 <- fast_fit(X_obs_km)
      cat(sprintf("max |coef| difference fast_fit vs ppm: %.2e | integral of fitted intensity %.2f (n = %d)\n",
                  max(abs(f0$beta - b[names(f0$beta)])), integral(lam_im_of(f0$beta)), N_A))
    }
    
    cat("\n== 2. M1 generative test (nsim ", NSIM, "; main session) ==\n", sep = "")
    t0 <- unname(proc.time()["elapsed"])
    R1 <- emp_m1(X_obs_km, NSIM, bw, rGv, SEED_E)
    cat(sprintf("done in %.0f s | usable simulations %d of %d\n", unname(proc.time()["elapsed"]) - t0, R1$usable, NSIM))
    cat(sprintf("ERL combined global envelope p (nstep 2): %.4f | nstep 1 (sensitivity): %.4f | deficiency at 0.05: %s\n",
                R1$p, R1$p_nstep1, if (is.na(R1$p)) "NA" else if (R1$p < TESTS$alpha) "YES" else "NO"))
    cat(sprintf("simulated n: mean %.1f sd %.1f (observed %d)\n", mean(R1$n), sd(R1$n), N_A))
    
    zL <- zsum(R1$ob$L, R1$sL); zG <- zsum(R1$ob$g, R1$sG)
    zfl <- max(abs(c(zL$flo, zG$flo)), na.rm = TRUE)
    rr <- c(rL, rGv); zr <- abs(c(zL$raw, zG$raw)); zr[!is.finite(zr)] <- NA
    ia <- which.max(zr)
    cat(sprintf("max |z| (floored sd) %.2f | raw max |z| %.2f at r = %.1f km on the %s curve\n",
                zfl, zr[ia], rr[ia], if (ia <= length(rL)) "L_inhom" else "g_inhom"))
    
    CV <- rbind(env_df("L_inhom_minus_r", rL, R1$ob$L, R1$sL), env_df("g_inhom", rGv, R1$ob$g, R1$sG))
    write.csv(CV, file.path(OUT_DIR, "42_m1_curves.csv"), row.names = FALSE)
    cat("\n-- Observed vs simulated (pointwise 95% band, display only; inference uses the global ERL test) --\n")
    show <- CV[round(CV$r_km, 6) %in% c(1, 2, 3, 4, 5, 7, 10), ]
    print(show, digits = 3, row.names = FALSE)
    
    cat("\n-- Pair counts: observed vs M1 simulations --\n")
    po <- vapply(PAIR_R, function(r) pairs_at(X_obs_km, r), numeric(1))
    PC <- data.frame(r_km = PAIR_R, observed = po, sim_mean = rowMeans(R1$pairs),
                     sim_lo = apply(R1$pairs, 1, quantile, 0.025), sim_hi = apply(R1$pairs, 1, quantile, 0.975),
                     share_sim_ge_obs = vapply(seq_along(PAIR_R), function(i) mean(R1$pairs[i, ] >= po[i]), numeric(1)))
    print(PC, digits = 3, row.names = FALSE)
    
    write.csv(data.frame(test = "M1", n_obs = N_A, nsim = NSIM, usable = R1$usable, p_erl_nstep2 = R1$p,
                         p_erl_nstep1 = R1$p_nstep1, zmax_floored = zfl, zraw_max = zr[ia], zraw_r_km = rr[ia],
                         zraw_curve = if (ia <= length(rL)) "L" else "g", seed = SEED_E,
                         date = format(Sys.time())), f_res, row.names = FALSE)
    
    cat("\n== Warning tally (main session) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "42_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")