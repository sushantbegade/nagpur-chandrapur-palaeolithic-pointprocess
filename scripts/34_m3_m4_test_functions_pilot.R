# =====================================================================
# 34_m3_m4_test_functions_pilot.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : define and time the M3 (LGCP) and M4 (Geyer) generative tests used
#           in recovery and empirical confirmatory analysis:
#             M3: kppm clik2 fit -> rLGCP (eps = 0.5, saveLambda) -> rpoint(n)
#             M4: Geyer MPL fit (rbord 0) -> simulate.ppm (p = 1, expand = 1, n.start = n)
#           statistic = frozen ERL joint test (L_inhom + g_inhom, nstep = 2),
#           inhomogeneous curves from the trend-only Poisson refit per pattern.
#           Pilot: N_PILOT synthetic patterns per mechanism, nsim = NSIM_PILOT.
# Outputs : outputs/34_report.txt, 34_pilot.csv, 34_sessionInfo.txt
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           This is a timing / sanity pilot, NOT a size estimate.
#           Conditioning on n (M3, M4) differs from the unconditional M1 test:
#           a documented design asymmetry. rLGCP eps = 0.5 per decision after
#           script 33. Restart R before running. Main session only (no workers).
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))

OUT_DIR <- dir_of("10_simulation_recovery", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEED7      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 23000
N_PILOT    <- 4
NSIM_PILOT <- 99
NREP_RMH   <- 5e5
LGCP_VAR   <- 1.5
LGCP_SCALE <- 8
EPS_LGCP   <- 0.5

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULT_G3 <- gm$mult[gm$cell == "G3"]; stopifnot(length(MULT_G3) == 1)

# ---- helpers ----
WARN <- new.env(); WARN$msgs <- character(0)
note_w  <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
fill_im <- function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
trend_im <- function(b) fill_im(Reduce(`+`, Map(function(a, nmc) a * unname(b[nmc]), covs, names(covs))) + unname(b["(Intercept)"]))

# generic ERL test of an observed pattern against a list of simulated patterns
curves_test <- function(Xobs, sims, bw) {
  f0 <- fast_fit(Xobs); ob <- fast_curves(Xobs, f0, bw)
  sL <- matrix(NA_real_, length(rL), length(sims)); sG <- matrix(NA_real_, length(rG), length(sims))
  for (k in seq_along(sims)) {
    Xs <- to_poly(sims[[k]]); fs <- fast_fit(Xs); cs <- fast_curves(Xs, fs, bw)
    sL[, k] <- cs$L; sG[, k] <- cs$g
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc))
}

m3_test <- function(Xobs, nsim, bw, seed) {
  if (nsim < 19) stop("ERL test at alpha 0.05 needs nsim >= 19")
  set.seed(seed)
  fit <- tryCatch(kppm(Xobs, trend = formula_E3, clusters = "LGCP", model = "exponential",
                       covariates = covs, method = "clik2"), error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA_real_, usable = 0L, sigma2 = NA_real_, alpha = NA_real_))
  s2 <- unname(fit$par["sigma2"]); al <- unname(fit$par["alpha"])
  eta <- trend_im(coef(fit)); n <- npoints(Xobs)
  S <- rLGCP(model = "exp", mu = eta - s2 / 2, param = list(var = s2, scale = al), win = W_km,
             eps = EPS_LGCP, saveLambda = TRUE, nsim = nsim, drop = TRUE)
  sims <- lapply(as_list(S), function(P) rpoint(n, f = attr(P, "Lambda"), win = W_km, forcewin = TRUE))
  r <- curves_test(Xobs, sims, bw)
  c(r, list(sigma2 = s2, alpha = al))
}

m4_test <- function(Xobs, r_km, sat, nsim, bw, seed) {
  if (nsim < 19) stop("ERL test at alpha 0.05 needs nsim >= 19")
  set.seed(seed)
  fit <- tryCatch(geyer_fit(Xobs, r = r_km, sat = sat, rbord = 0), error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA_real_, usable = 0L, gamma = NA_real_))
  gh <- exp(unname(tail(coef(fit), 1))); n <- npoints(Xobs)
  S <- simulate(fit, nsim = nsim, drop = TRUE, verbose = FALSE, progress = FALSE,
                control = rmhcontrol(nrep = NREP_RMH, p = 1, expand = 1), start = rmhstart(n.start = n))
  r <- curves_test(Xobs, as_list(S), bw)
  c(r, list(gamma = gh))
}

sim_geyer <- function(mult, gamma, r, sat, n_start, nsim, seed, nrep = NREP_RMH) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}

REPORT <- file.path(OUT_DIR, "34_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("34_m3_m4_test_functions_pilot | ", format(Sys.time()), "\n", sep = "")
    cat(sprintf("N_PILOT %d | nsim %d | rmh nrep %.0e | rLGCP eps %.2f\n\n", N_PILOT, NSIM_PILOT, NREP_RMH, EPS_LGCP))
    bw <- bw_for_n(RECOVERY_ENV$expected_n)
    mu0 <- Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta)) + truth$b0
    mu_fill <- fill_im(mu0)
    PIL <- data.frame()
    
    ## M3 pilot: L4 truth
    cat("== M3 test on L4 truth (sigma2 1.5, scale 8 km), 0.5 km generator grid ==\n")
    set.seed(SEED7 + 1)
    TL <- rLGCP(model = "exp", mu = mu_fill - LGCP_VAR / 2, param = list(var = LGCP_VAR, scale = LGCP_SCALE),
                win = W_km, eps = EPS_LGCP, nsim = N_PILOT, drop = TRUE)
    TLl <- lapply(as_list(TL), to_poly)
    for (k in seq_along(TLl)) {
      Xk <- TLl[[k]]; t0 <- unname(proc.time()["elapsed"])
      r <- tryCatch(m3_test(Xk, NSIM_PILOT, bw, SEED7 + 100 + k), error = function(e) { cat("  M3 FAILED:", conditionMessage(e), "\n"); NULL })
      sec <- unname(proc.time()["elapsed"]) - t0
      if (is.null(r)) next
      cat(sprintf("  rep %d: n %d | sigma2_hat %.3f alpha_hat %.2f | ERL p %.3f (usable sims %d) | %.0f s\n",
                  k, npoints(Xk), r$sigma2, r$alpha, r$p, r$usable, sec))
      PIL <- rbind(PIL, data.frame(test = "M3", truth = "L4", rep = k, n = npoints(Xk), est1 = r$sigma2, est2 = r$alpha,
                                   p = r$p, usable = r$usable, sec = sec))
    }
    
    ## M4 pilot: G3 truth, fitted at the true (r, sat)
    cat("\n== M4 test on G3 truth, fitted at r = 3, sat = 3 (rbord 0) ==\n")
    TG <- as_list(sim_geyer(MULT_G3, 2.5, 3, 3, 197, N_PILOT, SEED7 + 2))
    for (k in seq_along(TG)) {
      Xk <- to_poly(TG[[k]]); t0 <- unname(proc.time()["elapsed"])
      r <- tryCatch(m4_test(Xk, 3, 3, NSIM_PILOT, bw, SEED7 + 200 + k), error = function(e) { cat("  M4 FAILED:", conditionMessage(e), "\n"); NULL })
      sec <- unname(proc.time()["elapsed"]) - t0
      if (is.null(r)) next
      cat(sprintf("  rep %d: n %d | gamma_hat %.2f | ERL p %.3f (usable sims %d) | %.0f s\n",
                  k, npoints(Xk), r$gamma, r$p, r$usable, sec))
      PIL <- rbind(PIL, data.frame(test = "M4", truth = "G3", rep = k, n = npoints(Xk), est1 = r$gamma, est2 = NA_real_,
                                   p = r$p, usable = r$usable, sec = sec))
    }
    if (nrow(PIL)) {
      write.csv(PIL, file.path(OUT_DIR, "34_pilot.csv"), row.names = FALSE)
      cat("\n-- Timing summary (seconds per test at nsim =", NSIM_PILOT, ") --\n")
      print(aggregate(sec ~ test, data = PIL, FUN = function(x) c(mean = mean(x), max = max(x))), digits = 3)
    }
    cat("\n== Warning tally (whole run) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "34_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")