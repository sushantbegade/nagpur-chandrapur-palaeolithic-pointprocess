# =====================================================================
# 38_recovery_cells_log_and_smoke.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing;
#           writes R/recovery_cells.R and data_derived/recovery_obs_cache.rds)
# Purpose : (A) cache the O2/O3 offset layers as im objects (workers load no terra);
#           (B) append clarification rows to docs/protocol_deviation_log.csv and
#               update the NA-covariate row;
#           (C) write R/recovery_cells.R: truth generators for all 19 cells,
#               M1/M3/M4 tests with per-cell g-range, M4 profile search over
#               GEYER_R_KM x GEYER_SAT, hierarchy classification, rec_task;
#           (D) smoke test: 19 cells x N_SMOKE reps, nsim 19, 4 workers.
# Outputs : outputs/38_report.txt, 38_smoke.csv, 38_sessionInfo.txt,
#           protocol_deviation_log_before_38.csv; data_derived/recovery_obs_cache.rds;
#           R/recovery_cells.R; docs/protocol_deviation_log.csv (appended/amended)
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           The smoke test is NOT a size or power estimate (nsim 19, 2 reps).
#           logLik() of Geyer ppm fits as a profile criterion is UNVERIFIED against
#           documentation; the smoke test prints the selected (r, sat) per cell.
#           Restart R before running. terra is used with :: in the main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))
suppressPackageStartupMessages(library(parallel))

OUT_DIR    <- dir_of("10_simulation_recovery", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
TESTS_FILE <- file.path(REPO_DIR, "R", "recovery_tests.R"); stopifnot(file.exists(TESTS_FILE))
CELLS_FILE <- file.path(REPO_DIR, "R", "recovery_cells.R")
OBS_FILE   <- file.path(REPO_DIR, "data_derived", "recovery_obs_cache.rds")
N_WORK     <- 4
N_SMOKE    <- 2
NSIM_SMOKE <- 19

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULTS <- setNames(gm$mult, gm$cell); stopifnot(all(c("G1","G2","G3","G4","C1","C3") %in% names(MULTS)))
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")
env_dir <- file.path(REPO_DIR, "data_derived", "env", "500m")
to_im2 <- function(r) {
  m <- terra::as.matrix(r, wide = TRUE); m <- m[nrow(m):1, , drop = FALSE]
  im(m, xcol = terra::xFromCol(r, seq_len(ncol(r))) / 1000, yrow = rev(terra::yFromRow(r, seq_len(nrow(r)))) / 1000)
}
OBS <- list(logE_O2 = to_im2(terra::rast(file.path(env_dir, "logE_O2.tif"))),
            logE_O3 = to_im2(terra::rast(file.path(env_dir, "logE_O3.tif"))))

# ---- functions written to R/recovery_cells.R (sourced by workers and by script 39) ----
REC_NUMF <- c("sec", "n_raw", "n", "n_dropped", "p1", "p3", "p4", "zmax3", "zmax4", "usable1", "usable3", "usable4",
              "sigma2", "alpha", "gamma", "r_sel", "sat_sel", "sec1", "sec3", "sec4")
env_eta <- function(beta) Reduce(`+`, Map(function(a, b) a * b, covs, beta))
b0_for  <- function(eta, n) log(n / integral(exp(eta)[W_km, drop = FALSE]))
jitter_pat <- function(P) {
  probs <- c(FIELD = 64, AUTHOR_REPORTED = 73, GIS_MAP = 60) / 197
  n <- npoints(P); if (n == 0) return(P)
  cls <- sample(names(probs), n, replace = TRUE, prob = probs)
  s <- unname(JITTER_SIGMA_M[cls]) / 1000
  x <- P$x; y <- P$y
  for (k in seq_len(n)) {
    for (tr in seq_len(50)) {
      xn <- x[k] + rnorm(1, 0, s[k]); yn <- y[k] + rnorm(1, 0, s[k])
      if (inside.owin(xn, yn, W_km)) { x[k] <- xn; y[k] <- yn; break }
    }
  }
  ppp(x, y, window = W_km, check = FALSE)
}
gen_cell <- function(cc, seed) {
  row <- RECOVERY_CELLS[RECOVERY_CELLS$cell == cc, ]
  n <- row$n
  if (row$mechanism == "Geyer") return(to_poly(sim_geyer(MULTS[[cc]], row$gamma, row$r_km, row$sat, n, 1, seed)))
  if (cc %in% c("P1", "P2")) {
    tr <- make_truth(unname(RECOVERY_ENV$levels[row$env]), n)
    set.seed(seed); return(to_poly(rpoispp(tr$lamW)))
  }
  if (cc == "N1") {
    eta <- env_eta(truth$beta) + OBS$logE_O3; lam <- exp(eta + b0_for(eta, n))[W_km, drop = FALSE]
    set.seed(seed); return(to_poly(rpoispp(lam)))
  }
  if (cc == "N2") {
    tr1 <- make_truth(1, 100)
    set.seed(SEEDS$master + 5002); beta2 <- setNames(rnorm(length(covs), 0, RECOVERY_ENV$coef_sd), names(covs))
    eta2 <- env_eta(beta2); lam2 <- exp(eta2 + b0_for(eta2, 97))[W_km, drop = FALSE]
    set.seed(seed)
    A <- to_poly(rpoispp(tr1$lamW)); B <- to_poly(rpoispp(lam2))
    return(ppp(c(A$x, B$x), c(A$y, B$y), window = W_km, check = FALSE))
  }
  if (cc == "N3") { set.seed(seed); return(jitter_pat(to_poly(rpoispp(truth$lamW)))) }
  if (cc == "N4") {
    set.seed(SEEDS$master + 5001); b <- rnorm(2, 0, RECOVERY_ENV$coef_sd)
    eta <- env_eta(truth$beta) + b[1] * covs$z_lslope * covs$z_lslope + b[2] * covs$z_lslope * covs$lith_basalt
    lam <- exp(eta + b0_for(eta, n))[W_km, drop = FALSE]
    set.seed(seed); return(to_poly(rpoispp(lam)))
  }
  if (cc == "N5") {
    eta <- env_eta(truth$beta) + OBS$logE_O2; mu <- fill_im(eta + b0_for(eta, n))
    set.seed(seed)
    X <- to_poly(rLGCP(model = "exp", mu = mu - row$sigma2 / 2, param = list(var = row$sigma2, scale = row$scale_km),
                       win = W_km, eps = EPS_LGCP, nsim = 1, drop = TRUE))
    return(jitter_pat(X))
  }
  if (row$mechanism == "LGCP") {
    tr <- make_truth(unname(RECOVERY_ENV$levels[row$env]), n)
    mu <- fill_im(env_eta(tr$beta) + tr$b0)
    set.seed(seed)
    return(to_poly(rLGCP(model = "exp", mu = mu - row$sigma2 / 2, param = list(var = row$sigma2, scale = row$scale_km),
                         win = W_km, eps = EPS_LGCP, nsim = 1, drop = TRUE)))
  }
  stop("no generator for cell ", cc)
}
rG_for <- function(n) {
  half <- PCF_STOYAN_C / sqrt(n / area_km2)
  lower <- max(R_MIN_KM[["A"]], ceiling(half * 2) / 2)
  RGRID[RGRID >= lower & RGRID <= R_MAX_DECISION_KM]
}
fast_curves_r <- function(Xp, ftb, bw, rGv) {
  Lh <- do.call(Linhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "isotropic"), L_ARGS))
  gh <- do.call(pcfinhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "Ripley", bw = bw), G_ARGS_NS))
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(rGv, 6), round(gh$r, 6))])
}
curves_test_r <- function(Xobs, sims, bw, rGv) {
  f0 <- fast_fit(Xobs); ob <- fast_curves_r(Xobs, f0, bw, rGv)
  sL <- matrix(NA_real_, length(rL), length(sims)); sG <- matrix(NA_real_, length(rGv), length(sims))
  dsim <- integer(length(sims))
  for (k in seq_along(sims)) {
    Xr <- to_poly(sims[[k]]); Xs <- clean_pat(Xr); dsim[k] <- npoints(Xr) - npoints(Xs)
    fs <- fast_fit(Xs); cs <- fast_curves_r(Xs, fs, bw, rGv); sL[, k] <- cs$L; sG[, k] <- cs$g
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rGv, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  zmax <- NA_real_
  if (sum(okc) >= 2) {
    zL <- (ob$L - rowMeans(sL[, okc, drop = FALSE])) / apply(sL[, okc, drop = FALSE], 1, sd)
    zG <- (ob$g - rowMeans(sG[, okc, drop = FALSE])) / apply(sG[, okc, drop = FALSE], 1, sd)
    zz <- abs(c(zL, zG)); zz <- zz[is.finite(zz)]
    if (length(zz)) zmax <- max(zz)
  }
  list(p = erl_p(csL, csG, 2), usable = sum(okc), zmax = zmax, drop_sim = sum(dsim))
}
m1_test_r <- function(Xobs, nsim, bw, rGv, seed) {
  set.seed(seed)
  f0 <- fast_fit(Xobs); ob <- fast_curves_r(Xobs, f0, bw, rGv)
  lam_im <- lam_im_of(f0$beta)
  sL <- matrix(NA_real_, length(rL), nsim); sG <- matrix(NA_real_, length(rGv), nsim)
  for (k in seq_len(nsim)) {
    Xs <- to_poly(rpoispp(lam_im)); fs <- fast_fit(Xs); cs <- fast_curves_r(Xs, fs, bw, rGv)
    sL[, k] <- cs$L; sG[, k] <- cs$g
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rGv, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc))
}
m3_test_r <- function(Xobs, nsim, bw, rGv, seed) {
  set.seed(seed)
  fit <- tryCatch(kppm(Xobs, trend = formula_E3, clusters = "LGCP", model = "exponential",
                       covariates = covs, method = "clik2"), error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA_real_, usable = 0, zmax = NA_real_, sigma2 = NA_real_, alpha = NA_real_))
  s2 <- unname(fit$par["sigma2"]); al <- unname(fit$par["alpha"])
  eta <- trend_im(coef(fit)); n <- npoints(Xobs)
  S <- tryCatch(rLGCP(model = "exp", mu = eta - s2 / 2, param = list(var = s2, scale = al), win = W_km,
                      eps = EPS_LGCP, saveLambda = TRUE, nsim = nsim, drop = TRUE), error = function(e) NULL)
  if (is.null(S)) return(list(p = NA_real_, usable = 0, zmax = NA_real_, sigma2 = s2, alpha = al))
  sims <- lapply(as_list(S), function(P) rpoint(n, f = attr(P, "Lambda"), win = W_km, forcewin = TRUE))
  r <- curves_test_r(Xobs, sims, bw, rGv)
  c(r, list(sigma2 = s2, alpha = al))
}
m4_search_test <- function(Xobs, nsim, bw, rGv, seed) {
  set.seed(seed)
  grid <- expand.grid(r = GEYER_R_KM, sat = GEYER_SAT)
  fits <- lapply(seq_len(nrow(grid)), function(k)
    tryCatch(geyer_fit(Xobs, r = grid$r[k], sat = grid$sat[k], rbord = 0), error = function(e) NULL))
  ll <- vapply(fits, function(f) if (is.null(f)) NA_real_ else
    tryCatch(as.numeric(suppressWarnings(logLik(f))), error = function(e) NA_real_), numeric(1))
  if (all(is.na(ll))) return(list(p = NA_real_, usable = 0, zmax = NA_real_, gamma = NA_real_, r_sel = NA_real_, sat_sel = NA_real_))
  kb <- which.max(ll); fit <- fits[[kb]]
  gh <- exp(unname(tail(coef(fit), 1))); n <- npoints(Xobs)
  S <- simulate(fit, nsim = nsim, drop = TRUE, verbose = FALSE, progress = FALSE,
                control = rmhcontrol(nrep = NREP_RMH, p = 1, expand = 1), start = rmhstart(n.start = n))
  r <- curves_test_r(Xobs, as_list(S), bw, rGv)
  c(r, list(gamma = gh, r_sel = grid$r[kb], sat_sel = grid$sat[kb]))
}
run_hierarchy <- function(X, nsim, bw, rGv, seed) {
  R <- list(p1 = NA_real_, p3 = NA_real_, p4 = NA_real_, zmax3 = NA_real_, zmax4 = NA_real_,
            usable1 = NA_real_, usable3 = NA_real_, usable4 = NA_real_, sigma2 = NA_real_, alpha = NA_real_,
            gamma = NA_real_, r_sel = NA_real_, sat_sel = NA_real_, sec1 = NA_real_, sec3 = NA_real_, sec4 = NA_real_,
            class = "failed")
  t1 <- unname(proc.time()["elapsed"])
  a <- m1_test_r(X, nsim, bw, rGv, seed + 1)
  R$sec1 <- unname(proc.time()["elapsed"]) - t1; R$p1 <- a$p; R$usable1 <- a$usable
  if (is.na(a$p)) return(R)
  if (a$p >= TESTS$alpha) { R$class <- "Poisson"; return(R) }
  t3 <- unname(proc.time()["elapsed"])
  b <- m3_test_r(X, nsim, bw, rGv, seed + 2)
  R$sec3 <- unname(proc.time()["elapsed"]) - t3; R$p3 <- b$p; R$usable3 <- b$usable; R$zmax3 <- b$zmax
  R$sigma2 <- b$sigma2; R$alpha <- b$alpha
  if (is.na(b$p)) return(R)
  if (b$p >= TESTS$alpha) { R$class <- "LGCP"; return(R) }
  t4 <- unname(proc.time()["elapsed"])
  d <- m4_search_test(X, nsim, bw, rGv, seed + 3)
  R$sec4 <- unname(proc.time()["elapsed"]) - t4; R$p4 <- d$p; R$usable4 <- d$usable; R$zmax4 <- d$zmax
  R$gamma <- d$gamma; R$r_sel <- d$r_sel; R$sat_sel <- d$sat_sel
  if (is.na(d$p) || is.na(d$zmax) || is.na(b$zmax)) return(R)
  R$class <- if (d$p >= TESTS$alpha && d$zmax <= 0.5 * b$zmax) "Geyer" else "inadequate"
  R
}
rec_task <- function(job) {
  t0 <- unname(proc.time()["elapsed"])
  cc <- job$cell; i <- job$i; ci <- match(cc, RECOVERY_CELLS$cell)
  sb <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 100000
  n_t <- RECOVERY_CELLS$n[ci]; bw <- bw_for_n(n_t); rGv <- rG_for(n_t)
  out <- tryCatch({
    X0 <- gen_cell(cc, sb + 1000 * ci + i); X <- clean_pat(X0)
    res <- run_hierarchy(X, job$nsim, bw, rGv, sb + 500000 + 1000 * ci + i)
    c(list(n_raw = npoints(X0), n = npoints(X), n_dropped = npoints(X0) - npoints(X), err = NA_character_), res)
  }, error = function(e) list(class = "failed", err = conditionMessage(e)))
  c(list(cell = cc, i = i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
rec_rows_df <- function(L) do.call(rbind, lapply(L, function(x) {
  d <- data.frame(cell = as.character(x$cell), i = as.numeric(x$i), stringsAsFactors = FALSE)
  for (f in REC_NUMF) d[[f]] <- if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  d$class <- as.character(x$class); d$err <- as.character(x$err)
  d
}))
init_worker2 <- function(repo, cache, tests, cells, obs_cache, mults) {
  assign("REPO_DIR", repo, envir = globalenv())
  Sys.setenv(RECOVERY_CACHE = cache)
  source(file.path(repo, "R", "recovery_worker.R"), local = globalenv())
  source(file.path(repo, "R", "recovery_functions.R"), local = globalenv())
  source(tests, local = globalenv())
  source(cells, local = globalenv())
  assign("OBS", readRDS(obs_cache), envir = globalenv())
  assign("MULTS", mults, envir = globalenv())
  TRUE
}
dump(c("REC_NUMF", "env_eta", "b0_for", "jitter_pat", "gen_cell", "rG_for", "fast_curves_r", "curves_test_r", "m1_test_r",
       "m3_test_r", "m4_search_test", "run_hierarchy", "rec_task", "rec_rows_df", "init_worker2"), file = CELLS_FILE)
saveRDS(OBS, OBS_FILE)

WARN <- new.env(); WARN$msgs <- character(0)
note_w <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }

REPORT <- file.path(OUT_DIR, "38_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("38_recovery_cells_log_and_smoke | ", format(Sys.time()), "\n", sep = "")
    
    ## ---------------- A. offset cache ----------------
    cat("\n== A. Offset layers (500 m) vs covariate grid ==\n")
    for (nm in names(OBS)) {
      im1 <- OBS[[nm]]; c1 <- covs[[1]]
      same <- identical(dim(im1$v), dim(c1$v)) && isTRUE(all.equal(im1$xcol, c1$xcol)) && isTRUE(all.equal(im1$yrow, c1$yrow))
      nadiff <- sum(is.na(im1$v) != is.na(c1$v))
      cat(sprintf("  %s: dim %s, same grid as covariates: %s | NA-mask mismatches vs z_elev_p1: %d | range %.4f-%.4f\n",
                  nm, paste(dim(im1$v), collapse = "x"), same, nadiff, min(im1$v, na.rm = TRUE), max(im1$v, na.rm = TRUE)))
      if (!same) stop("offset grid differs from covariate grid; stop and paste output")
    }
    cat("saved", OBS_FILE, "| wrote", CELLS_FILE, "\n")
    cat("cell table rows:", nrow(RECOVERY_CELLS), "| Geyer multipliers available:", paste(names(MULTS), collapse = " "), "\n")
    cat("g-range lower bounds (km) at n = 197, 114, 75:", paste(vapply(c(197, 114, 75), function(n) min(rG_for(n)), 1), collapse = " "), "\n")
    
    ## ---------------- B. deviation log ----------------
    cat("\n== B. protocol_deviation_log.csv ==\n")
    D <- "2026-11-10"; AB <- "author-delegated (assistant decision; clarification, no frozen item changed)"
    ROWS <- data.frame(
      date = D,
      item = c("Clarification: floors activated (nsim 99, 100 reps per cell) with staged extension to 200 reps",
               "Clarification: operational definitions of recovery classification, power, false-positive rate, mechanism recovery and 'substantially reduces'",
               "Clarification: nuisance-cell construction and seed scheme for the recovery run"),
      reason = c("Script 36 (4 workers): M3 median 206 s and M4 157 s per test at nsim 99; arithmetic estimate for frozen counts (nsim 199, 200 reps) about 8-10 days. The freeze allows floors (nsim 99, 100 reps) if a timing pilot requires. The null pilot validated size at nsim 99. Stage 2 (reps 101-200, nsim 99, new seeds) only for cells whose primary Gate 6 metric has a Wilson 95% interval containing its threshold; rule fixed before any recovery output.",
                 "Per replicate: Poisson if M1 not rejected (ERL p >= 0.05); LGCP if M1 rejected and M3 not; Geyer if M1 and M3 rejected, M4 at the profile-selected (r, sat) (maximum pseudolikelihood over r in {2,4,7,10} km x sat in {2,4}, rbord 0) has ERL p >= 0.05 and max pointwise |z| over the joint L and g curves <= 50% of the M3 value; otherwise inadequate; failed fit = failed. Power = P(Geyer | Geyer truth); false positive = P(Geyer | non-Gibbs cell); mechanism recovery = P(correct class | truth). g lower bound per cell follows the frozen rule max(r_min, Stoyan half-width rounded up to 0.5 km): 2.0 km (n 197), 2.5 km (n 114), 3.0 km (n 75).",
                 "N1/N5 truths include the frozen logE_O3/logE_O2 layers (500 m), intercept recalibrated to expected n 197, and are fitted WITHOUT offset (primary pipeline is O1). N2 second pattern uses an independent N(0, 0.3^2) coefficient draw (seed MASTER + 5002; 100 + 97 expected points). N3/N5 jitter: Gaussian per axis, source-class proportions 64/73/60 of 197, sigma 10/25/250 m, redraw if outside window (50 tries). N4 extra terms seed MASTER + 5001. Seeds: truth MASTER + 10000 + 100000 + 1000 x cell index + rep; tests +500000; empirical stage offsets >= 30000."),
      effect_on_inference = c("nsim 99 gives ERL p resolution 0.01. Sequential extension can slightly inflate error rates; stages are reported separately and pooled. No archaeological model has been fitted.",
                              "Classification makes Geyer power conditional on M3 rejection: if LGCP absorbs Gibbs attraction, power is low and the paper pivots to residual spatial dependence (Gate 6 rule). No archaeological model has been fitted.",
                              "N-cell severity (N2 second coefficient set, jitter) is an assistant design choice; results are interpreted as robustness checks of the primary pipeline. No archaeological model has been fitted."),
      approved_by = AB, stringsAsFactors = FALSE)
    if (!identical(names(LOG0), LOG_COLS)) {
      cat("Existing log columns differ; NOT writing. Rows saved to 38_log_rows_for_manual_append.csv\n")
      write.csv(ROWS, file.path(OUT_DIR, "38_log_rows_for_manual_append.csv"), row.names = FALSE)
    } else {
      file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_38.csv"), overwrite = TRUE)
      LOG1 <- LOG0
      idx <- grep("^Clarification: NA-covariate points removed", LOG1$item)
      if (length(idx) == 1) {
        LOG1$reason[idx] <- paste("Script 34 M3 rep 4 failed: fast_fit dropped a data point (214 of 215). Script 35 found 0 NA-covariate pixel centres inside the window.",
                                  "Script 36: 2 of 1e6 uniform points in the window have NA covariates, at 0.001-0.025 km from the boundary (boundary pixels whose centre lies outside the polygon),",
                                  "so the mechanism is supported (n = 2 events). Cleaning dropped 1 point (M3) and 2 points (M4) over 40 x 99 simulations each.")
        cat("Amended NA-covariate row", idx, "\n")
      } else cat("NA-covariate row not found uniquely (match count", length(idx), "); not amended\n")
      NEW <- ROWS[!(ROWS$item %in% LOG1$item), , drop = FALSE]
      write.csv(rbind(LOG1, NEW), f_log, row.names = FALSE)
      cat(sprintf("Appended %d row(s); log now %d rows (backup in outputs).\n", nrow(NEW), nrow(LOG1) + nrow(NEW)))
    }
    
    ## ---------------- D. smoke test ----------------
    cat("\n== D. Smoke test: ", nrow(RECOVERY_CELLS), " cells x ", N_SMOKE, " reps, nsim ", NSIM_SMOKE, ", ", N_WORK, " workers ==\n", sep = "")
    cl <- makeCluster(N_WORK)
    ini <- clusterCall(cl, init_worker2, REPO_DIR, CACHE_FILE, TESTS_FILE, CELLS_FILE, OBS_FILE, MULTS)
    cat("workers initialised:", all(unlist(ini)), "\n")
    ex <- clusterCall(cl, function() c(exists("rec_task"), exists("gen_cell"), exists("OBS"), exists("MULTS"), "terra" %in% loadedNamespaces()))
    cat("worker check (rec_task, gen_cell, OBS, MULTS present; terra loaded):", paste(ex[[1]], collapse = " "), "\n")
    jobs <- unlist(lapply(RECOVERY_CELLS$cell, function(cc) lapply(900 + seq_len(N_SMOKE), function(i)
      list(cell = cc, i = i, nsim = NSIM_SMOKE))), recursive = FALSE)
    t0 <- unname(proc.time()["elapsed"])
    res <- parLapplyLB(cl, jobs, rec_task)
    wall <- unname(proc.time()["elapsed"]) - t0
    SM <- rec_rows_df(res)
    write.csv(SM, file.path(OUT_DIR, "38_smoke.csv"), row.names = FALSE)
    cat(sprintf("smoke wall %.0f s for %d replicates\n\n", wall, nrow(SM)))
    print(SM[, c("cell", "i", "class", "n_raw", "n", "n_dropped", "p1", "p3", "p4", "sigma2", "alpha", "gamma", "r_sel", "sat_sel", "sec", "err")],
          digits = 3, row.names = FALSE)
    cat("\n-- Per-cell means (n realised vs target; stage times in s at nsim", NSIM_SMOKE, ") --\n")
    AG <- do.call(rbind, lapply(split(SM, SM$cell), function(d) data.frame(
      cell = d$cell[1], target_n = RECOVERY_CELLS$n[RECOVERY_CELLS$cell == d$cell[1]], mean_n = mean(d$n),
      failed = sum(d$class == "failed"), sec_mean = mean(d$sec), sec1 = mean(d$sec1, na.rm = TRUE),
      sec3 = mean(d$sec3, na.rm = TRUE), sec4 = mean(d$sec4, na.rm = TRUE))))
    print(AG[match(RECOVERY_CELLS$cell, AG$cell), ], digits = 3, row.names = FALSE)
    cat("\nreplicates failed:", sum(SM$class == "failed"), "of", nrow(SM), "\n")
    if (any(SM$class == "failed")) { cat("errors:\n"); print(unique(SM[SM$class == "failed", c("cell", "err")]), row.names = FALSE) }
    cat("\n== Warning tally (main session only; worker warnings not captured) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink(); if (exists("cl")) try(stopCluster(cl), silent = TRUE) })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "38_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")