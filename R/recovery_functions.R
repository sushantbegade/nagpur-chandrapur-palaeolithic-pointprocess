
# recovery_functions.R -- written by 21_window_speed_test.R; sourced by recovery scripts and workers.
covsW <- lapply(covs, function(im) im[W_km, drop = FALSE])
area_km2 <- area.owin(W_km)
bw_for_n <- function(n) (PCF_STOYAN_C / sqrt(n / area_km2)) / sqrt(5)   # kernel sd; the Stoyan value is the half-width
L_ARGS    <- list(renormalise = TRUE, normpower = 1)
G_ARGS_NS <- list(kernel = "epanechnikov", zerocor = "convolution", divisor = "r", renormalise = TRUE, normpower = 2)
DUMMY_FIXED <- default.dummy(X_obs_km, eps = QUAD_KM)
truth <- make_truth()
make_quad <- function(Xp) quadscheme(Xp, dummy = DUMMY_FIXED, method = "grid")
fast_fit <- function(Xp) {
  Q <- make_quad(Xp)
  w <- w.quad(Q); isd <- is.data(Q); xq <- x.quad(Q); yq <- y.quad(Q)
  Z <- vapply(covs, function(im) lookup.im(im, xq, yq, naok = TRUE), numeric(length(xq)))
  ok <- stats::complete.cases(Z) & w > 0
  X <- cbind(`(Intercept)` = 1, Z)
  gfit <- suppressWarnings(stats::glm.fit(X[ok, , drop = FALSE], (isd / w)[ok], weights = w[ok], family = poisson(), intercept = FALSE))
  beta <- gfit$coefficients
  Xd <- X[isd & ok, , drop = FALSE]
  if (nrow(Xd) != npoints(Xp)) stop("data points dropped in fast_fit: ", nrow(Xd), " of ", npoints(Xp))
  list(beta = beta, lam_data = as.numeric(exp(Xd %*% beta)))
}
geyer_fit <- function(Xp, r, sat, rbord = 10)
  ppm(make_quad(Xp), trend = formula_E3, interaction = Geyer(r = r, sat = sat), covariates = covs, rbord = rbord)
lam_im_of <- function(beta) exp(beta[1] + Reduce(`+`, Map(function(a, b) a * b, covsW, beta[-1])))
fast_curves <- function(Xp, ftb, bw) {
  Lh <- do.call(Linhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "isotropic"), L_ARGS))
  gh <- do.call(pcfinhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "Ripley", bw = bw), G_ARGS_NS))
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(rG, 6), round(gh$r, 6))])
}
erl_p <- function(csL, csG, nstep) tryCatch(
  as.numeric(attr(GET::global_envelope_test(curve_sets = list(csL, csG), type = "erl", alpha = TESTS$alpha, nstep = nstep), "p")),
  error = function(e) NA_real_)
m1_test <- function(Xobs, nsim, bw, seed) {
  if (nsim < 19) stop("GET ERL test at alpha 0.05 needs nsim >= 19")
  set.seed(seed)
  f0 <- fast_fit(Xobs); ob <- fast_curves(Xobs, f0, bw)
  lam_im <- lam_im_of(f0$beta)
  sL <- matrix(NA_real_, length(rL), nsim); sG <- matrix(NA_real_, length(rG), nsim)
  for (k in seq_len(nsim)) {
    Xs <- to_poly(rpoispp(lam_im)); fs <- fast_fit(Xs); cs <- fast_curves(Xs, fs, bw)
    sL[, k] <- cs$L; sG[, k] <- cs$g
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), p_nstep1 = erl_p(csL, csG, 1), usable = sum(okc))
}
null_rep <- function(i, nsim) {
  t0 <- unname(proc.time()["elapsed"])
  set.seed(SEEDS$master + SEEDS$offsets["empirical"] + 1000 + i)
  X <- to_poly(rpoispp(truth$lamW))
  r <- m1_test(X, nsim, bw_for_n(RECOVERY_ENV$expected_n), SEEDS$master + SEEDS$offsets["empirical"] + 2000 + i)
  c(p = r$p, p_nstep1 = r$p_nstep1, usable = r$usable, n = npoints(X), sec = unname(proc.time()["elapsed"]) - t0)
}
time_probe <- function(npat) {
  ts <- tf <- tc <- numeric(npat)
  bw <- bw_for_n(RECOVERY_ENV$expected_n)
  for (i in seq_len(npat)) {
    set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 3000 + i)
    ts[i] <- system.time(X <- to_poly(rpoispp(truth$lamW)))["elapsed"]
    tf[i] <- system.time(f <- fast_fit(X))["elapsed"]
    tc[i] <- system.time(cv <- fast_curves(X, f, bw))["elapsed"]
  }
  c(sim = median(ts), fit = median(tf), curves = median(tc), vertices = sum(vapply(W_km$bdry, function(b) length(b$x), 1L)), area = area_km2)
}

