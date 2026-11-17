LGCP_VAR <-
1.5
LGCP_SCALE <-
8
EPS_LGCP <-
0.5
NREP_RMH <-
500000
as_list <-
function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
fill_im <-
function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
trend_im <-
function(b) fill_im(Reduce(`+`, Map(function(a, nmc) a * unname(b[nmc]), covs, names(covs))) + unname(b["(Intercept)"]))
lgcp_mu <-
function() fill_im(Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta)) + truth$b0)
na_cov_pts <-
function(P) {
  if (npoints(P) == 0) return(logical(0))
  Z <- do.call(cbind, lapply(covs, function(im) lookup.im(im, P$x, P$y, naok = TRUE)))
  rowSums(is.na(Z)) > 0
}
clean_pat <-
function(P) { bad <- na_cov_pts(P); if (any(bad)) P[!bad] else P }
sim_geyer <-
function(mult, gamma, r, sat, n_start, nsim, seed, nrep = NREP_RMH) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}
curves_test <-
function(Xobs, sims, bw) {
  Xo <- clean_pat(Xobs)
  f0 <- fast_fit(Xo); ob <- fast_curves(Xo, f0, bw)
  sL <- matrix(NA_real_, length(rL), length(sims)); sG <- matrix(NA_real_, length(rG), length(sims))
  dsim <- integer(length(sims))
  for (k in seq_along(sims)) {
    Xr <- to_poly(sims[[k]]); Xs <- clean_pat(Xr); dsim[k] <- npoints(Xr) - npoints(Xs)
    fs <- fast_fit(Xs); cs <- fast_curves(Xs, fs, bw)
    sL[, k] <- cs$L; sG[, k] <- cs$g
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc), drop_sim = sum(dsim), sims_dropped = sum(dsim > 0))
}
m3_test <-
function(Xobs, nsim, bw, seed) {
  if (nsim < 19) stop("ERL test at alpha 0.05 needs nsim >= 19")
  set.seed(seed)
  Xc <- clean_pat(Xobs)
  fit <- tryCatch(kppm(Xc, trend = formula_E3, clusters = "LGCP", model = "exponential",
                       covariates = covs, method = "clik2"), error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA_real_, usable = 0L, drop_sim = NA_integer_, sims_dropped = NA_integer_,
                                sigma2 = NA_real_, alpha = NA_real_))
  s2 <- unname(fit$par["sigma2"]); al <- unname(fit$par["alpha"])
  eta <- trend_im(coef(fit)); n <- npoints(Xc)
  S <- rLGCP(model = "exp", mu = eta - s2 / 2, param = list(var = s2, scale = al), win = W_km,
             eps = EPS_LGCP, saveLambda = TRUE, nsim = nsim, drop = TRUE)
  sims <- lapply(as_list(S), function(P) rpoint(n, f = attr(P, "Lambda"), win = W_km, forcewin = TRUE))
  r <- curves_test(Xc, sims, bw)
  c(r, list(sigma2 = s2, alpha = al))
}
m4_test <-
function(Xobs, r_km, sat, nsim, bw, seed) {
  if (nsim < 19) stop("ERL test at alpha 0.05 needs nsim >= 19")
  set.seed(seed)
  Xc <- clean_pat(Xobs)
  fit <- tryCatch(geyer_fit(Xc, r = r_km, sat = sat, rbord = 0), error = function(e) NULL)
  if (is.null(fit)) return(list(p = NA_real_, usable = 0L, drop_sim = NA_integer_, sims_dropped = NA_integer_, gamma = NA_real_))
  gh <- exp(unname(tail(coef(fit), 1))); n <- npoints(Xc)
  S <- simulate(fit, nsim = nsim, drop = TRUE, verbose = FALSE, progress = FALSE,
                control = rmhcontrol(nrep = NREP_RMH, p = 1, expand = 1), start = rmhstart(n.start = n))
  r <- curves_test(Xc, as_list(S), bw)
  c(r, list(gamma = gh))
}
null_task <-
function(job) {
  t0 <- unname(proc.time()["elapsed"])
  base <- SEEDS$master + unname(SEEDS$offsets["empirical"])
  bw <- bw_for_n(RECOVERY_ENV$expected_n)
  out <- tryCatch({
    if (job$test == "M3") {
      set.seed(base + 3000 + job$i)
      X <- clean_pat(to_poly(rLGCP(model = "exp", mu = lgcp_mu() - LGCP_VAR / 2,
                                   param = list(var = LGCP_VAR, scale = LGCP_SCALE), win = W_km,
                                   eps = EPS_LGCP, nsim = 1, drop = TRUE)))
      r <- m3_test(X, job$nsim, bw, base + 4000 + job$i)
      list(n = npoints(X), p = r$p, usable = r$usable, e1 = r$sigma2, e2 = r$alpha, dsim = r$drop_sim, err = NA_character_)
    } else {
      X <- to_poly(sim_geyer(job$mult, 2.5, 3, 3, 197, 1, base + 5000 + job$i))
      r <- m4_test(X, 3, 3, job$nsim, bw, base + 6000 + job$i)
      list(n = npoints(X), p = r$p, usable = r$usable, e1 = r$gamma, e2 = NA_real_, dsim = r$drop_sim, err = NA_character_)
    }
  }, error = function(e) list(n = NA_real_, p = NA_real_, usable = 0L, e1 = NA_real_, e2 = NA_real_,
                              dsim = NA_integer_, err = conditionMessage(e)))
  c(list(test = job$test, i = job$i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
