pairs_at <-
function(P, r) sum(closepaircounts(P, r)) / 2
jitter_obs <-
function(seed) {
  set.seed(seed)
  sig <- unname(JITTER_SIGMA_M[SITES_G$coordinate_source]) / 1000
  x <- X_obs_km$x; y <- X_obs_km$y; fb <- 0L; disp <- numeric(length(x))
  for (i in seq_along(x)) {
    ok <- FALSE
    for (tr in seq_len(50)) {
      xn <- x[i] + rnorm(1, 0, sig[i]); yn <- y[i] + rnorm(1, 0, sig[i])
      if (inside.owin(xn, yn, W_km)) {
        z <- vapply(covs, function(im) lookup.im(im, xn, yn, naok = TRUE), numeric(1))
        if (!anyNA(z)) { disp[i] <- sqrt((xn - x[i])^2 + (yn - y[i])^2); x[i] <- xn; y[i] <- yn; ok <- TRUE; break }
      }
    }
    if (!ok) fb <- fb + 1L
  }
  list(X = ppp(x, y, window = W_km, check = FALSE), fallbacks = fb, disp = disp)
}
jit_test <-
function(X, nsim, seed) {
  S <- mk("J", X, W_km, covs, DUMMY_FIXED)
  ft <- rob_fit(X, S); ob <- rob_curves(X, ft, S)
  lam <- exp(Reduce(`+`, Map(function(a, b) a * b, S$covW, ft$beta[-1])) + ft$beta[1])
  set.seed(seed)
  sL <- matrix(NA_real_, length(S$rL), nsim); sG <- matrix(NA_real_, length(S$rG), nsim)
  pr <- matrix(NA_real_, 3, nsim); nn <- rep(NA_real_, nsim)
  for (k in seq_len(nsim)) {
    res <- tryCatch({
      Xs <- rob_poly(rpoispp(lam), S); Xc <- rob_clean(Xs, S); fs <- rob_fit(Xc, S)
      list(cs = rob_curves(Xc, fs, S), n = npoints(Xc), pr = vapply(c(2, 3, 5), function(r) pairs_at(Xc, r), numeric(1)))
    }, error = function(e) NULL)
    if (!is.null(res)) { sL[, k] <- res$cs$L; sG[, k] <- res$cs$g; nn[k] <- res$n; pr[, k] <- res$pr }
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = S$rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = S$rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc), n_sim = mean(nn[okc]),
       pobs = vapply(c(2, 3, 5), function(r) pairs_at(X, r), numeric(1)), psim = rowMeans(pr[, okc, drop = FALSE]))
}
jit_task <-
function(job) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch({
    if (job$type == "jitter") { jj <- jitter_obs(job$seed_j); X <- jj$X; fb <- jj$fallbacks; md <- median(jj$disp) } else { X <- X_obs_km; fb <- 0L; md <- 0 }
    r <- jit_test(X, job$nsim, job$seed_t)
    list(n = npoints(X), p = r$p, usable = r$usable, n_sim = r$n_sim, po2 = r$pobs[1], po3 = r$pobs[2], po5 = r$pobs[3],
         ps2 = r$psim[1], ps3 = r$psim[2], ps5 = r$psim[3], fallbacks = fb, med_disp_km = md, err = NA_character_)
  }, error = function(e) list(err = conditionMessage(e)))
  c(list(type = job$type, i = job$i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
jrow_df <-
function(L) do.call(rbind, lapply(L, function(x) {
  g <- function(f) if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  data.frame(type = x$type, i = x$i, sec = x$sec, n = g("n"), p = g("p"), usable = g("usable"), n_sim = g("n_sim"),
             po2 = g("po2"), po3 = g("po3"), po5 = g("po5"), ps2 = g("ps2"), ps3 = g("ps3"), ps5 = g("ps5"),
             fallbacks = g("fallbacks"), med_disp_km = g("med_disp_km"),
             err = if (is.null(x$err)) NA_character_ else as.character(x$err), stringsAsFactors = FALSE)
}))
