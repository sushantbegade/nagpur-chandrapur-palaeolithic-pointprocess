jitter_obs2 <-
function(seed, sig_class_m) {
  set.seed(seed)
  sig <- unname(sig_class_m[SITES_G$coordinate_source]) / 1000
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
jit_task2 <-
function(job) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch({
    jj <- jitter_obs2(job$seed_j, job$sig); r <- jit_test(jj$X, job$nsim, job$seed_t)
    list(n = npoints(jj$X), p = r$p, usable = r$usable, n_sim = r$n_sim, po2 = r$pobs[1], po3 = r$pobs[2], po5 = r$pobs[3],
         ps2 = r$psim[1], ps3 = r$psim[2], ps5 = r$psim[3], fallbacks = jj$fallbacks, med_disp_km = median(jj$disp), err = NA_character_)
  }, error = function(e) list(err = conditionMessage(e)))
  c(list(scn = job$scn, i = job$i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
jrow2 <-
function(L) do.call(rbind, lapply(L, function(x) {
  g <- function(f) if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  data.frame(scn = x$scn, i = x$i, sec = x$sec, n = g("n"), p = g("p"), usable = g("usable"), n_sim = g("n_sim"),
             po2 = g("po2"), po3 = g("po3"), po5 = g("po5"), ps2 = g("ps2"), ps3 = g("ps3"), ps5 = g("ps5"),
             fallbacks = g("fallbacks"), med_disp_km = g("med_disp_km"), err = if (is.null(x$err)) NA_character_ else as.character(x$err), stringsAsFactors = FALSE)
}))
