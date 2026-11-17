GMIN_FLOOR <-
2
mk <-
function(name, X, W, cv, dum, Lmin = 1, corrL = "isotropic", corrG = "Ripley") {
  n <- npoints(X); a <- area.owin(W); half <- PCF_STOYAN_C / sqrt(n / a); gmin <- max(Lmin, GMIN_FLOOR, ceiling(half * 2) / 2)
  list(name = name, X = X, W = W, cov = cv, n = n, area = a, bw = half / sqrt(5),
       rL = RGRID[RGRID >= Lmin & RGRID <= R_MAX_DECISION_KM], rG = RGRID[RGRID >= gmin & RGRID <= R_MAX_DECISION_KM],
       corrL = corrL, corrG = corrG, dum = dum, covW = lapply(cv, function(im) im[W, drop = FALSE]))
}
cov_na <-
function(P, S) {
  if (npoints(P) == 0) return(logical(0))
  Z <- do.call(cbind, lapply(S$cov, function(im) lookup.im(im, P$x, P$y, naok = TRUE)))
  rowSums(is.na(Z)) > 0
}
rob_clean <-
function(P, S) { b <- cov_na(P, S); if (any(b)) P[!b] else P }
rob_poly <-
function(P, S) { Q <- P[S$W]; if (Window(Q)$type != "polygonal") stop("not polygonal"); Q }
rob_fit <-
function(Xp, S) {
  Q <- quadscheme(Xp, dummy = S$dum, method = "grid")
  w <- w.quad(Q); isd <- is.data(Q); xq <- x.quad(Q); yq <- y.quad(Q)
  Z <- vapply(S$cov, function(im) lookup.im(im, xq, yq, naok = TRUE), numeric(length(xq)))
  ok <- stats::complete.cases(Z) & w > 0
  X1 <- cbind(`(Intercept)` = 1, Z)
  g <- suppressWarnings(stats::glm.fit(X1[ok, , drop = FALSE], (isd / w)[ok], weights = w[ok], family = poisson(), intercept = FALSE))
  Xd <- X1[isd & ok, , drop = FALSE]
  if (nrow(Xd) != npoints(Xp)) stop("data points dropped in rob_fit: ", nrow(Xd), " of ", npoints(Xp))
  list(beta = g$coefficients, lam_data = as.numeric(exp(Xd %*% g$coefficients)))
}
rob_curves <-
function(Xp, ft, S) {
  Lh <- do.call(Linhom, c(list(Xp, lambda = ft$lam_data, r = RGRID, correction = S$corrL), L_ARGS))
  gh <- do.call(pcfinhom, c(list(Xp, lambda = ft$lam_data, r = RGRID, correction = S$corrG, bw = S$bw), G_ARGS_NS))
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(S$rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(S$rG, 6), round(gh$r, 6))])
}
make_setting <-
function(nm) {
  if (nm == "R02_datasetC") return(mk(nm, X_obs_km[SITES_G$dataset_C], W_km, covs, DUMMY_FIXED, Lmin = 2))
  if (nm == "R12_quad05") return(mk(nm, X_obs_km, W_km, covs, default.dummy(X_obs_km, eps = 0.5)))
  stop("unknown setting ", nm)
}
chunk_run <-
function(S, nsim, seed) {
  set.seed(seed)
  ft <- rob_fit(S$X, S)
  lam <- exp(Reduce(`+`, Map(function(a, b) a * b, S$covW, ft$beta[-1])) + ft$beta[1])
  sL <- matrix(NA_real_, length(S$rL), nsim); sG <- matrix(NA_real_, length(S$rG), nsim)
  for (k in seq_len(nsim)) {
    res <- tryCatch({ Xs <- rob_poly(rpoispp(lam), S); Xc <- rob_clean(Xs, S); fs <- rob_fit(Xc, S); rob_curves(Xc, fs, S) }, error = function(e) NULL)
    if (!is.null(res)) { sL[, k] <- res$L; sG[, k] <- res$g }
  }
  list(sL = sL, sG = sG)
}
chunk_task <-
function(job) {
  t0 <- unname(proc.time()["elapsed"])
  S <- make_setting(job$setting); r <- chunk_run(S, job$nsim, job$seed)
  c(list(setting = job$setting, ch = job$ch, sec = unname(proc.time()["elapsed"]) - t0), r)
}
