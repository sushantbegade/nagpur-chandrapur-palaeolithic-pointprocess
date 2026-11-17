REC_NUMF <-
c("sec", "n_raw", "n", "n_dropped", "p1", "p3", "p4", "zmax3", 
"zmax4", "zraw3", "zraw4", "zarg_r3", "zarg_r4", "zarg_cv3", 
"zarg_cv4", "usable1", "usable3", "usable4", "sigma2", "alpha", 
"gamma", "r_sel", "sat_sel", "sec1", "sec3", "sec4")
num1 <-
function(x) if (is.null(x)) NA_real_ else as.numeric(x)
zsum <-
function(obs, S, floorf = 0.25) {
  m <- rowMeans(S); s <- apply(S, 1, sd)
  med <- median(s[is.finite(s) & s > 0]); if (!is.finite(med)) med <- 1
  sf <- pmax(s, floorf * med)
  list(raw = (obs - m) / s, flo = (obs - m) / sf)
}
env_eta <-
function(beta) Reduce(`+`, Map(function(a, b) a * b, covs, beta))
b0_for <-
function(eta, n) log(n / integral(exp(eta)[W_km, drop = FALSE]))
jitter_pat <-
function(P) {
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
gen_cell <-
function(cc, seed) {
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
rG_for <-
function(n) {
  half <- PCF_STOYAN_C / sqrt(n / area_km2)
  lower <- max(R_MIN_KM[["A"]], ceiling(half * 2) / 2)
  RGRID[RGRID >= lower & RGRID <= R_MAX_DECISION_KM]
}
fast_curves_r <-
function(Xp, ftb, bw, rGv) {
  Lh <- do.call(Linhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "isotropic"), L_ARGS))
  gh <- do.call(pcfinhom, c(list(Xp, lambda = ftb$lam_data, r = RGRID, correction = "Ripley", bw = bw), G_ARGS_NS))
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(rGv, 6), round(gh$r, 6))])
}
curves_test_r <-
function(Xobs, sims, bw, rGv) {
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
  zmax <- NA_real_; zraw <- NA_real_; zr <- NA_real_; zcv <- NA_real_
  if (sum(okc) >= 2) {
    zL <- zsum(ob$L, sL[, okc, drop = FALSE]); zG <- zsum(ob$g, sG[, okc, drop = FALSE])
    fl <- abs(c(zL$flo, zG$flo)); rw <- abs(c(zL$raw, zG$raw)); rr <- c(rL, rGv)
    fl[!is.finite(fl)] <- NA; rw[!is.finite(rw)] <- NA
    if (any(!is.na(fl))) zmax <- max(fl, na.rm = TRUE)
    if (any(!is.na(rw))) { ia <- which.max(rw); zraw <- rw[ia]; zr <- rr[ia]; zcv <- if (ia <= length(rL)) 1 else 2 }
  }
  list(p = erl_p(csL, csG, 2), usable = sum(okc), zmax = zmax, zraw = zraw, zarg_r = zr, zarg_cv = zcv, drop_sim = sum(dsim))
}
m1_test_r <-
function(Xobs, nsim, bw, rGv, seed) {
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
m3_test_r <-
function(Xobs, nsim, bw, rGv, seed) {
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
m4_search_test <-
function(Xobs, nsim, bw, rGv, seed) {
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
run_hierarchy <-
function(X, nsim, bw, rGv, seed, force = FALSE) {
  R <- list(p1 = NA_real_, p3 = NA_real_, p4 = NA_real_, zmax3 = NA_real_, zmax4 = NA_real_,
            zraw3 = NA_real_, zraw4 = NA_real_, zarg_r3 = NA_real_, zarg_r4 = NA_real_, zarg_cv3 = NA_real_, zarg_cv4 = NA_real_,
            usable1 = NA_real_, usable3 = NA_real_, usable4 = NA_real_, sigma2 = NA_real_, alpha = NA_real_,
            gamma = NA_real_, r_sel = NA_real_, sat_sel = NA_real_, sec1 = NA_real_, sec3 = NA_real_, sec4 = NA_real_,
            class = "failed_M1")
  t1 <- unname(proc.time()["elapsed"])
  a <- tryCatch(m1_test_r(X, nsim, bw, rGv, seed + 1), error = function(e) list(p = NA_real_, usable = 0))
  R$sec1 <- unname(proc.time()["elapsed"]) - t1; R$p1 <- num1(a$p); R$usable1 <- num1(a$usable)
  if (is.na(R$p1)) return(R)
  if (!force && R$p1 >= TESTS$alpha) { R$class <- "Poisson"; return(R) }
  R$class <- "failed_M3"
  t3 <- unname(proc.time()["elapsed"])
  b <- tryCatch(m3_test_r(X, nsim, bw, rGv, seed + 2), error = function(e) list(p = NA_real_, usable = 0))
  R$sec3 <- unname(proc.time()["elapsed"]) - t3; R$p3 <- num1(b$p); R$usable3 <- num1(b$usable)
  R$zmax3 <- num1(b$zmax); R$zraw3 <- num1(b$zraw); R$zarg_r3 <- num1(b$zarg_r); R$zarg_cv3 <- num1(b$zarg_cv)
  R$sigma2 <- num1(b$sigma2); R$alpha <- num1(b$alpha)
  if (is.na(R$p3)) return(R)
  if (!force && R$p3 >= TESTS$alpha) { R$class <- "LGCP"; return(R) }
  R$class <- "failed_M4"
  t4 <- unname(proc.time()["elapsed"])
  d <- tryCatch(m4_search_test(X, nsim, bw, rGv, seed + 3), error = function(e) list(p = NA_real_, usable = 0))
  R$sec4 <- unname(proc.time()["elapsed"]) - t4; R$p4 <- num1(d$p); R$usable4 <- num1(d$usable)
  R$zmax4 <- num1(d$zmax); R$zraw4 <- num1(d$zraw); R$zarg_r4 <- num1(d$zarg_r); R$zarg_cv4 <- num1(d$zarg_cv)
  R$gamma <- num1(d$gamma); R$r_sel <- num1(d$r_sel); R$sat_sel <- num1(d$sat_sel)
  if (is.na(R$p4) || is.na(R$zmax4) || is.na(R$zmax3)) return(R)
  R$class <- if (R$p1 >= TESTS$alpha) "Poisson" else if (R$p3 >= TESTS$alpha) "LGCP" else
    if (R$p4 >= TESTS$alpha && R$zmax4 <= 0.5 * R$zmax3) "Geyer" else "inadequate"
  R
}
rec_task <-
function(job) {
  t0 <- unname(proc.time()["elapsed"])
  cc <- job$cell; i <- job$i; ci <- match(cc, RECOVERY_CELLS$cell)
  sb <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 100000
  n_t <- RECOVERY_CELLS$n[ci]; bw <- bw_for_n(n_t); rGv <- rG_for(n_t)
  out <- tryCatch({
    X0 <- gen_cell(cc, sb + 1000 * ci + i); X <- clean_pat(X0)
    res <- run_hierarchy(X, job$nsim, bw, rGv, sb + 500000 + 1000 * ci + i, force = isTRUE(job$force))
    c(list(n_raw = npoints(X0), n = npoints(X), n_dropped = npoints(X0) - npoints(X), err = NA_character_), res)
  }, error = function(e) list(class = "failed_gen", err = conditionMessage(e)))
  c(list(cell = cc, i = i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
rec_rows_df <-
function(L) do.call(rbind, lapply(L, function(x) {
  d <- data.frame(cell = as.character(x$cell), i = as.numeric(x$i), stringsAsFactors = FALSE)
  for (f in REC_NUMF) d[[f]] <- if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  d$class <- as.character(x$class); d$err <- as.character(x$err)
  d
}))
init_worker2 <-
function(repo, cache, tests, cells, obs_cache, mults) {
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
