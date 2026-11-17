CONT <-
c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", 
"z_ldist_riv5")
LITH <-
c("lith_basalt", "lith_sandstone", "lith_carbonate", "lith_other"
)
as_list <-
function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
fill_im <-
function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
make_setting_m3 <-
function(nm) {
  if (nm == "R01_redep") return(mk(nm, X_obs_km[!SITES_G$sens_redeposited], W_km, covs, DUMMY_FIXED))
  if (nm == "R06_eroded5") { Wer <- erosion(W_km, 5); Xe <- X_obs_km[Wer]; return(mk(nm, Xe, Wer, covs, default.dummy(Xe, eps = QUAD_KM))) }
  if (nm == "R08_E2") return(mk(nm, X_obs_km, W_km, covs[c(CONT, LITH)], DUMMY_FIXED))
  if (nm == "R09_E4") return(mk(nm, X_obs_km, W_km, covs[c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", LITH)], DUMMY_FIXED))
  stop("unknown setting ", nm)
}
m3_run <-
function(S, nsim, seed, chunk = 100) {
  set.seed(seed)
  fm <- as.formula(paste("~", paste(names(S$cov), collapse = " + ")))
  fit <- kppm(S$X, trend = fm, clusters = "LGCP", model = "exponential", covariates = S$cov, method = "clik2")
  s2 <- unname(fit$par["sigma2"]); al <- unname(fit$par["alpha"]); b <- coef(fit)
  eta <- Reduce(`+`, Map(function(a, nmc) a * unname(b[nmc]), S$cov, names(S$cov))) + unname(b["(Intercept)"])
  mu <- fill_im(eta); n <- npoints(S$X)
  ft <- rob_fit(S$X, S); ob <- rob_curves(S$X, ft, S)
  sL <- matrix(NA_real_, length(S$rL), nsim); sG <- matrix(NA_real_, length(S$rG), nsim); nn <- rep(NA_real_, nsim)
  k0 <- 0
  while (k0 < nsim) {
    m <- min(chunk, nsim - k0)
    Sl <- as_list(rLGCP(model = "exp", mu = mu - s2 / 2, param = list(var = s2, scale = al), win = S$W, eps = 0.5,
                        saveLambda = TRUE, nsim = m, drop = TRUE))
    for (j in seq_len(m)) {
      res <- tryCatch({
        Xs <- rob_poly(rpoint(n, f = attr(Sl[[j]], "Lambda"), win = S$W, forcewin = TRUE), S)
        Xc <- rob_clean(Xs, S); fs <- rob_fit(Xc, S); list(cs = rob_curves(Xc, fs, S), n = npoints(Xc))
      }, error = function(e) NULL)
      if (!is.null(res)) { sL[, k0 + j] <- res$cs$L; sG[, k0 + j] <- res$cs$g; nn[k0 + j] <- res$n }
    }
    k0 <- k0 + m
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = S$rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = S$rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc), sigma2 = s2, alpha = al, n_mean = mean(nn[okc]), n = n)
}
m3_task <-
function(job) {
  t0 <- unname(proc.time()["elapsed"])
  S <- make_setting_m3(job$setting)
  r <- tryCatch(m3_run(S, job$nsim, job$seed), error = function(e) list(err = conditionMessage(e)))
  c(list(setting = job$setting, sec = unname(proc.time()["elapsed"]) - t0), r)
}
m3_row <-
function(x) {
  g <- function(f) if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  data.frame(setting = x$setting, n = g("n"), nsim = g("usable"), p = g("p"), sigma2 = g("sigma2"), alpha = g("alpha"),
             n_sim_mean = g("n_mean"), sec = g("sec"), err = if (is.null(x$err)) NA_character_ else as.character(x$err), stringsAsFactors = FALSE)
}
