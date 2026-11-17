# =====================================================================
# 29_lgcp_bias_ncond_expand_probe.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : (A) LGCP sigma2 over-estimation (script 28): generator field
#               variance vs 'var'; kppm clik2 vs mincon on identical patterns;
#               clik2 at 4x intensity (finite-sample test);
#           (B) conditional-on-n LGCP simulation: simulate.kppm(n.cond,
#               covariates = covs) and direct rLGCP(n.cond) from fitted trend;
#           (C) Geyer conditional simulation: default window expansion vs
#               rmhcontrol(expand = 1): n constancy, outside-W points, pair ratio.
# Outputs : outputs/29_report.txt, 29_lgcp_bias_probe.csv,
#           29_geyer_expand_probe.csv, 29_sessionInfo.txt
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           mincon is a DIAGNOSTIC only; the frozen LGCP estimator is clik2.
#           rLGCP saveLambda attribute, n.cond semantics and expand = 1
#           behaviour are UNVERIFIED against documentation; this script tests them.
#           Restart R before running. Main session only (no workers).
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

SEED0      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 6000   # same 8 G3 patterns as scripts 26-28
SEED2      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 9000
NREP_PROD  <- 5e5
N_LGCP     <- 20
N_FIT      <- 10
N_FIT_HI   <- 5
HI_MULT    <- 4
N_PAT_GY   <- 4
N_SIM      <- 10
N26        <- c(209, 174, 197, 180, 183, 190, 199, 174)   # VERIFIED n of the 8 patterns

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULT_G3 <- gm$mult[gm$cell == "G3"]; stopifnot(length(MULT_G3) == 1)

# ---- helpers ----
try_ <- function(label, expr) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      cat(sprintf("  [%s] warning: %s\n", label, conditionMessage(w))); invokeRestart("muffleWarning") }),
    error = function(e) { cat(sprintf("[%s] FAILED: %s\n", label, conditionMessage(e))); NULL })
  cat(sprintf("[%s] %.1f s\n", label, unname(proc.time()["elapsed"]) - t0))
  out
}
quiet   <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NULL)
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
qs      <- function(x) quantile(x, c(0.25, 0.75), na.rm = TRUE)
fill_im <- function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
sim_geyer <- function(mult, gamma, r, sat, n_start, nsim, seed, nrep = NREP_PROD) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}
gen_lgcp <- function(v, s, shift, nsim, seed, save) {
  set.seed(seed)
  rLGCP(model = "exp", mu = mu_fill + shift - v / 2, param = list(var = v, scale = s),
        win = W_km, nsim = nsim, saveLambda = save, drop = TRUE)
}
fit_lgcp <- function(X, method)
  quiet(kppm(X, trend = formula_E3, clusters = "LGCP", model = "exponential", covariates = covs, method = method))
pars_of  <- function(f) if (is.null(f)) c(NA_real_, NA_real_) else unname(c(f$par["sigma2"], f$par["alpha"]))
summ <- function(label, s2, al, v, s)
  cat(sprintf("  %-24s ok %d/%d | sigma2 median %.3f (IQR %.3f-%.3f) truth %.2f | alpha median %.2f (IQR %.2f-%.2f) truth %.0f\n",
              label, sum(!is.na(s2)), length(s2), median(s2, na.rm = TRUE), qs(s2)[1], qs(s2)[2], v,
              median(al, na.rm = TRUE), qs(al)[1], qs(al)[2], s))
rows_of <- function(cell, setting, nn, s2, al, v, s)
  data.frame(cell = cell, setting = setting, rep = seq_along(s2), n = nn[seq_along(s2)],
             sigma2_hat = s2, alpha_hat = al, var_true = v, scale_true = s)
npairs3  <- function(P) sum(closepaircounts(P, 3)) / 2
n_raw    <- function(S) vapply(as_list(S), npoints, numeric(1))
n_poly   <- function(S) vapply(as_list(S), function(P) npoints(to_poly(P)), numeric(1))
pairs_of <- function(S) vapply(as_list(S), function(P) npairs3(to_poly(P)), numeric(1))
out_of   <- function(S) sum(vapply(as_list(S), function(P) npoints(P) - sum(inside.owin(P$x, P$y, W_km)), numeric(1)))
sim_cond <- function(fit, X, nsim, seed, ctrl_name) {
  set.seed(seed)
  tryCatch(suppressWarnings(simulate(fit, nsim = nsim, drop = TRUE, verbose = FALSE, progress = FALSE,
                                     control = switch(ctrl_name,
                                                      default = rmhcontrol(nrep = NREP_PROD, p = 1),
                                                      expand1 = rmhcontrol(nrep = NREP_PROD, p = 1, expand = 1)),
                                     start = rmhstart(n.start = npoints(X)))),
           error = function(e) { cat("  simulate FAILED:", conditionMessage(e), "\n"); NULL })
}

REPORT <- file.path(OUT_DIR, "29_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("29_lgcp_bias_ncond_expand_probe | ", format(Sys.time()), "\n", sep = "")
  cat("G3 multiplier", MULT_G3, "| production nrep", NREP_PROD, "\n\n")
  
  ## ---------------- A. LGCP sigma2 bias ----------------
  cat("== A. LGCP sigma2 over-estimation ==\n")
  mu0 <- Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta)) + truth$b0
  mu_fill <- fill_im(mu0)
  LG <- list(L1 = c(var = 0.5, scale = 3), L4 = c(var = 1.5, scale = 8))
  LGR <- data.frame(); KP <- NULL
  for (nm in names(LG)) {
    v <- unname(LG[[nm]]["var"]); s <- unname(LG[[nm]]["scale"])
    cat(sprintf("\n-- %s: var %.1f, scale %.0f km --\n", nm, v, s))
    PL1 <- try_(sprintf("rLGCP %s x1, %d sims, saveLambda", nm, N_LGCP), gen_lgcp(v, s, 0, N_LGCP, SEED2 + 10, TRUE))
    if (is.null(PL1)) next
    vv <- vapply(as_list(PL1), function(P) tryCatch(
      var(as.vector((log(attr(P, "Lambda")) - (mu_fill - v / 2))[W_km, drop = TRUE])),
      error = function(e) NA_real_), numeric(1))
    cat(sprintf("  generator field variance inside W (pixel sample var): mean %.3f, truth %.2f, usable %d of %d\n",
                mean(vv, na.rm = TRUE), v, sum(!is.na(vv)), length(vv)))
    PLl <- lapply(as_list(PL1), to_poly); nn <- vapply(PLl, npoints, numeric(1))
    cat(sprintf("  x1 mean n %.1f sd %.1f\n", mean(nn), sd(nn)))
    nf <- min(N_FIT, length(PLl))
    s2c <- rep(NA_real_, nf); alc <- s2c; s2m <- s2c; alm <- s2c
    for (k in seq_len(nf)) {
      fc <- fit_lgcp(PLl[[k]], "clik2");  pc <- pars_of(fc); s2c[k] <- pc[1]; alc[k] <- pc[2]
      fm <- fit_lgcp(PLl[[k]], "mincon"); pm <- pars_of(fm); s2m[k] <- pm[1]; alm[k] <- pm[2]
      if (nm == "L4" && is.null(KP) && !is.null(fc)) KP <- list(fit = fc, X = PLl[[k]])
    }
    summ("clik2, x1", s2c, alc, v, s)
    summ("mincon (diagnostic), x1", s2m, alm, v, s)
    LGR <- rbind(LGR, rows_of(nm, "clik2_x1", nn, s2c, alc, v, s), rows_of(nm, "mincon_x1", nn, s2m, alm, v, s))
    PH <- try_(sprintf("rLGCP %s x%d intensity, %d sims", nm, HI_MULT, N_FIT_HI),
               gen_lgcp(v, s, log(HI_MULT), N_FIT_HI, SEED2 + 20, FALSE))
    if (!is.null(PH)) {
      PHl <- lapply(as_list(PH), to_poly); nh <- vapply(PHl, npoints, numeric(1))
      cat(sprintf("  x%d mean n %.1f\n", HI_MULT, mean(nh)))
      s2h <- rep(NA_real_, length(PHl)); alh <- s2h
      for (k in seq_along(PHl)) {
        th <- try_(sprintf("clik2 x%d fit %d", HI_MULT, k), fit_lgcp(PHl[[k]], "clik2"))
        ph <- pars_of(th); s2h[k] <- ph[1]; alh[k] <- ph[2]
      }
      summ(sprintf("clik2, x%d", HI_MULT), s2h, alh, v, s)
      LGR <- rbind(LGR, rows_of(nm, sprintf("clik2_x%d", HI_MULT), nh, s2h, alh, v, s))
    }
  }
  if (nrow(LGR)) write.csv(LGR, file.path(OUT_DIR, "29_lgcp_bias_probe.csv"), row.names = FALSE)
  
  ## ---------------- B. LGCP n.cond ----------------
  cat("\n== B. Conditional-on-n LGCP simulation (first L4 clik2 fit) ==\n")
  if (!is.null(KP)) {
    n_obs <- npoints(KP$X); p_obs <- npairs3(KP$X)
    cat(sprintf("n_obs %d, pairs <= 3 km observed %d\n", n_obs, p_obs))
    Sa <- try_("simulate.kppm n.cond + covariates", simulate(KP$fit, nsim = 20, n.cond = n_obs, covariates = covs,
                                                             drop = TRUE, verbose = FALSE))
    if (!is.null(Sa)) { na <- n_poly(Sa)
    cat(sprintf("  (a) n range %d-%d, constant %s, outside W %d, pairs ratio %.3f\n",
                min(na), max(na), length(unique(na)) == 1, out_of(Sa), mean(pairs_of(Sa)) / p_obs)) }
    b <- coef(KP$fit); pr <- KP$fit$par
    eta_hat <- fill_im(Reduce(`+`, Map(function(a, nmc) a * unname(b[nmc]), covs, names(covs))) + unname(b["(Intercept)"]))
    set.seed(SEED2 + 30)
    Sb <- try_("direct rLGCP n.cond from fitted trend", rLGCP(model = "exp", mu = eta_hat - unname(pr["sigma2"]) / 2,
                                                              param = list(var = unname(pr["sigma2"]), scale = unname(pr["alpha"])), win = W_km,
                                                              nsim = 20, n.cond = n_obs, drop = TRUE))
    if (!is.null(Sb)) { nb <- n_poly(Sb)
    cat(sprintf("  (b) n range %d-%d, constant %s, outside W %d, pairs ratio %.3f\n",
                min(nb), max(nb), length(unique(nb)) == 1, out_of(Sb), mean(pairs_of(Sb)) / p_obs)) }
  } else cat("  no L4 clik2 fit available\n")
  
  ## ---------------- C. Geyer expansion ----------------
  cat("\n== C. Geyer conditional simulation: default expansion vs expand = 1 (rbord 0, r 3, sat 3) ==\n")
  P8 <- as_list(sim_geyer(MULT_G3, 2.5, 3, 3, 197, 8, SEED0 + 1))
  cat("8 G3 patterns identical to scripts 26-28:",
      identical(as.numeric(vapply(P8, function(P) npoints(to_poly(P)), numeric(1))), as.numeric(N26)), "\n")
  GR <- data.frame()
  for (k in seq_len(N_PAT_GY)) {
    Xk <- to_poly(P8[[k]]); nk <- npoints(Xk); pk <- npairs3(Xk)
    fk <- quiet(geyer_fit(Xk, r = 3, sat = 3, rbord = 0))
    if (is.null(fk)) { cat(sprintf("  pattern %d: fit FAILED\n", k)); next }
    for (cn in c("default", "expand1")) {
      S <- sim_cond(fk, Xk, N_SIM, SEED2 + 1000 * k + 500 * (cn == "expand1"), cn)
      if (is.null(S)) next
      nr <- n_raw(S); np <- n_poly(S); pp <- pairs_of(S)
      GR <- rbind(GR, data.frame(k = k, ctrl = cn, n_obs = nk, n_raw_min = min(nr), n_raw_max = max(nr),
                                 n_poly_min = min(np), n_poly_max = max(np), outside = out_of(S),
                                 n_ratio = mean(np) / nk, pairs_obs = pk, pairs_mean = mean(pp), pair_ratio = mean(pp) / pk))
      cat(sprintf("  pat %d %-8s: n_obs %d | raw n %d-%d, in-W n %d-%d, points outside W %d | pairs<=3km obs %d, sim %.1f, ratio %.3f\n",
                  k, cn, nk, min(nr), max(nr), min(np), max(np), out_of(S), pk, mean(pp), mean(pp) / pk))
    }
  }
  if (nrow(GR)) {
    write.csv(GR, file.path(OUT_DIR, "29_geyer_expand_probe.csv"), row.names = FALSE)
    cat("\n-- Summary by control (n fixed in all sims = every raw n equals n_obs) --\n")
    AG <- aggregate(cbind(n_ratio, pair_ratio) ~ ctrl, data = GR, FUN = mean)
    AG$n_fixed <- vapply(AG$ctrl, function(cc) with(GR[GR$ctrl == cc, ], all(n_raw_min == n_obs & n_raw_max == n_obs)), logical(1))
    print(AG, digits = 3, row.names = FALSE)
  }
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "29_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")