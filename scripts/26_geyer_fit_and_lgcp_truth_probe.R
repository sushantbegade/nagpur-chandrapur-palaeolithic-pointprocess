# =====================================================================
# 26_geyer_fit_and_lgcp_truth_probe.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : (A) diagnose why vcov() of the Geyer ppm fit is singular (script 25);
#           (B) frequency of NA coefficients / failed vcov over 8 G3 patterns;
#           (C) simulate.ppm(Geyer fit) mean n vs observed n at production
#               rmh length (mixing / border behaviour);
#           (D) smoke-test the LGCP truth generator (rLGCP) and kppm
#               parameter recovery for cells L1 and L4.
# Outputs : outputs/26_report.txt, 26_geyer_fit_diag.csv, 26_sim_vs_obs.csv,
#           26_lgcp_probe.csv, 26_sessionInfo.txt
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           The rLGCP mean parametrisation (mu = log(lamW) - var/2) is an
#           INFERENCE; section D tests it through the mean n.
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

SEED_PROBE  <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 5000   # same base as script 25
SEED0       <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 6000
RBORD       <- 10
NREP_PROD   <- 5e5
N_FITS      <- 8
N_MIX_FITS  <- 4
N_SIM_FIT   <- 10
N_LGCP      <- 20
N_LGCP_FIT  <- 10

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
quiet <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) NULL)
vcov_ok <- function(fit) {
  V <- quiet(vcov(fit))
  !is.null(V) && !anyNA(V) && all(is.finite(diag(V)))
}
sim_geyer <- function(mult, gamma, r, sat, n_start, nsim, seed, nrep = NREP_PROD) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)

REPORT <- file.path(OUT_DIR, "26_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("26_geyer_fit_and_lgcp_truth_probe | ", format(Sys.time()), "\n", sep = "")
  cat("rbord", RBORD, "km | production rmh nrep", NREP_PROD, "| G3 multiplier", MULT_G3, "\n\n")
  
  ## ---------------- A. script-25 pattern ----------------
  cat("== A. Script-25 pattern (G3 truth, same seed) ==\n")
  XG <- sim_geyer(MULT_G3, 2.5, 3, 3, 197, 1, SEED_PROBE + 3)
  cat("n =", npoints(XG), "(script 25 gave 177)\n")
  fitG <- geyer_fit(XG, r = 3, sat = 3, rbord = RBORD)
  fitP <- ppm(make_quad(XG), trend = formula_E3, covariates = covs)
  bG <- coef(fitG); bP <- coef(fitP)
  print(data.frame(term = names(bG), Geyer = unname(bG), Poisson = unname(bP[names(bG)])),
        digits = 4, row.names = FALSE)
  cat("NA coefficients in Geyer fit:", sum(is.na(bG)), "\n")
  
  Q   <- make_quad(XG)
  xq  <- x.quad(Q); yq <- y.quad(Q); dq <- is.data(Q)
  bdq <- bdist.points(ppp(xq, yq, window = W_km, check = FALSE))
  act <- bdq >= RBORD
  Zq  <- vapply(covs, function(im) lookup.im(im, xq, yq, naok = TRUE), numeric(length(xq)))
  Xq  <- cbind(`(Intercept)` = 1, Zq)
  okr <- act & stats::complete.cases(Xq)
  cat(sprintf("Quadrature points %d (data %d). Active at rbord %g km: %d (data %d). Active rows with NA covariates: %d\n",
              length(xq), sum(dq), RBORD, sum(act), sum(act & dq), sum(act & !stats::complete.cases(Xq))))
  qa <- qr(Xq[okr, , drop = FALSE])
  cat("Design rank over active rows:", qa$rank, "of", ncol(Xq), "\n")
  if (qa$rank < ncol(Xq))
    cat("Aliased columns:", paste(colnames(Xq)[qa$pivot[-seq_len(qa$rank)]], collapse = ", "), "\n")
  ind <- grep("^(lith_|geom_)", colnames(Xq), value = TRUE)
  print(data.frame(column = ind,
                   all_quad    = colMeans(Xq[, ind, drop = FALSE], na.rm = TRUE),
                   active_quad = colMeans(Xq[okr, ind, drop = FALSE]),
                   active_data = colMeans(Xq[okr & dq, ind, drop = FALSE])),
        digits = 3, row.names = FALSE)
  
  cat("\n-- vcov variants on the Section A Geyer fit --\n")
  for (v in list(list("default", list()), list("fine=TRUE", list(fine = TRUE)),
                 list("hessian=TRUE", list(hessian = TRUE)))) {
    Vv <- try_(paste("vcov", v[[1]]), do.call(vcov, c(list(fitG), v[[2]])))
    cat(sprintf("  vcov %-14s usable: %s\n", v[[1]], !is.null(Vv) && !anyNA(Vv) && all(is.finite(diag(Vv)))))
  }
  cat("\n-- same pattern, other rbord values --\n")
  for (rb in c(0, 5)) {
    fb <- try_(sprintf("geyer_fit rbord=%g", rb), geyer_fit(XG, r = 3, sat = 3, rbord = rb))
    if (!is.null(fb))
      cat(sprintf("  rbord %g: NA coef %d, log gamma %.3f, vcov usable %s\n",
                  rb, sum(is.na(coef(fb))), unname(tail(coef(fb), 1)), vcov_ok(fb)))
  }
  
  ## ---------------- B. 8 truth patterns ----------------
  cat("\n== B. Eight independent G3 truth patterns, fitted at truth (r = 3, sat = 3) ==\n")
  P8 <- as_list(sim_geyer(MULT_G3, 2.5, 3, 3, 197, N_FITS, SEED0 + 1))
  R8 <- data.frame()
  for (k in seq_along(P8)) {
    Xk <- to_poly(P8[[k]])
    fk <- geyer_fit(Xk, r = 3, sat = 3, rbord = RBORD)
    bk <- coef(fk)
    R8 <- rbind(R8, data.frame(k = k, n = npoints(Xk), n_active = sum(bdist.points(Xk) >= RBORD),
                               na_coef = sum(is.na(bk)), gamma_hat = exp(unname(bk["Interaction"])),
                               vcov_ok = vcov_ok(fk)))
  }
  print(R8, digits = 4, row.names = FALSE)
  write.csv(R8, file.path(OUT_DIR, "26_geyer_fit_diag.csv"), row.names = FALSE)
  
  ## ---------------- C. simulate.ppm vs observed n ----------------
  cat("\n== C. simulate.ppm from the fitted Geyer model: mean n vs observed n ==\n")
  MIX <- data.frame()
  for (k in seq_len(N_MIX_FITS)) {
    Xk <- to_poly(P8[[k]]); fk <- geyer_fit(Xk, r = 3, sat = 3, rbord = RBORD)
    for (nr in if (k == 1) c(1e4, 1e5, NREP_PROD) else NREP_PROD) {
      t0 <- unname(proc.time()["elapsed"])
      set.seed(SEED0 + 100 * k + round(log10(nr)))
      S <- quiet(simulate(fk, nsim = N_SIM_FIT, drop = TRUE, verbose = FALSE, progress = FALSE,
                          control = rmhcontrol(nrep = nr)))
      ns <- if (is.null(S)) NA_real_ else vapply(as_list(S), function(P) npoints(to_poly(P)), numeric(1))
      MIX <- rbind(MIX, data.frame(fit = k, nrep = nr, n_obs = npoints(Xk), mean_sim = mean(ns),
                                   sd_sim = sd(ns), ratio = mean(ns) / npoints(Xk),
                                   sec = unname(proc.time()["elapsed"]) - t0))
      cat(sprintf("  fit %d nrep %.0e: n_obs %d, sim mean %.1f sd %.1f ratio %.3f (%.0f s)\n",
                  k, nr, npoints(Xk), mean(ns), sd(ns), mean(ns) / npoints(Xk),
                  unname(proc.time()["elapsed"]) - t0))
    }
  }
  p5 <- MIX[MIX$nrep == NREP_PROD, ]
  cat(sprintf("Pooled ratio at production length: mean %.3f, se %.3f (n fits %d)\n",
              mean(p5$ratio, na.rm = TRUE), sd(p5$ratio, na.rm = TRUE) / sqrt(sum(!is.na(p5$ratio))), nrow(p5)))
  write.csv(MIX, file.path(OUT_DIR, "26_sim_vs_obs.csv"), row.names = FALSE)
  
  ## ---------------- D. LGCP truth generator ----------------
  cat("\n== D. LGCP truth generator (rLGCP) and kppm recovery ==\n")
  cat("rLGCP formals:", paste(names(formals(rLGCP)), collapse = ", "), "\n")
  LG <- list(L1 = c(var = 0.5, scale = 3), L4 = c(var = 1.5, scale = 8))
  LGR <- data.frame()
  for (nm in names(LG)) {
    v <- unname(LG[[nm]]["var"]); s <- unname(LG[[nm]]["scale"])
    PL <- NULL; used <- NA_character_
    for (mod in c("exp", "exponential")) {
      if (is.null(PL)) {
        set.seed(SEED0 + 200)
        PL <- try_(sprintf("rLGCP %s model=%s", nm, mod),
                   rLGCP(model = mod, mu = log(truth$lamW) - v / 2, var = v, scale = s,
                         win = W_km, nsim = N_LGCP, drop = TRUE))
        if (!is.null(PL)) used <- mod
      }
    }
    if (is.null(PL)) { cat(sprintf("  %s: generator failed with both model names\n", nm)); next }
    PLl <- lapply(as_list(PL), to_poly)
    nn  <- vapply(PLl, npoints, numeric(1))
    cat(sprintf("  %s (var %.1f, scale %.0f km, model '%s'): mean n %.1f sd %.1f min %d max %d (project_state: 192-203 at n = 197)\n",
                nm, v, s, used, mean(nn), sd(nn), min(nn), max(nn)))
    s2 <- al <- rep(NA_real_, min(N_LGCP_FIT, length(PLl)))
    for (k in seq_along(s2)) {
      fk <- quiet(kppm(PLl[[k]], trend = formula_E3, clusters = "LGCP", model = "exponential",
                       covariates = covs, method = "clik2"))
      if (!is.null(fk)) { s2[k] <- unname(fk$par["sigma2"]); al[k] <- unname(fk$par["alpha"]) }
    }
    cat(sprintf("     kppm fits ok %d of %d | sigma2 hat: median %.3f (IQR %.3f-%.3f), truth %.2f | alpha hat: median %.2f km (IQR %.2f-%.2f), truth %.0f\n",
                sum(!is.na(s2)), length(s2), median(s2, na.rm = TRUE),
                quantile(s2, 0.25, na.rm = TRUE), quantile(s2, 0.75, na.rm = TRUE), v,
                median(al, na.rm = TRUE), quantile(al, 0.25, na.rm = TRUE), quantile(al, 0.75, na.rm = TRUE), s))
    LGR <- rbind(LGR, data.frame(cell = nm, rep = seq_along(s2), n = nn[seq_along(s2)],
                                 sigma2_hat = s2, alpha_hat = al, var_true = v, scale_true = s, model = used))
  }
  if (nrow(LGR)) write.csv(LGR, file.path(OUT_DIR, "26_lgcp_probe.csv"), row.names = FALSE)
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "26_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")