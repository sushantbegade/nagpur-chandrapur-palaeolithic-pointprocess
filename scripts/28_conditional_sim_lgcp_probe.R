# =====================================================================
# 28_conditional_sim_lgcp_probe.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : (A) LGCP truth generator on a NA-filled full-frame mu (cells L1, L4):
#               mean n, kppm (clik2) recovery of sigma2 / alpha, and
#               simulate.kppm(n.cond) conditional simulation;
#           (B) conditional-on-n Geyer simulation (rmhcontrol p = 1) at rbord 0
#               on 4 G3 truth patterns: n constant, inside W, pair-count
#               agreement, and mixing (uniform start vs data start; 5e5 vs 2e6).
# Outputs : outputs/28_report.txt, 28_lgcp_probe.csv, 28_geyer_cond_probe.csv,
#           28_sessionInfo.txt
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           rLGCP argument form param = list(var, scale) and the fixed-n
#           behaviour of rmhcontrol(p = 1) / simulate.kppm(n.cond) are
#           UNVERIFIED against documentation; this script tests them.
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

SEED0      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 6000   # same as scripts 26/27 -> same 8 patterns
SEED1      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 8000
NREP_PROD  <- 5e5
NREP_LONG  <- 2e6
RBORD      <- 0
N_PAT_GY   <- 4
N_SIM      <- 10
N_LGCP     <- 20
N_LGCP_FIT <- 10
N26        <- c(209, 174, 197, 180, 183, 190, 199, 174)   # VERIFIED n of the 8 patterns (scripts 26, 27)

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULT_G3 <- gm$mult[gm$cell == "G3"]; stopifnot(length(MULT_G3) == 1)

# ---- helpers (no assignments inside function bodies) ----
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
sim_geyer <- function(mult, gamma, r, sat, n_start, nsim, seed, nrep = NREP_PROD) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}
sim_cond <- function(fit, X, nsim, seed, nrep, start_type) {
  set.seed(seed)
  tryCatch(suppressWarnings(simulate(fit, nsim = nsim, drop = TRUE, verbose = FALSE, progress = FALSE,
                                     control = rmhcontrol(nrep = nrep, p = 1),
                                     start = if (start_type == "uniform") rmhstart(n.start = npoints(X)) else rmhstart(x.start = X))),
           error = function(e) { cat("  simulate FAILED:", conditionMessage(e), "\n"); NULL })
}
npairs3  <- function(P) sum(closepaircounts(P, 3)) / 2
n_of     <- function(S) vapply(as_list(S), function(P) npoints(to_poly(P)), numeric(1))
pairs_of <- function(S) vapply(as_list(S), function(P) npairs3(to_poly(P)), numeric(1))
in_of    <- function(S) vapply(as_list(S), function(P) sum(inside.owin(P$x, P$y, W_km)) - npoints(P), numeric(1))

REPORT <- file.path(OUT_DIR, "28_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("28_conditional_sim_lgcp_probe | ", format(Sys.time()), "\n", sep = "")
  cat("G3 multiplier", MULT_G3, "| rbord", RBORD, "| production nrep", NREP_PROD, "| long nrep", NREP_LONG, "\n")
  cat("closepaircounts formals:", paste(names(formals(closepaircounts)), collapse = ", "), "\n\n")
  
  ## ---------------- 0. regenerate patterns ----------------
  cat("== 0. Regenerate 8 G3 truth patterns (same seed as scripts 26/27) ==\n")
  P8 <- as_list(sim_geyer(MULT_G3, 2.5, 3, 3, 197, 8, SEED0 + 1))
  n8 <- vapply(P8, function(P) npoints(to_poly(P)), numeric(1))
  cat("Identical to script 26/27:", identical(as.numeric(n8), as.numeric(N26)), "\n")
  
  ## ---------------- A. LGCP ----------------
  cat("\n== A. LGCP truth generator on NA-filled full-frame mu ==\n")
  eta_full <- Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta))
  mu0 <- eta_full + truth$b0
  na_pix <- is.na(mu0$v)
  cat(sprintf("mu0 pixels: %d, NA: %d (%.1f%%)\n", length(na_pix), sum(na_pix), 100 * mean(na_pix)))
  cat(sprintf("Integral exp(mu0) over W_km: %.3f (expect %d)\n",
              integral(exp(mu0)[W_km, drop = FALSE]), RECOVERY_ENV$expected_n))
  fillv <- median(mu0$v[!na_pix])
  mat   <- mu0$v
  mat[is.na(mat)] <- fillv
  mu_fill <- im(mat, xcol = mu0$xcol, yrow = mu0$yrow)
  cat("mu_fill: NA left", anyNA(mu_fill$v), "| Frame(mu_fill) covers Frame(W_km):",
      is.subset.owin(Frame(W_km), Frame(mu_fill)), "\n")
  cat("rLGCP formals:", paste(names(formals(rLGCP)), collapse = ", "), "\n")
  
  LG <- list(L1 = c(var = 0.5, scale = 3), L4 = c(var = 1.5, scale = 8))
  LGR <- data.frame(); KP <- NULL
  for (nm in names(LG)) {
    v <- unname(LG[[nm]]["var"]); s <- unname(LG[[nm]]["scale"])
    PL <- NULL; used <- NA_character_
    for (wn in c("W_km", "frame")) {
      if (is.null(PL)) {
        WIN <- if (wn == "W_km") W_km else Frame(mu_fill)
        set.seed(SEED1 + 1)
        t1 <- try_(sprintf("rLGCP %s 1 sim, win=%s (timing)", nm, wn),
                   rLGCP(model = "exp", mu = mu_fill - v / 2, param = list(var = v, scale = s),
                         win = WIN, nsim = 1, drop = TRUE))
        if (!is.null(t1)) {
          set.seed(SEED1 + 2)
          PL <- try_(sprintf("rLGCP %s %d sims, win=%s", nm, N_LGCP, wn),
                     rLGCP(model = "exp", mu = mu_fill - v / 2, param = list(var = v, scale = s),
                           win = WIN, nsim = N_LGCP, drop = TRUE))
          if (!is.null(PL)) used <- wn
        }
      }
    }
    if (is.null(PL)) { cat(sprintf("  %s: generator failed in both windows\n", nm)); next }
    PLl <- lapply(as_list(PL), to_poly)
    nn  <- vapply(PLl, npoints, numeric(1))
    cat(sprintf("  %s (var %.1f, scale %.0f km, win '%s'): mean n %.1f sd %.1f min %d max %d (project_state: 192-203 at n = 197)\n",
                nm, v, s, used, mean(nn), sd(nn), min(nn), max(nn)))
    nf <- min(N_LGCP_FIT, length(PLl))
    s2 <- rep(NA_real_, nf); al <- rep(NA_real_, nf)
    for (k in seq_len(nf)) {
      fk <- quiet(kppm(PLl[[k]], trend = formula_E3, clusters = "LGCP", model = "exponential",
                       covariates = covs, method = "clik2"))
      if (!is.null(fk)) { s2[k] <- unname(fk$par["sigma2"]); al[k] <- unname(fk$par["alpha"]) }
      if (!is.null(fk) && is.null(KP) && nm == "L4") KP <- list(fit = fk, X = PLl[[k]])
    }
    cat(sprintf("     kppm ok %d of %d | sigma2 hat median %.3f (IQR %.3f-%.3f), truth %.2f | alpha hat median %.2f (IQR %.2f-%.2f), truth %.0f | sigma2 hat < 0.01: %d\n",
                sum(!is.na(s2)), nf, median(s2, na.rm = TRUE), quantile(s2, 0.25, na.rm = TRUE),
                quantile(s2, 0.75, na.rm = TRUE), v, median(al, na.rm = TRUE),
                quantile(al, 0.25, na.rm = TRUE), quantile(al, 0.75, na.rm = TRUE), s,
                sum(s2 < 0.01, na.rm = TRUE)))
    LGR <- rbind(LGR, data.frame(cell = nm, rep = seq_len(nf), n = nn[seq_len(nf)],
                                 sigma2_hat = s2, alpha_hat = al, var_true = v, scale_true = s, win = used))
  }
  if (nrow(LGR)) write.csv(LGR, file.path(OUT_DIR, "28_lgcp_probe.csv"), row.names = FALSE)
  
  cat("\n-- simulate.kppm unconditional vs n.cond (first L4 kppm fit) --\n")
  if (!is.null(KP)) {
    n_obs <- npoints(KP$X)
    Su <- try_("simulate(kppm) unconditional", simulate(KP$fit, nsim = 20, drop = TRUE, verbose = FALSE))
    if (!is.null(Su)) { nu <- n_of(Su); cat(sprintf("  unconditional: n_obs %d, sim mean %.1f sd %.1f\n", n_obs, mean(nu), sd(nu))) }
    Sc <- try_("simulate(kppm) n.cond", simulate(KP$fit, nsim = 20, n.cond = n_obs, drop = TRUE, verbose = FALSE))
    if (!is.null(Sc)) {
      nc <- n_of(Sc)
      cat(sprintf("  n.cond = %d: sim n range %d-%d (constant: %s); points outside W: %d\n",
                  n_obs, min(nc), max(nc), length(unique(nc)) == 1, sum(in_of(Sc))))
    }
  } else cat("  no kppm fit available for the n.cond test\n")
  
  ## ---------------- B. conditional Geyer ----------------
  cat("\n== B. Conditional-on-n Geyer simulation (rmhcontrol p = 1), fit at r = 3, sat = 3, rbord 0 ==\n")
  GR <- data.frame()
  for (k in seq_len(N_PAT_GY)) {
    Xk <- to_poly(P8[[k]]); nk <- npoints(Xk); pk <- npairs3(Xk)
    fk <- quiet(geyer_fit(Xk, r = 3, sat = 3, rbord = RBORD))
    if (is.null(fk)) { cat(sprintf("  pattern %d: fit FAILED\n", k)); next }
    for (st in c("uniform", "data")) {
      for (nr in c(NREP_PROD, NREP_LONG)) {
        t0 <- unname(proc.time()["elapsed"])
        S  <- sim_cond(fk, Xk, N_SIM, SEED1 + 1000 * k + (if (st == "data") 500 else 0) + round(log10(nr)), nr, st)
        if (is.null(S)) next
        ns <- n_of(S); pp <- pairs_of(S)
        GR <- rbind(GR, data.frame(k = k, start = st, nrep = nr, n_obs = nk, n_min = min(ns), n_max = max(ns),
                                   outside = sum(in_of(S)), pairs_obs = pk, pairs_mean = mean(pp),
                                   pairs_sd = sd(pp), pair_ratio = mean(pp) / pk,
                                   sec = unname(proc.time()["elapsed"]) - t0))
        cat(sprintf("  pat %d start %-7s nrep %.0e: n_obs %d, sim n %d-%d | pairs<=3km obs %d, sim mean %.1f sd %.1f ratio %.3f (%.0f s)\n",
                    k, st, nr, nk, min(ns), max(ns), pk, mean(pp), sd(pp), mean(pp) / pk,
                    unname(proc.time()["elapsed"]) - t0))
      }
    }
  }
  if (nrow(GR)) {
    write.csv(GR, file.path(OUT_DIR, "28_geyer_cond_probe.csv"), row.names = FALSE)
    cat("\n-- Summary by start x nrep (mean pair ratio over patterns; n fixed in all sims:",
        all(GR$n_min == GR$n_obs & GR$n_max == GR$n_obs), ") --\n")
    AG <- aggregate(pair_ratio ~ start + nrep, data = GR, FUN = mean)
    print(AG, digits = 3, row.names = FALSE)
  }
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "28_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")