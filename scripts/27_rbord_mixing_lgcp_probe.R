# =====================================================================
# 27_rbord_mixing_lgcp_probe.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : (A) compare Geyer fits at rbord 0 / 5 / 10 km on the 8 G3 truth
#               patterns of script 26: simulated-n / observed-n ratio, gamma_hat,
#               coefficient separation, default vcov usability;
#           (B) mixing test on the worst pattern: mean simulated n against
#               start value (low / observed / high) and chain length (5e5 / 2e6);
#           (C) LGCP truth generator on the full-frame mu (cells L1, L4) and
#               kppm (clik2) parameter recovery.
# Outputs : outputs/27_report.txt, 27_geyer_rbord_probe.csv,
#           27_mixing_probe.csv, 27_lgcp_probe.csv, 27_sessionInfo.txt
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           Acceptance rule (set before running): choose the rbord with the
#           smallest mean |log(sim n / obs n)| if it is <= 0.10.
#           The rLGCP mean parametrisation mu = log(lambda) - var/2 is an
#           INFERENCE; section C tests it through the mean n.
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

SEED0      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 6000   # same as script 26 -> same 8 patterns
SEED1      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 7000
NREP_PROD  <- 5e5
NREP_LONG  <- 2e6
N_PAT      <- 8
N_SIM      <- 10
RB_SET     <- c(0, 5, 10)
N_LGCP     <- 20
N_LGCP_FIT <- 10
N26        <- c(209, 174, 197, 180, 183, 190, 199, 174)   # n of the 8 patterns in script 26 (VERIFIED output)
ACCEPT_LOG <- 0.10

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
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
sim_geyer <- function(mult, gamma, r, sat, n_start, nsim, seed, nrep = NREP_PROD) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}
sim_fit <- function(fit, nsim, seed, nrep, n_start = NULL) {
  set.seed(seed)
  a <- list(nsim = nsim, drop = TRUE, verbose = FALSE, progress = FALSE, control = rmhcontrol(nrep = nrep))
  if (!is.null(n_start)) a$start <- rmhstart(n.start = n_start)
  tryCatch(suppressWarnings(do.call(simulate, c(list(fit), a))),
           error = function(e) { cat("  simulate FAILED:", conditionMessage(e), "\n"); NULL })
}
n_sims <- function(S) if (is.null(S)) NA_real_ else vapply(as_list(S), function(P) npoints(to_poly(P)), numeric(1))

REPORT <- file.path(OUT_DIR, "27_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("27_rbord_mixing_lgcp_probe | ", format(Sys.time()), "\n", sep = "")
  cat("G3 multiplier", MULT_G3, "| production nrep", NREP_PROD, "| long nrep", NREP_LONG, "\n\n")
  
  ## ---------------- 0. regenerate the 8 patterns ----------------
  cat("== 0. Regenerate 8 G3 truth patterns (same seed as script 26) ==\n")
  P8 <- as_list(sim_geyer(MULT_G3, 2.5, 3, 3, 197, N_PAT, SEED0 + 1))
  n8 <- vapply(P8, function(P) npoints(to_poly(P)), numeric(1))
  cat("n now :", paste(n8, collapse = " "), "\n")
  cat("n s.26:", paste(N26, collapse = " "), "\n")
  cat("Identical to script 26:", identical(as.numeric(n8), as.numeric(N26)), "\n")
  
  ## ---------------- A. rbord comparison ----------------
  cat("\n== A. Geyer fit at r = 3, sat = 3: rbord 0 / 5 / 10; simulate 10 patterns per fit at 5e5 steps ==\n")
  RES <- data.frame()
  for (k in seq_along(P8)) {
    Xk <- to_poly(P8[[k]])
    for (rb in RB_SET) {
      fk <- quiet(geyer_fit(Xk, r = 3, sat = 3, rbord = rb))
      if (is.null(fk)) { cat(sprintf("  pattern %d rbord %g: fit FAILED\n", k, rb)); next }
      bk <- coef(fk)
      S <- sim_fit(fk, N_SIM, SEED1 + 100 * k + rb, NREP_PROD)
      ns <- n_sims(S)
      RES <- rbind(RES, data.frame(k = k, rbord = rb, n_obs = npoints(Xk),
                                   n_active = sum(bdist.points(Xk) >= rb),
                                   na_coef = sum(is.na(bk)),
                                   max_abs_coef = max(abs(bk[-1]), na.rm = TRUE),
                                   gamma_hat = exp(unname(bk["Interaction"])),
                                   vcov_ok = vcov_ok(fk),
                                   sim_mean = mean(ns), sim_sd = sd(ns),
                                   ratio = mean(ns) / npoints(Xk)))
      cat(sprintf("  pattern %d rbord %2g: n_obs %d active %d | max|coef| %.1f gamma_hat %.2f vcov_ok %s | sim mean %.1f sd %.1f ratio %.3f\n",
                  k, rb, npoints(Xk), sum(bdist.points(Xk) >= rb), max(abs(bk[-1]), na.rm = TRUE),
                  exp(unname(bk["Interaction"])), RES$vcov_ok[nrow(RES)], mean(ns), sd(ns), mean(ns) / npoints(Xk)))
    }
  }
  write.csv(RES, file.path(OUT_DIR, "27_geyer_rbord_probe.csv"), row.names = FALSE)
  AGG <- do.call(rbind, lapply(split(RES, RES$rbord), function(d)
    data.frame(rbord = d$rbord[1], n_fits = nrow(d), mean_ratio = mean(d$ratio), sd_ratio = sd(d$ratio),
               mean_abs_log_ratio = mean(abs(log(d$ratio))), mean_gamma_hat = mean(d$gamma_hat),
               sd_gamma_hat = sd(d$gamma_hat), frac_vcov_ok = mean(d$vcov_ok),
               n_separated = sum(d$max_abs_coef > 8), mean_n_active = mean(d$n_active))))
  cat("\n-- Summary by rbord (separated = max|non-intercept coef| > 8; threshold is an INFERENCE) --\n")
  print(AGG, digits = 3, row.names = FALSE)
  best <- AGG[which.min(AGG$mean_abs_log_ratio), ]
  cat(sprintf("Acceptance rule (<= %.2f): best rbord %g, mean |log ratio| %.3f -> %s\n",
              ACCEPT_LOG, best$rbord, best$mean_abs_log_ratio,
              if (best$mean_abs_log_ratio <= ACCEPT_LOG) "ACCEPTED" else "NOT MET"))
  
  ## ---------------- B. mixing / start-value test ----------------
  cat("\n== B. Mixing test on the pattern with the largest ratio at rbord 10 ==\n")
  R10 <- RES[RES$rbord == 10, ]
  kmax <- R10$k[which.max(R10$ratio)]
  Xw <- to_poly(P8[[kmax]]); nw <- npoints(Xw)
  cat(sprintf("Worst pattern at rbord 10: k = %d, n_obs = %d, ratio %.3f\n", kmax, nw, max(R10$ratio)))
  MIXR <- data.frame()
  for (rb in c(10, 0)) {
    fw <- quiet(geyer_fit(Xw, r = 3, sat = 3, rbord = rb))
    if (is.null(fw)) next
    for (st in c(low = 100, obs = nw, high = 400)) {
      for (nr in c(NREP_PROD, NREP_LONG)) {
        t0 <- unname(proc.time()["elapsed"])
        S <- sim_fit(fw, N_SIM, SEED1 + 5000 + rb * 7 + st + round(log10(nr)), nr, n_start = st)
        ns <- n_sims(S)
        MIXR <- rbind(MIXR, data.frame(rbord = rb, start_n = st, nrep = nr, mean_sim = mean(ns),
                                       sd_sim = sd(ns), ratio = mean(ns) / nw))
        cat(sprintf("  rbord %2g start %3d nrep %.0e: sim mean %.1f sd %.1f ratio %.3f (%.0f s)\n",
                    rb, st, nr, mean(ns), sd(ns), mean(ns) / nw, unname(proc.time()["elapsed"]) - t0))
      }
    }
  }
  write.csv(MIXR, file.path(OUT_DIR, "27_mixing_probe.csv"), row.names = FALSE)
  
  ## ---------------- C. LGCP truth generator ----------------
  cat("\n== C. LGCP truth generator on the full-frame mu, and kppm recovery ==\n")
  eta_full <- Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta))
  mu0 <- eta_full + truth$b0
  cat("mu0 has NA:", anyNA(mu0$v), "\n")
  cat(sprintf("Integral of exp(mu0) over W_km: %.3f (expect %d)\n",
              integral(exp(mu0)[W_km, drop = FALSE]), RECOVERY_ENV$expected_n))
  LG <- list(L1 = c(var = 0.5, scale = 3), L4 = c(var = 1.5, scale = 8))
  LGR <- data.frame()
  for (nm in names(LG)) {
    v <- unname(LG[[nm]]["var"]); s <- unname(LG[[nm]]["scale"])
    PL <- NULL; used <- NA_character_
    for (mod in c("exp", "exponential")) {
      if (is.null(PL)) {
        set.seed(SEED1 + 9000)
        PL <- try_(sprintf("rLGCP %s model=%s", nm, mod),
                   rLGCP(model = mod, mu = mu0 - v / 2, var = v, scale = s, win = W_km,
                         nsim = N_LGCP, drop = TRUE))
        if (!is.null(PL)) used <- mod
      }
    }
    if (is.null(PL)) { cat(sprintf("  %s: generator failed with both model names\n", nm)); next }
    PLl <- lapply(as_list(PL), to_poly)
    nn <- vapply(PLl, npoints, numeric(1))
    cat(sprintf("  %s (var %.1f, scale %.0f km, model '%s'): mean n %.1f sd %.1f min %d max %d (project_state: 192-203 at n = 197)\n",
                nm, v, s, used, mean(nn), sd(nn), min(nn), max(nn)))
    nf <- min(N_LGCP_FIT, length(PLl))
    s2 <- rep(NA_real_, nf); al <- rep(NA_real_, nf)
    for (k in seq_len(nf)) {
      fk <- quiet(kppm(PLl[[k]], trend = formula_E3, clusters = "LGCP", model = "exponential",
                       covariates = covs, method = "clik2"))
      if (!is.null(fk)) { s2[k] <- unname(fk$par["sigma2"]); al[k] <- unname(fk$par["alpha"]) }
    }
    cat(sprintf("     kppm ok %d of %d | sigma2 hat median %.3f (IQR %.3f-%.3f), truth %.2f | alpha hat median %.2f (IQR %.2f-%.2f), truth %.0f | sigma2 hat < 0.01: %d\n",
                sum(!is.na(s2)), nf, median(s2, na.rm = TRUE), quantile(s2, 0.25, na.rm = TRUE),
                quantile(s2, 0.75, na.rm = TRUE), v, median(al, na.rm = TRUE),
                quantile(al, 0.25, na.rm = TRUE), quantile(al, 0.75, na.rm = TRUE), s,
                sum(s2 < 0.01, na.rm = TRUE)))
    LGR <- rbind(LGR, data.frame(cell = nm, rep = seq_len(nf), n = nn[seq_len(nf)],
                                 sigma2_hat = s2, alpha_hat = al, var_true = v, scale_true = s, model = used))
  }
  if (nrow(LGR)) write.csv(LGR, file.path(OUT_DIR, "27_lgcp_probe.csv"), row.names = FALSE)
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "27_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")