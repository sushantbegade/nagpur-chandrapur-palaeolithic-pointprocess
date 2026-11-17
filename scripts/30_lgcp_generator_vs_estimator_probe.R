# =====================================================================
# 30_lgcp_generator_vs_estimator_probe.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : (A) separate LGCP generator from clik2 estimator: realised field
#               mean/variance; pair-correlation at 4x intensity vs theory
#               exp(var * exp(-r/scale)) using the TRUE intensity; clik2 with
#               rmax = default / 10 / 20 km on the same patterns;
#           (B) conditional-on-n LGCP via saved Lambda + rpoint(): n constancy,
#               points outside W, timing, pairs per n^2 vs truth patterns;
#           (C) Geyer conditional simulation: expand = 1 alone vs
#               expand = 1 + window = W_km (8 patterns).
# Outputs : outputs/30_report.txt, 30_lgcp_diag.csv, 30_lgcp_ncond_probe.csv,
#           30_geyer_window_probe.csv, 30_sessionInfo.txt
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           rLGCP saveLambda attribute, rpoint(f = im, forcewin) and
#           simulate.ppm(window = ) behaviour are UNVERIFIED against
#           documentation; this script tests them. Warnings are tallied at the end.
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

SEED0     <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 6000    # same 8 G3 patterns as scripts 26-29
SEED3     <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 15000
NREP_PROD <- 5e5
HI_MULT   <- 4
N_HI      <- 6
N_GEN     <- 3
N_TRUTH   <- 6
N_NCOND   <- 20
N_PAT_GY  <- 8
N_SIM     <- 15
N26       <- c(209, 174, 197, 180, 183, 190, 199, 174)   # VERIFIED n of the 8 patterns
RG        <- seq(0, 15, by = 0.5)
R_SHOW    <- c(1, 2, 3, 5, 8, 12)
BW_PCF    <- 0.5
PCF_ARGS  <- list(kernel = "epanechnikov", zerocor = "convolution", divisor = "r", renormalise = FALSE)

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULT_G3 <- gm$mult[gm$cell == "G3"]; stopifnot(length(MULT_G3) == 1)

# ---- helpers ----
WARN <- new.env(); WARN$msgs <- character(0)
note_w  <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
try_ <- function(label, expr) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch(withCallingHandlers(expr, warning = note_w),
                  error = function(e) { cat(sprintf("[%s] FAILED: %s\n", label, conditionMessage(e))); NULL })
  cat(sprintf("[%s] %.1f s\n", label, unname(proc.time()["elapsed"]) - t0))
  out
}
quiet   <- function(expr) tryCatch(withCallingHandlers(expr, warning = note_w), error = function(e) NULL)
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
fit_lgcp2 <- function(X, rmax) {
  a <- list(X, trend = formula_E3, clusters = "LGCP", model = "exponential", covariates = covs, method = "clik2")
  if (!is.null(rmax)) a$rmax <- rmax
  quiet(do.call(kppm, a))
}
pars_of <- function(f) if (is.null(f)) c(NA_real_, NA_real_) else unname(c(f$par["sigma2"], f$par["alpha"]))
summ <- function(label, s2, al, v, s)
  cat(sprintf("  %-26s ok %d/%d | sigma2 median %.3f (IQR %.3f-%.3f) truth %.2f | alpha median %.2f (IQR %.2f-%.2f) truth %.0f\n",
              label, sum(!is.na(s2)), length(s2), median(s2, na.rm = TRUE), qs(s2)[1], qs(s2)[2], v,
              median(al, na.rm = TRUE), qs(al)[1], qs(al)[2], s))
pcf_show <- function(P, lam) tryCatch({
  gh <- do.call(pcfinhom, c(list(to_poly(P), lambda = lam, r = RG, correction = "Ripley", bw = BW_PCF), PCF_ARGS))
  col_of(gh, c("iso", "trans"))[match(round(R_SHOW, 6), round(gh$r, 6))]
}, error = function(e) rep(NA_real_, length(R_SHOW)))
npairs3  <- function(P) sum(closepaircounts(P, 3)) / 2
n_raw    <- function(S) vapply(as_list(S), npoints, numeric(1))
n_poly   <- function(S) vapply(as_list(S), function(P) npoints(to_poly(P)), numeric(1))
pairs_of <- function(S) vapply(as_list(S), function(P) npairs3(to_poly(P)), numeric(1))
out_of   <- function(S) sum(vapply(as_list(S), function(P) npoints(P) - sum(inside.owin(P$x, P$y, W_km)), numeric(1)))
pn2_of   <- function(S) vapply(as_list(S), function(P) { Q <- to_poly(P); npairs3(Q) / npoints(Q)^2 }, numeric(1))
sim_cond2 <- function(fit, X, nsim, seed, mode) {
  set.seed(seed)
  a <- list(fit, nsim = nsim, drop = TRUE, verbose = FALSE, progress = FALSE,
            control = rmhcontrol(nrep = NREP_PROD, p = 1, expand = 1),
            start = rmhstart(n.start = npoints(X)))
  if (mode == "window") a$window <- W_km
  tryCatch(withCallingHandlers(do.call(simulate, a), warning = note_w),
           error = function(e) { cat("  simulate FAILED (", mode, "): ", conditionMessage(e), "\n", sep = ""); NULL })
}

REPORT <- file.path(OUT_DIR, "30_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("30_lgcp_generator_vs_estimator_probe | ", format(Sys.time()), "\n", sep = "")
  cat("G3 multiplier", MULT_G3, "| production nrep", NREP_PROD, "\n\n")
  
  mu0 <- Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta)) + truth$b0
  mu_fill <- fill_im(mu0)
  lam4 <- exp(mu_fill + log(HI_MULT))          # true expected intensity at x4
  LG <- list(L1 = c(var = 0.5, scale = 3), L4 = c(var = 1.5, scale = 8))
  
  ## ---------------- A. generator vs estimator ----------------
  cat("== A. Generator vs estimator ==\n")
  PCFR <- data.frame(); CLR <- data.frame()
  for (nm in names(LG)) {
    v <- unname(LG[[nm]]["var"]); s <- unname(LG[[nm]]["scale"])
    cat(sprintf("\n-- %s: var %.1f, scale %.0f km --\n", nm, v, s))
    
    ## A1 saved field
    P1 <- try_(sprintf("rLGCP %s %d sims saveLambda", nm, N_GEN), gen_lgcp(v, s, 0, N_GEN, SEED3 + 10, TRUE))
    if (!is.null(P1)) {
      Pl <- as_list(P1)
      cat("  attributes of sim 1:", paste(names(attributes(Pl[[1]])), collapse = ", "), "\n")
      cat("  class of attr 'Lambda':", paste(class(attr(Pl[[1]], "Lambda")), collapse = ","), "\n")
      uu <- t(vapply(Pl, function(P) tryCatch({
        u <- (log(attr(P, "Lambda")) - mu_fill)[W_km, drop = TRUE]
        c(mean(u), var(u))
      }, error = function(e) { cat("  field check FAILED:", conditionMessage(e), "\n"); c(NA_real_, NA_real_) }), numeric(2)))
      cat(sprintf("  realised field u = log Lambda - mu inside W: mean %.3f (expect %.3f), pixel variance %.3f (expect about %.2f)\n",
                  mean(uu[, 1], na.rm = TRUE), -v / 2, mean(uu[, 2], na.rm = TRUE), v))
    }
    
    ## A2 pair-correlation at x4 vs theory
    PH <- try_(sprintf("rLGCP %s x%d, %d sims", nm, HI_MULT, N_HI), gen_lgcp(v, s, log(HI_MULT), N_HI, SEED3 + 20, FALSE))
    if (!is.null(PH)) {
      PHl <- as_list(PH)
      cat(sprintf("  x%d mean n %.1f\n", HI_MULT, mean(n_poly(PHl))))
      G <- t(vapply(PHl, function(P) pcf_show(P, lam4), numeric(length(R_SHOW))))
      th <- exp(v * exp(-R_SHOW / s))
      cat("  pcfinhom with TRUE intensity (renormalise = FALSE), pooled over sims:\n")
      print(data.frame(r_km = R_SHOW, theory = th, obs_mean = colMeans(G, na.rm = TRUE),
                       obs_sd = apply(G, 2, sd, na.rm = TRUE), ratio = colMeans(G, na.rm = TRUE) / th),
            digits = 3, row.names = FALSE)
      PCFR <- rbind(PCFR, data.frame(cell = nm, r_km = R_SHOW, theory = th, obs_mean = colMeans(G, na.rm = TRUE),
                                     obs_sd = apply(G, 2, sd, na.rm = TRUE)))
      
      ## A3 clik2 rmax variants on the same x4 patterns
      PHp <- lapply(PHl, to_poly)
      for (rmx in list(NULL, 10, 20)) {
        lab <- if (is.null(rmx)) "default" else as.character(rmx)
        s2 <- rep(NA_real_, length(PHp)); al <- s2
        for (k in seq_along(PHp)) { pr <- pars_of(fit_lgcp2(PHp[[k]], rmx)); s2[k] <- pr[1]; al[k] <- pr[2] }
        summ(sprintf("clik2 x%d rmax=%s", HI_MULT, lab), s2, al, v, s)
        CLR <- rbind(CLR, data.frame(cell = nm, setting = paste0("clik2_x", HI_MULT, "_rmax_", lab),
                                     rep = seq_along(s2), sigma2_hat = s2, alpha_hat = al, var_true = v, scale_true = s))
      }
    }
  }
  write.csv(rbind(transform(PCFR, kind = "pcf", setting = NA_character_, rep = NA_integer_,
                            sigma2_hat = NA_real_, alpha_hat = NA_real_, var_true = NA_real_, scale_true = NA_real_),
                  transform(CLR, kind = "clik2", r_km = NA_real_, theory = NA_real_, obs_mean = NA_real_, obs_sd = NA_real_))[0, 0],
            file.path(OUT_DIR, "30_lgcp_diag.csv"), row.names = FALSE)   # placeholder header overwritten below
  if (nrow(PCFR)) write.csv(PCFR, file.path(OUT_DIR, "30_lgcp_diag.csv"), row.names = FALSE)
  if (nrow(CLR))  write.csv(CLR, file.path(OUT_DIR, "30_lgcp_diag_clik2.csv"), row.names = FALSE)
  
  ## ---------------- B. conditional LGCP ----------------
  cat("\n== B. Conditional-on-n LGCP: saved Lambda + rpoint ==\n")
  vB <- 1.5; sB <- 8
  PB <- try_("rLGCP L4 x1 truth patterns", gen_lgcp(vB, sB, 0, N_TRUTH, SEED3 + 30, FALSE))
  KP <- NULL
  if (!is.null(PB)) {
    PBl <- lapply(as_list(PB), to_poly)
    pn2_truth <- vapply(PBl, function(P) npairs3(P) / npoints(P)^2, numeric(1))
    for (k in seq_along(PBl)) { fk <- fit_lgcp2(PBl[[k]], NULL); if (!is.null(fk)) { KP <- list(fit = fk, X = PBl[[k]]); break } }
    cat(sprintf("  truth patterns (L4): pairs<=3km / n^2 x1e4: mean %.3f (range %.3f-%.3f)\n",
                1e4 * mean(pn2_truth), 1e4 * min(pn2_truth), 1e4 * max(pn2_truth)))
  }
  if (!is.null(KP)) {
    n_obs <- npoints(KP$X); p_obs <- npairs3(KP$X)
    b <- coef(KP$fit); pr <- KP$fit$par
    eta_hat <- fill_im(Reduce(`+`, Map(function(a, nmc) a * unname(b[nmc]), covs, names(covs))) + unname(b["(Intercept)"]))
    cat(sprintf("  fitted pattern: n %d, pairs<=3km %d, sigma2_hat %.3f, alpha_hat %.2f (truth %.1f, %.0f)\n",
                n_obs, p_obs, pr["sigma2"], pr["alpha"], vB, sB))
    cat(sprintf("  observed pairs/n^2 x1e4: %.3f\n", 1e4 * p_obs / n_obs^2))
    
    set.seed(SEED3 + 40)
    Sd <- try_("(a) direct rLGCP n.cond, 5 sims", rLGCP(model = "exp", mu = eta_hat - unname(pr["sigma2"]) / 2,
                                                        param = list(var = unname(pr["sigma2"]), scale = unname(pr["alpha"])), win = W_km,
                                                        nsim = 5, n.cond = n_obs, drop = TRUE))
    if (!is.null(Sd)) cat(sprintf("  (a) n %d-%d, outside W %d, pairs/n^2 x1e4 %.3f\n",
                                  min(n_raw(Sd)), max(n_raw(Sd)), out_of(Sd), 1e4 * mean(pn2_of(Sd))))
    
    set.seed(SEED3 + 50)
    SU <- try_(sprintf("(b1) unconditional rLGCP saveLambda, %d sims", N_NCOND),
               rLGCP(model = "exp", mu = eta_hat - unname(pr["sigma2"]) / 2,
                     param = list(var = unname(pr["sigma2"]), scale = unname(pr["alpha"])), win = W_km,
                     nsim = N_NCOND, saveLambda = TRUE, drop = TRUE))
    if (!is.null(SU)) {
      SUl <- as_list(SU)
      cat(sprintf("  (b1) unconditional: n mean %.1f sd %.1f, pairs/n^2 x1e4 %.3f\n",
                  mean(n_poly(SUl)), sd(n_poly(SUl)), 1e4 * mean(pn2_of(SUl))))
      set.seed(SEED3 + 60)
      R1 <- try_("(b2) rpoint test, 1 sim", rpoint(n_obs, f = attr(SUl[[1]], "Lambda"), win = W_km, forcewin = TRUE))
      if (!is.null(R1)) {
        set.seed(SEED3 + 61)
        RP <- try_(sprintf("(b2) rpoint conditional, %d sims", N_NCOND),
                   lapply(SUl, function(P) rpoint(n_obs, f = attr(P, "Lambda"), win = W_km, forcewin = TRUE)))
        if (!is.null(RP)) {
          cat(sprintf("  (b2) conditional: n %d-%d, outside W %d, window type %s, pairs/n^2 x1e4 %.3f\n",
                      min(n_raw(RP)), max(n_raw(RP)), out_of(RP), Window(RP[[1]])$type, 1e4 * mean(pn2_of(RP))))
          write.csv(data.frame(rep = seq_along(RP), n = n_raw(RP), outside = vapply(RP, function(P) npoints(P) - sum(inside.owin(P$x, P$y, W_km)), numeric(1)),
                               pn2 = pn2_of(RP), n_obs = n_obs, p_obs = p_obs),
                    file.path(OUT_DIR, "30_lgcp_ncond_probe.csv"), row.names = FALSE)
        }
      }
    }
  } else cat("  no L4 clik2 fit available\n")
  
  ## ---------------- C. Geyer window ----------------
  cat("\n== C. Geyer conditional simulation: expand = 1 vs expand = 1 + window = W_km (rbord 0, r 3, sat 3) ==\n")
  P8 <- as_list(sim_geyer(MULT_G3, 2.5, 3, 3, 197, 8, SEED0 + 1))
  cat("8 G3 patterns identical to scripts 26-29:",
      identical(as.numeric(vapply(P8, function(P) npoints(to_poly(P)), numeric(1))), as.numeric(N26)), "\n")
  GR <- data.frame()
  for (k in seq_len(N_PAT_GY)) {
    Xk <- to_poly(P8[[k]]); nk <- npoints(Xk); pk <- npairs3(Xk)
    fk <- quiet(geyer_fit(Xk, r = 3, sat = 3, rbord = 0))
    if (is.null(fk)) { cat(sprintf("  pattern %d: fit FAILED\n", k)); next }
    for (md in c("expand1", "window")) {
      S <- sim_cond2(fk, Xk, N_SIM, SEED3 + 1000 * k + 500 * (md == "window"), md)
      if (is.null(S)) next
      nr <- n_raw(S); np <- n_poly(S); pp <- pairs_of(S)
      GR <- rbind(GR, data.frame(k = k, mode = md, n_obs = nk, n_raw_min = min(nr), n_raw_max = max(nr),
                                 n_poly_min = min(np), n_poly_max = max(np), outside = out_of(S),
                                 n_ratio = mean(np) / nk, pairs_obs = pk, pairs_mean = mean(pp), pair_ratio = mean(pp) / pk))
      cat(sprintf("  pat %d %-8s: n_obs %d | raw n %d-%d, in-W n %d-%d, outside W %d | pairs<=3km obs %d, sim %.1f, ratio %.3f\n",
                  k, md, nk, min(nr), max(nr), min(np), max(np), out_of(S), pk, mean(pp), mean(pp) / pk))
    }
  }
  if (nrow(GR)) {
    write.csv(GR, file.path(OUT_DIR, "30_geyer_window_probe.csv"), row.names = FALSE)
    cat("\n-- Summary by mode --\n")
    AG <- do.call(rbind, lapply(split(GR, GR$mode), function(d)
      data.frame(mode = d$mode[1], n_patterns = nrow(d), n_ratio = mean(d$n_ratio), pair_ratio = mean(d$pair_ratio),
                 sd_pair_ratio = sd(d$pair_ratio), outside_total = sum(d$outside),
                 raw_n_fixed = all(d$n_raw_min == d$n_obs & d$n_raw_max == d$n_obs))))
    print(AG, digits = 3, row.names = FALSE)
  }
  
  ## ---------------- warnings tally ----------------
  cat("\n== Warning tally (all sections) ==\n")
  if (length(WARN$msgs)) {
    tw <- sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE)
    print(head(tw, 10))
  } else cat("none\n")
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "30_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")