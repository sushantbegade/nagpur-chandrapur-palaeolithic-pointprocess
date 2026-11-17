# =====================================================================
# 31_generator_control_and_geyer_decomposition.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : (A) LGCP generator check on ALIGNED grids (field mean/variance,
#               mean of Lambda vs exp(mu)); pcfinhom at 4x intensity with the
#               TRUE intensity for a POISSON CONTROL (expect ratio 1) and for
#               L1, L4 (theory exp(var * exp(-r/scale)));
#           (B) Geyer conditional simulation from the FITTED model vs the TRUE
#               model on the 8 G3 patterns: simulation window, outside-W
#               points, pairs<=3km / n^2 vs the observed pattern.
# Outputs : outputs/31_report.txt, 31_lgcp_gen.csv, 31_geyer_decomp.csv,
#           31_sessionInfo.txt
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           Whole run is under withCallingHandlers: all warnings are tallied.
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

SEED0     <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 6000    # same 8 G3 patterns as scripts 26-30
SEED4     <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 17000
NREP_PROD <- 5e5
HI_MULT   <- 4
N_FIELD   <- 10
N_PCF     <- 20
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
  out <- tryCatch(expr, error = function(e) { cat(sprintf("[%s] FAILED: %s\n", label, conditionMessage(e))); NULL })
  cat(sprintf("[%s] %.1f s\n", label, unname(proc.time()["elapsed"]) - t0))
  out
}
quiet   <- function(expr) tryCatch(expr, error = function(e) NULL)
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
fill_im <- function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
sim_geyer <- function(mult, gamma, r, sat, n_start, nsim, seed, nrep = NREP_PROD) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}
sim_geyer_cond <- function(mult, gamma, r, sat, n, nsim, seed) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n), control = rmhcontrol(nrep = NREP_PROD, p = 1, expand = 1),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}
sim_fit_cond <- function(fit, n, nsim, seed) {
  set.seed(seed)
  simulate(fit, nsim = nsim, drop = TRUE, verbose = FALSE, progress = FALSE,
           control = rmhcontrol(nrep = NREP_PROD, p = 1, expand = 1), start = rmhstart(n.start = n))
}
gen_lgcp <- function(v, s, shift, nsim, seed, save) {
  set.seed(seed)
  rLGCP(model = "exp", mu = mu_fill + shift - v / 2, param = list(var = v, scale = s),
        win = W_km, nsim = nsim, saveLambda = save, drop = TRUE)
}
field_stats <- function(P, v) tryCatch({
  Lam <- attr(P, "Lambda"); df <- as.data.frame(Lam)
  keep <- inside.owin(df$x, df$y, W_km)
  lam  <- df[[3]][keep]
  muv  <- lookup.im(mu_fill, df$x[keep], df$y[keep], naok = TRUE)
  u    <- log(lam) - muv
  c(mean(u, na.rm = TRUE), var(u, na.rm = TRUE), mean(lam, na.rm = TRUE) / mean(exp(muv), na.rm = TRUE))
}, error = function(e) { cat("  field check FAILED:", conditionMessage(e), "\n"); c(NA_real_, NA_real_, NA_real_) })
pcf_show <- function(P, lam) tryCatch({
  gh <- do.call(pcfinhom, c(list(to_poly(P), lambda = lam, r = RG, correction = "Ripley", bw = BW_PCF), PCF_ARGS))
  col_of(gh, c("iso", "trans"))[match(round(R_SHOW, 6), round(gh$r, 6))]
}, error = function(e) rep(NA_real_, length(R_SHOW)))
npairs3 <- function(P) sum(closepaircounts(P, 3)) / 2
n_raw   <- function(S) vapply(as_list(S), npoints, numeric(1))
n_poly  <- function(S) vapply(as_list(S), function(P) npoints(to_poly(P)), numeric(1))
pn2_of  <- function(S) vapply(as_list(S), function(P) { Q <- to_poly(P); npairs3(Q) / npoints(Q)^2 }, numeric(1))
out_of  <- function(S) sum(vapply(as_list(S), function(P) npoints(P) - sum(inside.owin(P$x, P$y, W_km)), numeric(1)))

REPORT <- file.path(OUT_DIR, "31_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("31_generator_control_and_geyer_decomposition | ", format(Sys.time()), "\n", sep = "")
    cat("G3 multiplier", MULT_G3, "| production nrep", NREP_PROD, "\n\n")
    
    mu0 <- Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta)) + truth$b0
    mu_fill <- fill_im(mu0)
    lam4 <- exp(mu_fill + log(HI_MULT))
    LG <- list(L1 = c(var = 0.5, scale = 3), L4 = c(var = 1.5, scale = 8))
    
    ## ---------------- A. generator ----------------
    cat("== A. LGCP generator on aligned grids; pcf with TRUE intensity (x", HI_MULT, ") ==\n", sep = "")
    GEN <- data.frame()
    ## A0 Poisson control
    set.seed(SEED4 + 5)
    PO <- try_(sprintf("Poisson control, %d sims", N_PCF), lapply(seq_len(N_PCF), function(i) to_poly(rpoispp(lam4))))
    GP <- NULL
    if (!is.null(PO)) {
      cat(sprintf("  Poisson control mean n %.1f\n", mean(vapply(PO, npoints, numeric(1)))))
      GP <- t(vapply(PO, function(P) pcf_show(P, lam4), numeric(length(R_SHOW))))
      print(data.frame(r_km = R_SHOW, theory = 1, obs_mean = colMeans(GP, na.rm = TRUE),
                       se = apply(GP, 2, sd, na.rm = TRUE) / sqrt(N_PCF),
                       ratio = colMeans(GP, na.rm = TRUE)), digits = 3, row.names = FALSE)
      GEN <- rbind(GEN, data.frame(cell = "POISSON", r_km = R_SHOW, theory = 1, obs_mean = colMeans(GP, na.rm = TRUE),
                                   se = apply(GP, 2, sd, na.rm = TRUE) / sqrt(N_PCF)))
    }
    for (nm in names(LG)) {
      v <- unname(LG[[nm]]["var"]); s <- unname(LG[[nm]]["scale"])
      cat(sprintf("\n-- %s: var %.1f, scale %.0f km --\n", nm, v, s))
      P1 <- try_(sprintf("rLGCP %s x1, %d sims saveLambda", nm, N_FIELD), gen_lgcp(v, s, 0, N_FIELD, SEED4 + 10, TRUE))
      if (!is.null(P1)) {
        FS <- t(vapply(as_list(P1), field_stats, numeric(3), v = v))
        cat(sprintf("  realised u = log Lambda - mu inside W: mean %.3f (expect %.3f) | pixel var %.3f (expect about %.2f) | mean Lambda / mean exp(mu) %.3f (expect 1)\n",
                    mean(FS[, 1], na.rm = TRUE), -v / 2, mean(FS[, 2], na.rm = TRUE), v, mean(FS[, 3], na.rm = TRUE)))
      }
      PH <- try_(sprintf("rLGCP %s x%d, %d sims", nm, HI_MULT, N_PCF), gen_lgcp(v, s, log(HI_MULT), N_PCF, SEED4 + 20, FALSE))
      if (!is.null(PH)) {
        PHl <- as_list(PH)
        cat(sprintf("  x%d mean n %.1f (expect %d)\n", HI_MULT, mean(n_poly(PHl)), HI_MULT * RECOVERY_ENV$expected_n))
        G <- t(vapply(PHl, function(P) pcf_show(P, lam4), numeric(length(R_SHOW))))
        th <- exp(v * exp(-R_SHOW / s))
        print(data.frame(r_km = R_SHOW, theory = th, obs_mean = colMeans(G, na.rm = TRUE),
                         se = apply(G, 2, sd, na.rm = TRUE) / sqrt(N_PCF),
                         ratio = colMeans(G, na.rm = TRUE) / th,
                         log_excess = log(colMeans(G, na.rm = TRUE)) - log(th)), digits = 3, row.names = FALSE)
        GEN <- rbind(GEN, data.frame(cell = nm, r_km = R_SHOW, theory = th, obs_mean = colMeans(G, na.rm = TRUE),
                                     se = apply(G, 2, sd, na.rm = TRUE) / sqrt(N_PCF)))
      }
    }
    if (nrow(GEN)) write.csv(GEN, file.path(OUT_DIR, "31_lgcp_gen.csv"), row.names = FALSE)
    
    ## ---------------- B. Geyer decomposition ----------------
    cat("\n== B. Geyer conditional simulation: FITTED vs TRUE model (rbord 0, r 3, sat 3, p = 1, expand = 1) ==\n")
    P8 <- as_list(sim_geyer(MULT_G3, 2.5, 3, 3, 197, 8, SEED0 + 1))
    cat("8 G3 patterns identical to scripts 26-30:",
        identical(as.numeric(vapply(P8, function(P) npoints(to_poly(P)), numeric(1))), as.numeric(N26)), "\n")
    GR <- data.frame()
    for (k in seq_len(N_PAT_GY)) {
      Xk <- to_poly(P8[[k]]); nk <- npoints(Xk); pk <- npairs3(Xk) / nk^2
      fk <- quiet(geyer_fit(Xk, r = 3, sat = 3, rbord = 0))
      if (is.null(fk)) { cat(sprintf("  pattern %d: fit FAILED\n", k)); next }
      gh <- exp(unname(tail(coef(fk), 1)))
      for (src in c("fitted", "true")) {
        S <- tryCatch(
          if (src == "fitted") sim_fit_cond(fk, nk, N_SIM, SEED4 + 1000 * k)
          else sim_geyer_cond(MULT_G3, 2.5, 3, 3, nk, N_SIM, SEED4 + 1000 * k + 500),
          error = function(e) { cat("  simulate FAILED (", src, "): ", conditionMessage(e), "\n", sep = ""); NULL })
        if (is.null(S)) next
        Sl <- as_list(S); W1 <- Window(Sl[[1]])
        nr <- n_raw(Sl); np <- n_poly(Sl); pp <- pn2_of(Sl)
        GR <- rbind(GR, data.frame(k = k, source = src, n_obs = nk, gamma_fit = gh,
                                   sim_win_type = W1$type, sim_win_area_ratio = area.owin(W1) / area.owin(W_km),
                                   n_raw_min = min(nr), n_raw_max = max(nr), outside = out_of(Sl),
                                   n_ratio = mean(np) / nk, pn2_obs_e4 = 1e4 * pk, pn2_sim_e4 = 1e4 * mean(pp),
                                   pair_ratio = mean(pp) / pk))
        cat(sprintf("  pat %d %-6s: gamma_hat %.2f | sim window %s (area ratio %.3f) | raw n %d-%d, in-W n ratio %.3f, outside W %d | pairs/n^2 x1e4 obs %.2f sim %.2f ratio %.3f\n",
                    k, src, gh, W1$type, area.owin(W1) / area.owin(W_km), min(nr), max(nr), mean(np) / nk,
                    out_of(Sl), 1e4 * pk, 1e4 * mean(pp), mean(pp) / pk))
      }
    }
    if (nrow(GR)) {
      write.csv(GR, file.path(OUT_DIR, "31_geyer_decomp.csv"), row.names = FALSE)
      cat("\n-- Summary by source --\n")
      AG <- do.call(rbind, lapply(split(GR, GR$source), function(d)
        data.frame(source = d$source[1], n_patterns = nrow(d), pair_ratio = mean(d$pair_ratio),
                   sd_pair_ratio = sd(d$pair_ratio), n_ratio = mean(d$n_ratio), outside_total = sum(d$outside),
                   mean_win_area_ratio = mean(d$sim_win_area_ratio), mean_gamma_fit = mean(d$gamma_fit))))
      print(AG, digits = 3, row.names = FALSE)
    }
    
    ## ---------------- warnings ----------------
    cat("\n== Warning tally (whole run) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "31_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")