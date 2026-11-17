# =============================================================================
# 22_probe_nstep_Ccells.R
# Stage: 10_simulation_recovery | closes the open questions left by script 21
# Purpose: (A) exact vs simplified window per-pattern cost, measured in the main session;
#   (B) diagnose why nstep = 1 gave no p-value; (C) Geyer truth multipliers for the
#   n = 75 and n = 114 cells (C1, C3); (D) LGCP truth mean n for n = 75 and 114 (C2, C4).
# SYNTHETIC patterns only. Requires R/recovery_setup.R (script 17) and
# R/recovery_functions.R as written by script 21.
# Outputs (10_simulation_recovery/outputs/): probe_report.txt, geyer_multipliers_C.csv,
#   lgcp_mean_n_C.csv, sessionInfo_22.txt
# Not executed by assistant. RESTART R FIRST, then run top to bottom (about 10 min).
# WARNING: first line clears the objects in your R workspace (not files).
# =============================================================================

rm(list = ls(envir = globalenv()), envir = globalenv()); invisible(gc())

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
PROBE_N <- 10L; NSIM_NSTEP <- 19L
MAX_IT <- 6L; TOL <- 0.04; RMH_NREP <- 5e5; RMH_BATCH <- 8L; LGCP_NSIM <- 20L

while (sink.number() > 0) sink()
source(REPO_CONFIG)
source(file.path(REPO_DIR, "R", "recovery_setup.R"))
func_f <- file.path(REPO_DIR, "R", "recovery_functions.R")
stopifnot(file.exists(func_f), any(grepl("time_probe", readLines(func_f))))
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
out <- dir_of("10_simulation_recovery"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)

bnd <- st_transform(st_make_valid(st_read(file.path(RAW_DIR, man$path[man$layer == "Boundary"]), quiet = TRUE)), CRS_PROJ)
win_sf <- st_union(st_geometry(bnd))
simp <- st_union(st_make_valid(st_simplify(win_sf, preserveTopology = TRUE, dTolerance = 100)))
W_simpl <- rescale(as.owin(simp), 1000, c("kilometre", "kilometres"))
W_exact <- W_km; X_exact <- X_obs_km
stopifnot(W_simpl$type == "polygonal", all(inside.owin(X_exact$x, X_exact$y, W_simpl)))

sink(file.path(out, "probe_report.txt"), split = TRUE)
tryCatch({
  
  sec("A window probe", {
    cat("==== A. PER-PATTERN COST: exact vs simplified window (main session, ", PROBE_N, " patterns each) ====\n", sep = "")
    source(func_f)
    pe <- time_probe(PROBE_N)
    W_km <- W_simpl; X_obs_km <- X_exact[W_simpl]
    source(func_f)
    ps <- time_probe(PROBE_N)
    W_km <- W_exact; X_obs_km <- X_exact
    source(func_f)
    f <- function(p) sprintf("vertices %d | area %.2f | simulate %.3f s + fit %.3f s + curves %.3f s = %.3f s", p["vertices"], p["area"], p["sim"], p["fit"], p["curves"], p["sim"] + p["fit"] + p["curves"])
    cat("exact      :", f(pe), "\nsimplified :", f(ps), "\n")
    cat(sprintf("per-pattern speed-up of the simplified window: x%.2f\n", sum(pe[c("sim", "fit", "curves")]) / sum(ps[c("sim", "fit", "curves")])))
  })
  
  sec("B nstep", {
    cat("\n==== B. WHY nstep = 1 RETURNED NO P-VALUE ====\n")
    set.seed(SEEDS$master + SEEDS$offsets["empirical"] + 77)
    Xo <- to_poly(rpoispp(truth$lamW)); bw <- bw_for_n(RECOVERY_ENV$expected_n)
    f0 <- fast_fit(Xo); ob <- fast_curves(Xo, f0, bw); lam_im <- lam_im_of(f0$beta)
    sL <- matrix(NA_real_, length(rL), NSIM_NSTEP); sG <- matrix(NA_real_, length(rG), NSIM_NSTEP)
    for (k in seq_len(NSIM_NSTEP)) { Xs <- to_poly(rpoispp(lam_im)); fs <- fast_fit(Xs); cs <- fast_curves(Xs, fs, bw); sL[, k] <- cs$L; sG[, k] <- cs$g }
    csL <- GET::create_curve_set(list(r = rL, obs = ob$L, sim_m = sL)); csG <- GET::create_curve_set(list(r = rG, obs = ob$g, sim_m = sG))
    for (ns in 1:2) {
      res <- tryCatch(GET::global_envelope_test(curve_sets = list(csL, csG), type = "erl", alpha = TESTS$alpha, nstep = ns), error = function(e) e)
      if (inherits(res, "error")) cat(sprintf("nstep = %d: ERROR: %s\n", ns, conditionMessage(res))) else cat(sprintf("nstep = %d: p = %s\n", ns, format(attr(res, "p"))))
    }
  })
  
  sec("C Geyer C cells", {
    cat("\n==== C. GEYER TRUTH MULTIPLIERS FOR n = 75 AND n = 114 (gamma 2.5, r 3 km, sat 3; 5e5 steps, batches of ", RMH_BATCH, ") ====\n", sep = "")
    mk <- function(g, r, m) rmhmodel(cif = "geyer", par = list(beta = 1, gamma = g, r = r, sat = 3), w = W_km, trend = truth$lamW * m)
    cells <- data.frame(cell = c("C1", "C3"), n_target = c(75, 114))
    cells$start <- 0.1659 * (cells$n_target / 197)^(1 / 1.5)
    res_g <- list()
    for (k in seq_len(nrow(cells))) {
      m <- cells$start[k]; hist <- NULL; ok <- FALSE
      for (it in seq_len(MAX_IT)) {
        set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 1000 + k * 10 + it)
        PP <- rmh(mk(2.5, 3, m), start = list(n.start = cells$n_target[k]), control = list(nrep = RMH_NREP, expand = 1), nsim = RMH_BATCH, verbose = FALSE)
        ns <- sapply(PP, function(z) npoints(z[W_km])); nb <- mean(ns)
        hist <- rbind(hist, data.frame(cell = cells$cell[k], target = cells$n_target[k], iter = it, mult = round(m, 4), mean_n = round(nb, 1), sd_n = round(sd(ns), 1)))
        if (abs(nb / cells$n_target[k] - 1) < TOL) { ok <- TRUE; break }
        m <- m * (cells$n_target[k] / nb)^(1 / 1.5)
      }
      print(hist, row.names = FALSE); res_g[[k]] <- cbind(hist[nrow(hist), ], converged = ok)
    }
    gm <- do.call(rbind, res_g); cat("Final multipliers:\n"); print(gm, row.names = FALSE)
    write.csv(gm, file.path(out, "geyer_multipliers_C.csv"), row.names = FALSE)
  })
  
  sec("D LGCP C cells", {
    cat("\n==== D. LGCP TRUTH MEAN n FOR n = 75 AND 114 (sigma2 1.5, scale 8 km) ====\n")
    rows <- list()
    for (tg in c(75, 114)) {
      tr <- make_truth(n_target = tg)
      set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 2000 + tg)
      PL <- rLGCP("exp", mu = log(tr$lamW) - 1.5 / 2, param = list(var = 1.5, scale = 8), saveLambda = FALSE, dimyx = c(504, 362), nsim = LGCP_NSIM)
      ns <- sapply(PL, function(z) npoints(z[W_km]))
      rows[[length(rows) + 1]] <- data.frame(target = tg, mean_n = round(mean(ns), 1), sd_n = round(sd(ns), 1), min_n = min(ns), max_n = max(ns))
    }
    lm_ <- do.call(rbind, rows); print(lm_, row.names = FALSE); write.csv(lm_, file.path(out, "lgcp_mean_n_C.csv"), row.names = FALSE)
  })
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_22.txt"))
cat("\nDONE. Outputs in:", out, "\n")