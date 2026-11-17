# =====================================================================
# 24_geyer_multipliers_simplified.R   (v3: fixes rmhcontrol/rmhstart check)
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (OVERWRITES v2 of this file)
# Purpose : (1) re-check Geyer truth multipliers G1-G4 (n = 197) on the
#           100 m simplified window (recovery_cache_simpl.rds);
#           (2) calibrate C-cell Geyer multipliers C1 (n = 75) and
#           C3 (n = 114) with the corrected update exponent (2.83);
#           (3) independent verification run for every multiplier.
# Outputs : outputs/geyer_calibration_log_simpl.csv  (every iteration)
#           outputs/geyer_multipliers_simpl.csv      (final + verification)
#           outputs/24_report.txt, outputs/24_sessionInfo.txt
# Notes   : Does not touch frozen items. Geyer truth parametrisation
#           (beta = 1, trend = truth$lamW * mult) is INFERRED from the
#           multiplier scaling; iteration 0 of G1-G4 tests it.
#           Signature checks are informational; the smoke test is the gate.
#           Restart R before running. Main session only (no workers).
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))     # freeze params, cache, make_truth, to_poly
source(file.path(REPO_DIR, "R", "recovery_functions.R"))  # defines truth, make_quad, etc.

# ---- settings (calibration tooling only; not frozen design items) ----
EXPONENT   <- 2.83     # n ~ mult^2.83 (project_state)
NREP       <- 5e5      # rmh steps, as for G1-G4 calibration
NSIM_ITER  <- 40       # patterns per calibration iteration
NSIM_VER   <- 100      # patterns per independent verification
MAX_IT     <- 6        # max updates after iteration 0
TOL_REL    <- 0.05     # accept if |mean n / target - 1| <= 5%
SEED_BASE  <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 4000

OUT_DIR <- dir_of("10_simulation_recovery", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- open all inputs BEFORE sink ----
f_old <- file.path(OUT_DIR, "geyer_multipliers.csv")
G_OLD_FALLBACK <- c(G1 = 0.5815, G2 = 0.3323, G3 = 0.1659, G4 = 0.0714)  # project_state values
old_src <- "fallback (project_state values)"
g_old <- G_OLD_FALLBACK
if (file.exists(f_old)) {
  d_old <- read.csv(f_old, stringsAsFactors = FALSE)
  if (all(c("cell", "mult") %in% names(d_old)) && all(names(G_OLD_FALLBACK) %in% d_old$cell)) {
    g_old <- setNames(d_old$mult[match(names(G_OLD_FALLBACK), d_old$cell)], names(G_OLD_FALLBACK))
    old_src <- f_old
  }
}
C_START <- c(C1 = 0.1193, C3 = 0.1383)   # stated start values (unverified)

# ---- cell table: all Geyer cells from the frozen design ----
RC <- RECOVERY_CELLS[RECOVERY_CELLS$mechanism == "Geyer", c("cell", "n", "gamma", "r_km", "sat")]
RC$mult0 <- ifelse(RC$cell %in% names(g_old), g_old[RC$cell], C_START[RC$cell])
stopifnot(!anyNA(RC$mult0), nrow(RC) == 6)

# ---- pure helpers (no assignments inside) ----
sim_geyer <- function(mult, gamma, r, sat, n_start, nsim, seed, nrep = NREP) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
n_of <- function(pats) vapply(as_list(pats), function(P) npoints(to_poly(P)), numeric(1))
show_formals <- function(generic, cls) {
  mth <- utils::getS3method(generic, cls, optional = TRUE)
  if (is.null(mth)) { cat(sprintf("S3 method %s.%s: NOT FOUND\n", generic, cls)); return(invisible(NULL)) }
  cat(sprintf("%s.%s formals: %s\n", generic, cls, paste(names(formals(mth)), collapse = ", ")))
}

# ---- main, under sink ----
REPORT <- file.path(OUT_DIR, "24_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("24_geyer_multipliers_simplified v3 | ", format(Sys.time()), "\n", sep = "")
  cat("Old multipliers source:", old_src, "\n\n")
  
  ## 0. environment checks
  cat("== 0. Checks ==\n")
  nv <- sum(vapply(W_km$bdry, function(b) length(b$x), 1L))
  cat(sprintf("Window: type %s, vertices %d, area %.1f km2 (expect simplified: ~1,378 vertices, ~21,262 km2)\n",
              W_km$type, nv, area.owin(W_km)))
  stopifnot(W_km$type == "polygonal", nv < 5000)
  ins <- inside.owin(X_obs_km$x, X_obs_km$y, W_km)
  cat(sprintf("Observed sites inside simplified window: %d of %d (expect 197 of 197)\n", sum(ins), length(ins)))
  cat(sprintf("Truth lamW integral (should be %d): %.3f\n", RECOVERY_ENV$expected_n, integral(truth$lamW)))
  
  cat("-- signatures (informational; generics use '...', arguments live in default methods) --\n")
  show_formals("rmh", "default")
  show_formals("rmhcontrol", "default")
  show_formals("rmhstart", "default")
  show_formals("rmhmodel", "default")
  
  cat("Cross-check start values from G3 via exponent 1/2.83:\n")
  cat(sprintf("  n = 75 : %.4f (stated 0.1193)   n = 114: %.4f (stated 0.1383)\n",
              g_old["G3"] * (75 / 197)^(1 / EXPONENT), g_old["G3"] * (114 / 197)^(1 / EXPONENT)))
  
  cat("Smoke test: 1 pattern, G1 parameters, 1e4 steps ... ")
  tsm <- system.time(sm <- tryCatch(
    sim_geyer(RC$mult0[1], RC$gamma[1], RC$r_km[1], RC$sat[1], 197, 1, SEED_BASE, nrep = 1e4),
    error = function(e) { cat("\nSMOKE TEST FAILED: ", conditionMessage(e), "\n", sep = ""); stop("smoke test failed; paste output") }
  ))["elapsed"]
  cat(sprintf("ok, class %s, n = %d, %.1f s\n", class(sm)[1], n_of(sm), tsm))
  
  ## 1. calibration loop (top level, no replicate wrappers)
  LOG <- data.frame()
  FINAL <- data.frame()
  for (ci in seq_len(nrow(RC))) {
    cell <- RC$cell[ci]; target <- RC$n[ci]
    gam <- RC$gamma[ci]; rr <- RC$r_km[ci]; ss <- RC$sat[ci]
    mult <- RC$mult0[ci]; converged <- FALSE
    cat(sprintf("\n== %s: n target %d, gamma %.2f, r %.0f km, sat %d, start mult %.5f ==\n",
                cell, target, gam, rr, ss, mult))
    for (it in 0:MAX_IT) {
      t0 <- unname(proc.time()["elapsed"])
      pats <- sim_geyer(mult, gam, rr, ss, target, NSIM_ITER, SEED_BASE + ci * 100 + it)
      ns <- n_of(pats)
      m <- mean(ns); s <- sd(ns); se <- s / sqrt(length(ns)); rel <- m / target - 1
      ok <- abs(rel) <= TOL_REL
      cat(sprintf("  it %d  mult %.5f  mean n %.1f  sd %.1f  se %.1f  rel %+.3f  %s  (%.0f s)\n",
                  it, mult, m, s, se, rel, if (ok) "ACCEPT" else "update",
                  unname(proc.time()["elapsed"]) - t0))
      LOG <- rbind(LOG, data.frame(cell = cell, iter = it, mult = mult, mean_n = m, sd_n = s,
                                   se_n = se, rel = rel, accepted = ok))
      if (ok) { converged <- TRUE; break }
      fac <- min(max((target / max(m, 1))^(1 / EXPONENT), 0.7), 1.4)
      mult <- mult * fac
    }
    # independent verification, fresh seeds, no further update
    tv <- unname(proc.time()["elapsed"])
    pv <- sim_geyer(mult, gam, rr, ss, target, NSIM_VER, SEED_BASE + ci * 100 + 90)
    nv_ <- n_of(pv)
    cat(sprintf("  VERIFY mult %.5f  mean n %.1f  sd %.1f  se %.1f  rel %+.3f  min %d  max %d  (%.0f s)\n",
                mult, mean(nv_), sd(nv_), sd(nv_) / sqrt(length(nv_)), mean(nv_) / target - 1,
                min(nv_), max(nv_), unname(proc.time()["elapsed"]) - tv))
    FINAL <- rbind(FINAL, data.frame(cell = cell, target = target, gamma = gam, r_km = rr, sat = ss,
                                     mult = mult, start_mult = RC$mult0[ci], iter_converged = converged,
                                     verify_mean_n = mean(nv_), verify_sd_n = sd(nv_),
                                     verify_se_n = sd(nv_) / sqrt(length(nv_)),
                                     verify_rel = mean(nv_) / target - 1,
                                     verify_min = min(nv_), verify_max = max(nv_),
                                     nsim_verify = length(nv_), window = "simplified 100 m"))
    write.csv(LOG, file.path(OUT_DIR, "geyer_calibration_log_simpl.csv"), row.names = FALSE)   # checkpoint
    write.csv(FINAL, file.path(OUT_DIR, "geyer_multipliers_simpl.csv"), row.names = FALSE)
  }
  
  ## 2. summary
  cat("\n== SUMMARY ==\n")
  print(FINAL[, c("cell", "target", "start_mult", "mult", "iter_converged",
                  "verify_mean_n", "verify_sd_n", "verify_rel")], row.names = FALSE, digits = 5)
  cat("\nIteration-0 result for G1-G4 (old multiplier, simplified window) is the parametrisation check.\n")
  cat("Files: geyer_calibration_log_simpl.csv, geyer_multipliers_simpl.csv\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "24_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")