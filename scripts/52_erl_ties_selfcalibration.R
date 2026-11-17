# =====================================================================
# 52_erl_ties_selfcalibration.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW)
# Purpose : EXPLORATORY validity check of the combined ERL p-value (type erl, nstep 2) by
#           exchangeability, using the 999 saved M1 simulations (43_m1_sims.rds):
#           (A) GET version / arguments and the p-value interval of the real observed pattern;
#           (B) self-calibration: 500 simulated curves as pseudo-observed vs N others
#               (N = 99, 199, 499, 998): rejection rate P(p < 0.05) vs the ideal exact-MC level
#               (N+1 = 100: 0.040; 200: 0.045; 500: 0.048; 999: 0.049), plus tie-interval width;
#           (C) real observed vs random subsets (100 each): p and interval by N.
# Outputs : outputs/52_report.txt, 52_selfcal.csv, 52_sessionInfo.txt
# Notes   : Simulated curves are i.i.d. from the fitted null, hence exchangeable: a valid
#           Monte Carlo p gives P(p < 0.05) = (ceil(0.05 (N+1)) - 1) / (N + 1).
#           Reading rules (set before running): VALID if the observed rejection rate lies within
#           the binomial 95% interval around the ideal level at every N; ANTICONSERVATIVE if
#           above it at N = 99 or 199; CONSERVATIVE if below it. Tie intervals are reported if
#           GET exposes attr 'p_interval'. Changes no reported result.
#           Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))

OUT_DIR <- dir_of("11_generative_gof", "outputs")
f_sims <- file.path(OUT_DIR, "43_m1_sims.rds"); stopifnot(file.exists(f_sims))
SM <- readRDS(f_sims); stopifnot(ncol(SM$sL) == 999, ncol(SM$sG) == 999)
SEED_C <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 90000
NPS <- 500; NLEV <- c(99, 199, 499, 998); NREAL <- 100

cs_for <- function(oL, oG, idx) list(
  GET::create_curve_set(list(r = SM$rL, obs = oL, sim_m = SM$sL[, idx, drop = FALSE])),
  GET::create_curve_set(list(r = SM$rG, obs = oG, sim_m = SM$sG[, idx, drop = FALSE])))
get_p <- function(cs) {
  r <- tryCatch(GET::global_envelope_test(curve_sets = cs, type = "erl", alpha = 0.05, nstep = 2), error = function(e) NULL)
  if (is.null(r)) return(c(p = NA_real_, plo = NA_real_, pup = NA_real_))
  pi <- attr(r, "p_interval")
  c(p = as.numeric(attr(r, "p")), plo = if (is.null(pi)) NA_real_ else as.numeric(pi)[1], pup = if (is.null(pi)) NA_real_ else as.numeric(pi)[2])
}
ideal <- function(N) (ceiling(0.05 * (N + 1)) - 1) / (N + 1)

REPORT <- file.path(OUT_DIR, "52_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("52_erl_ties_selfcalibration (EXPLORATORY) | ", format(Sys.time()), "\n", sep = "")
  
  cat("\n== A. GET and the real observed pattern vs all 999 ==\n")
  cat("GET version:", as.character(packageVersion("GET")), "\n")
  cat("global_envelope_test arguments:", paste(names(formals(GET::global_envelope_test)), collapse = ", "), "\n")
  r0 <- tryCatch(GET::global_envelope_test(curve_sets = cs_for(SM$ob$L, SM$ob$g, seq_len(999)), type = "erl", alpha = 0.05, nstep = 2), error = function(e) NULL)
  cat("attributes of the result:", if (is.null(r0)) "test failed" else paste(names(attributes(r0)), collapse = ", "), "\n")
  a0 <- get_p(cs_for(SM$ob$L, SM$ob$g, seq_len(999)))
  cat(sprintf("real observed vs 999: p %.4f | interval %s\n", a0["p"], if (is.na(a0["plo"])) "not exposed" else sprintf("%.4f-%.4f", a0["plo"], a0["pup"])))
  
  cat("\n== B. Self-calibration (", NPS, " pseudo-observed simulated curves) ==\n", sep = "")
  set.seed(SEED_C)
  pj <- sample.int(999, NPS)
  SC <- data.frame()
  for (N in NLEV) {
    t0 <- unname(proc.time()["elapsed"])
    R <- t(vapply(pj, function(j) { idx <- sample(setdiff(seq_len(999), j), N); get_p(cs_for(SM$sL[, j], SM$sG[, j], idx)) }, numeric(3)))
    SC <- rbind(SC, data.frame(level = "pseudo", N = N, rep = seq_len(NPS), p = R[, 1], plo = R[, 2], pup = R[, 3]))
    k <- sum(R[, 1] < 0.05, na.rm = TRUE); n <- sum(!is.na(R[, 1])); id <- ideal(N); ci <- binom.test(k, n)$conf.int
    cat(sprintf("  N = %3d (%.0f s): P(p < 0.05) = %.3f (95%% %.3f-%.3f) | ideal exact-MC %.3f | mean p %.3f | share p <= 0.10 %.3f%s\n",
                N, unname(proc.time()["elapsed"]) - t0, k / n, ci[1], ci[2], id, mean(R[, 1], na.rm = TRUE), mean(R[, 1] <= 0.10, na.rm = TRUE),
                if (all(is.na(R[, 2]))) " | no tie interval" else sprintf(" | P(upper end < 0.05) %.3f | mean interval width %.4f", mean(R[, 3] < 0.05, na.rm = TRUE), mean(R[, 3] - R[, 2], na.rm = TRUE))))
  }
  
  cat("\n== C. Real observed vs random subsets (", NREAL, " per N) ==\n", sep = "")
  for (N in c(99, 199, 499)) {
    R <- t(vapply(seq_len(NREAL), function(b) get_p(cs_for(SM$ob$L, SM$ob$g, sample.int(999, N))), numeric(3)))
    SC <- rbind(SC, data.frame(level = "real", N = N, rep = seq_len(NREAL), p = R[, 1], plo = R[, 2], pup = R[, 3]))
    cat(sprintf("  N = %3d: mean p %.4f | %s\n", N, mean(R[, 1], na.rm = TRUE),
                if (all(is.na(R[, 2]))) "no tie interval" else sprintf("mean interval %.4f-%.4f (mean width %.4f) | share of upper ends < 0.05: %.3f", mean(R[, 2]), mean(R[, 3]), mean(R[, 3] - R[, 2]), mean(R[, 3] < 0.05))))
  }
  write.csv(SC, file.path(OUT_DIR, "52_selfcal.csv"), row.names = FALSE)
  
  cat("\n== Verdict (pre-set rule) ==\n")
  for (N in NLEV) {
    d <- SC[SC$level == "pseudo" & SC$N == N & !is.na(SC$p), ]; k <- sum(d$p < 0.05); ci <- binom.test(k, nrow(d))$conf.int; id <- ideal(N)
    cat(sprintf("  N = %3d: %s\n", N, if (id >= ci[1] && id <= ci[2]) "VALID (ideal level inside the 95% interval)" else if (k / nrow(d) > id) "ANTICONSERVATIVE" else "CONSERVATIVE"))
  }
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "52_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")