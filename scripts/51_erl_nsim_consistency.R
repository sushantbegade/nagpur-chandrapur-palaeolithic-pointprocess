# =====================================================================
# 51_erl_nsim_consistency.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW)
# Purpose : EXPLORATORY check of the scale of the ERL p-value across nsim, using the 999
#           saved M1 simulations (43_m1_sims.rds): random subsets of 99 / 199 / 499
#           simulations (without replacement) -> combined ERL p (nstep 2), L-only and g-only
#           ERL p; compared with the full-set p and with the script-50 baseline replicates
#           (nsim 199, nominal coordinates, fresh simulations).
# Outputs : outputs/51_report.txt, 51_nsim_consistency.csv, 51_sessionInfo.txt
# Notes   : Subsets come from ONE fixed simulation set, so they show how the p-value
#           function behaves with nsim on that set; they do not reproduce fresh-simulation
#           variance. Changes no reported result. Reading rules (set before running):
#           consistent if |median p(199 subsets) - p(999)| <= 0.01 and the share of 199-subsets
#           with p < 0.05 is within 0.05 of the script-50 baseline retention.
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
f_r42 <- file.path(OUT_DIR, "42_m1_result.csv"); stopifnot(file.exists(f_r42))
f_b50 <- file.path(OUT_DIR, "50_jitter_replicates.csv"); stopifnot(file.exists(f_b50))
SM <- readRDS(f_sims); R42 <- read.csv(f_r42, stringsAsFactors = FALSE); B50 <- read.csv(f_b50, stringsAsFactors = FALSE)
stopifnot(all(c("rL", "rG", "ob", "sL", "sG") %in% names(SM)), ncol(SM$sL) == 999, ncol(SM$sG) == 999)
SEED_C <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 80000
SPEC <- data.frame(N = c(99, 199, 499), B = c(300, 300, 200))

csets <- function(idx) list(
  L = GET::create_curve_set(list(r = SM$rL, obs = SM$ob$L, sim_m = SM$sL[, idx, drop = FALSE])),
  G = GET::create_curve_set(list(r = SM$rG, obs = SM$ob$g, sim_m = SM$sG[, idx, drop = FALSE])))
p_one <- function(cs) tryCatch(as.numeric(attr(GET::global_envelope_test(curve_sets = cs, type = "erl", alpha = 0.05), "p")), error = function(e) NA_real_)
p_sub <- function(idx) {
  cs <- csets(idx)
  c(comb = erl_p(cs$L, cs$G, 2), L = p_one(cs$L), g = p_one(cs$G))
}

REPORT <- file.path(OUT_DIR, "51_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("51_erl_nsim_consistency (EXPLORATORY) | ", format(Sys.time()), "\n", sep = "")
  full <- p_sub(seq_len(999))
  cat(sprintf("full set (999): combined p %.4f (script 42: %.4f) | L-only %.4f | g-only %.4f\n", full["comb"], R42$p_erl_nstep2, full["L"], full["g"]))
  stopifnot(isTRUE(abs(full["comb"] - R42$p_erl_nstep2) < 1e-9))
  set.seed(SEED_C)
  ROWS <- data.frame()
  for (i in seq_len(nrow(SPEC))) {
    N <- SPEC$N[i]; B <- SPEC$B[i]; t0 <- unname(proc.time()["elapsed"])
    P <- t(vapply(seq_len(B), function(b) p_sub(sample.int(999, N)), numeric(3)))
    ROWS <- rbind(ROWS, data.frame(N = N, b = seq_len(B), comb = P[, 1], L = P[, 2], g = P[, 3]))
    cat(sprintf("  N = %3d: %d subsets done in %.0f s\n", N, B, unname(proc.time()["elapsed"]) - t0))
  }
  write.csv(ROWS, file.path(OUT_DIR, "51_nsim_consistency.csv"), row.names = FALSE)
  q <- function(x) quantile(x, c(0.1, 0.5, 0.9), na.rm = TRUE)
  SUMM <- do.call(rbind, lapply(split(ROWS, ROWS$N), function(d) data.frame(
    N = d$N[1], subsets = nrow(d), comb_mean = mean(d$comb, na.rm = TRUE), comb_q10 = q(d$comb)[1], comb_median = q(d$comb)[2], comb_q90 = q(d$comb)[3],
    comb_share_lt05 = mean(d$comb < 0.05, na.rm = TRUE), L_median = q(d$L)[2], L_share_lt05 = mean(d$L < 0.05, na.rm = TRUE),
    g_median = q(d$g)[2], g_share_lt05 = mean(d$g < 0.05, na.rm = TRUE), n_na = sum(is.na(d$comb)))))
  cat("\n== Subset distributions (full-set values: combined 0.031, L 0.022, g 0.026) ==\n")
  print(SUMM, digits = 3, row.names = FALSE)
  bs <- B50[B50$type == "base" & !is.na(B50$p), ]; jt <- B50[B50$type == "jitter" & !is.na(B50$p), ]
  cat(sprintf("\nscript-50 baseline (fresh simulations, nsim 199, %d replicates): mean p %.4f, median %.4f, share p < 0.05 %.3f\n", nrow(bs), mean(bs$p), median(bs$p), mean(bs$p < 0.05)))
  cat(sprintf("script-50 jitter (nsim 199, %d replicates): mean p %.4f, median %.4f, share p < 0.05 %.3f\n", nrow(jt), mean(jt$p), median(jt$p), mean(jt$p < 0.05)))
  s199 <- SUMM[SUMM$N == 199, ]
  ok1 <- abs(s199$comb_median - full["comb"]) <= 0.01; ok2 <- abs(s199$comb_share_lt05 - mean(bs$p < 0.05)) <= 0.05
  cat(sprintf("\nconsistency rule: |median(199 subsets) - p(999)| = %.4f (<= 0.01: %s) | share < 0.05 at N = 199: %.3f vs baseline retention %.3f (within 0.05: %s)\n",
              abs(s199$comb_median - full["comb"]), ok1, s199$comb_share_lt05, mean(bs$p < 0.05), ok2))
  cat("verdict:", if (ok1 && ok2) "CONSISTENT: no nsim-dependent p scale detected on this set" else "NOT CONSISTENT: p scale differs with nsim and/or code path; investigate before comparing nsim-99/199 and nsim-999 results", "\n")
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "51_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")