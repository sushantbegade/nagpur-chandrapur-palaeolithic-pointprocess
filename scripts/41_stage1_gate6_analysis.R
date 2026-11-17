# =====================================================================
# 41_stage1_gate6_analysis.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : single pre-set analysis of Stage 1 (19 cells x 100 reps, nsim 99):
#           class counts; Gate 6 metrics per cell (power, false-positive rate,
#           mechanism recovery) with Wilson 95% CIs and tri-state verdicts against
#           the frozen thresholds; Stage 2 trigger (any AMBIG metric); stage reach;
#           stratification by realised n; estimator distributions; sensitivity of
#           the Geyer class to the z-max definition; timing.
# Outputs : outputs/41_report.txt, 41_class_counts.csv, 41_gate6_cells.csv,
#           41_stage2_cells.csv, 41_sessionInfo.txt
# Notes   : Reads 40_stage1_all.csv and R/freeze_v1_params.R only. No simulation.
#           Gate 6 itself needs cells bracketing the EMPIRICAL estimates, which do not
#           exist yet; this script gives per-cell and pooled verdicts, not the gate.
#           Verdict rules (set before running): PASS if the Wilson 95% interval lies
#           wholly on the passing side of the threshold, FAIL if wholly on the failing
#           side, AMBIG if it contains the threshold. Descriptive metrics get no verdict.
#           Restart R before running.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

OUT_DIR <- dir_of("10_simulation_recovery", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ---- open all inputs BEFORE sink ----
f_all <- file.path(OUT_DIR, "40_stage1_all.csv"); stopifnot(file.exists(f_all))
f_par <- file.path(REPO_DIR, "R", "freeze_v1_params.R"); stopifnot(file.exists(f_par))
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
source(f_par)
if (!grepl("^FROZEN", FREEZE_STATUS)) stop("Freeze v1 is not FROZEN. Aborting.")
A   <- read.csv(f_all, stringsAsFactors = FALSE)
LOG <- read.csv(f_log, stringsAsFactors = FALSE)
RC  <- RECOVERY_CELLS
A_THR <- DECISION$power_min; F_THR <- DECISION$fpr_max; M_THR <- DECISION$mechanism_recovery_min
ALPHA <- TESTS$alpha

# ---- helpers ----
wilson <- function(k, n) if (n == 0) c(NA_real_, NA_real_) else as.numeric(suppressWarnings(prop.test(k, n, correct = FALSE)$conf.int))
tri_ge <- function(lo, hi, thr) if (is.na(lo)) NA_character_ else if (lo >= thr) "PASS" else if (hi < thr) "FAIL" else "AMBIG"
tri_le <- function(lo, hi, thr) if (is.na(lo)) NA_character_ else if (hi <= thr) "PASS" else if (lo > thr) "FAIL" else "AMBIG"
metric_row <- function(cell, truth, metric, k, n, thr, rule) {
  ci <- wilson(k, n)
  v <- if (rule == "ge") tri_ge(ci[1], ci[2], thr) else if (rule == "le") tri_le(ci[1], ci[2], thr) else "descriptive"
  data.frame(cell = cell, truth = truth, metric = metric, k = k, n = n, est = if (n) k / n else NA_real_,
             lo = ci[1], hi = ci[2], threshold = thr, rule = rule, verdict = v, stringsAsFactors = FALSE)
}
qf <- function(x, p) if (sum(!is.na(x))) as.numeric(quantile(x, p, na.rm = TRUE)) else NA_real_

REPORT <- file.path(OUT_DIR, "41_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("41_stage1_gate6_analysis | ", format(Sys.time()), "\n", sep = "")
  cat(sprintf("thresholds: power >= %.2f | false-positive <= %.2f | mechanism recovery >= %.2f | alpha %.2f\n", A_THR, F_THR, M_THR, ALPHA))
  
  ## ---------------- 0. integrity ----------------
  cat("\n== 0. Integrity ==\n")
  cat(sprintf("rows %d | cells %d | replicates per cell: %s\n", nrow(A), length(unique(A$cell)),
              paste(unique(table(A$cell)), collapse = ",")))
  cat("duplicate (cell, i):", sum(duplicated(A[, c("cell", "i")])), "| failed replicates:", sum(grepl("^failed", A$class)),
      "| error texts:", sum(!is.na(A$err) & A$err != ""), "\n")
  cat("classes present:", paste(names(table(A$class)), collapse = ", "), "\n")
  cat(sprintf("deviation log: %d rows; Stage 1 / sd-floor rows present: %s / %s\n", nrow(LOG),
              any(grepl("Stage 1 recovery run", LOG$item)), any(grepl("sd floor", LOG$item))))
  
  ## ---------------- 1. class counts ----------------
  cat("\n== 1. Class counts per cell (valid replicates) ==\n")
  CC <- do.call(rbind, lapply(RC$cell, function(cc) {
    d <- A[A$cell == cc & !grepl("^failed", A$class), ]
    data.frame(cell = cc, truth = RC$mechanism[RC$cell == cc], n_valid = nrow(d),
               Poisson = sum(d$class == "Poisson"), LGCP = sum(d$class == "LGCP"),
               Geyer = sum(d$class == "Geyer"), inadequate = sum(d$class == "inadequate"))
  }))
  print(CC, row.names = FALSE)
  write.csv(CC, file.path(OUT_DIR, "41_class_counts.csv"), row.names = FALSE)
  
  ## ---------------- 2. Gate 6 metrics ----------------
  cat("\n== 2. Gate 6 metrics per cell (Wilson 95% CI; verdict vs frozen threshold) ==\n")
  G6 <- data.frame()
  for (cc in RC$cell) {
    d <- A[A$cell == cc & !grepl("^failed", A$class), ]; n <- nrow(d)
    tr <- RC$mechanism[RC$cell == cc]
    kg <- sum(d$class == "Geyer")
    if (tr == "Geyer") G6 <- rbind(G6, metric_row(cc, tr, "power = P(Geyer)", kg, n, A_THR, "ge"))
    else G6 <- rbind(G6, metric_row(cc, tr, "false positive = P(Geyer)", kg, n, F_THR, "le"))
    if (tr == "LGCP") G6 <- rbind(G6, metric_row(cc, tr, "recovery = P(LGCP)", sum(d$class == "LGCP"), n, M_THR, "ge"))
    if (tr == "Poisson") G6 <- rbind(G6, metric_row(cc, tr, "P(Poisson) (descriptive)", sum(d$class == "Poisson"), n, NA_real_, "none"))
  }
  print(G6, digits = 3, row.names = FALSE)
  write.csv(G6, file.path(OUT_DIR, "41_gate6_cells.csv"), row.names = FALSE)
  cat("\nverdict counts:\n"); print(table(G6$verdict))
  
  cat("\n-- Pooled false-positive rate P(Geyer) over non-Gibbs cells --\n")
  for (grp in list(list("Poisson-truth cells", RC$cell[RC$mechanism == "Poisson"]),
                   list("LGCP-truth cells", RC$cell[RC$mechanism == "LGCP"]),
                   list("all non-Gibbs cells", RC$cell[RC$mechanism != "Geyer"]))) {
    d <- A[A$cell %in% grp[[2]] & !grepl("^failed", A$class), ]
    ci <- wilson(sum(d$class == "Geyer"), nrow(d))
    cat(sprintf("  %-22s %d/%d = %.3f (Wilson %.3f-%.3f)\n", grp[[1]], sum(d$class == "Geyer"), nrow(d),
                sum(d$class == "Geyer") / nrow(d), ci[1], ci[2]))
  }
  
  ## ---------------- 3. Stage 2 trigger ----------------
  cat("\n== 3. Stage 2 trigger (any AMBIG metric) ==\n")
  amb <- unique(G6$cell[!is.na(G6$verdict) & G6$verdict == "AMBIG"])
  S2 <- data.frame(cell = amb, stringsAsFactors = FALSE)
  cat("cells to extend to reps 101-200 (nsim 99):", if (length(amb)) paste(amb, collapse = " ") else "none", "\n")
  cat("cells with FAIL metrics (not extended; reported as failed):",
      paste(unique(G6$cell[!is.na(G6$verdict) & G6$verdict == "FAIL"]), collapse = " "), "\n")
  write.csv(S2, file.path(OUT_DIR, "41_stage2_cells.csv"), row.names = FALSE)
  
  ## ---------------- 4. stage reach ----------------
  cat("\n== 4. Stage reach and conditional test behaviour ==\n")
  SR <- do.call(rbind, lapply(RC$cell, function(cc) {
    d <- A[A$cell == cc & !grepl("^failed", A$class), ]
    r3 <- d[!is.na(d$p3), ]; r4 <- d[!is.na(d$p4), ]
    data.frame(cell = cc, M1_reject = mean(d$p1 < ALPHA), reach_M3 = nrow(r3),
               M3_reject_given_reach = if (nrow(r3)) mean(r3$p3 < ALPHA) else NA_real_,
               reach_M4 = nrow(r4), M4_ok_given_reach = if (nrow(r4)) mean(r4$p4 >= ALPHA) else NA_real_,
               Geyer_given_reach = if (nrow(r4)) mean(r4$class == "Geyer") else NA_real_)
  }))
  print(SR, digits = 3, row.names = FALSE)
  
  ## ---------------- 5. stratification by realised n (Geyer cells) ----------------
  cat("\n== 5. P(Geyer) by realised-n tercile within Geyer cells ==\n")
  for (cc in RC$cell[RC$mechanism == "Geyer"]) {
    d <- A[A$cell == cc & !grepl("^failed", A$class), ]
    br <- unique(quantile(d$n, c(0, 1/3, 2/3, 1)))
    if (length(br) < 3) { cat(cc, ": too few distinct n values\n"); next }
    grp <- cut(d$n, br, include.lowest = TRUE)
    TT <- do.call(rbind, lapply(levels(grp), function(l) { x <- d[grp == l, ]
    data.frame(cell = cc, n_range = l, reps = nrow(x), mean_n = mean(x$n), P_Geyer = mean(x$class == "Geyer"),
               P_LGCP = mean(x$class == "LGCP")) }))
    print(TT, digits = 3, row.names = FALSE)
  }
  
  ## ---------------- 6. estimators ----------------
  cat("\n== 6a. M3 estimates in LGCP-truth cells (replicates that reached M3) ==\n")
  EL <- do.call(rbind, lapply(RC$cell[RC$mechanism == "LGCP"], function(cc) {
    d <- A[A$cell == cc & !is.na(A$sigma2), ]
    data.frame(cell = cc, reached = nrow(d), s2_med = qf(d$sigma2, .5), s2_q25 = qf(d$sigma2, .25), s2_q75 = qf(d$sigma2, .75),
               al_med = qf(d$alpha, .5), al_q25 = qf(d$alpha, .25), al_q75 = qf(d$alpha, .75),
               share_alpha_gt10 = if (nrow(d)) mean(d$alpha > 10) else NA_real_, share_s2_lt001 = if (nrow(d)) mean(d$sigma2 < 0.01) else NA_real_)
  }))
  print(EL, digits = 3, row.names = FALSE)
  cat("\n== 6b. M4 estimates in Geyer-truth cells (replicates that reached M4) ==\n")
  for (cc in RC$cell[RC$mechanism == "Geyer"]) {
    d <- A[A$cell == cc & !is.na(A$gamma), ]
    cat(sprintf("%s (truth gamma %.1f, r %.0f, sat %.0f): reached %d | gamma_hat median %.3f (IQR %.3f-%.3f)\n", cc,
                RC$gamma[RC$cell == cc], RC$r_km[RC$cell == cc], RC$sat[RC$cell == cc], nrow(d),
                qf(d$gamma, .5), qf(d$gamma, .25), qf(d$gamma, .75)))
    if (nrow(d)) { cat("   selected r (km):", paste(names(table(d$r_sel)), table(d$r_sel), sep = "=", collapse = "  "),
                       "| selected sat:", paste(names(table(d$sat_sel)), table(d$sat_sel), sep = "=", collapse = "  "), "\n") }
  }
  
  ## ---------------- 7. sensitivity of the Geyer class ----------------
  cat("\n== 7. P(Geyer) under alternative M4 criteria (all valid replicates) ==\n")
  SE <- do.call(rbind, lapply(RC$cell, function(cc) {
    d <- A[A$cell == cc & !grepl("^failed", A$class), ]
    reach <- !is.na(d$p4)
    prim <- d$class == "Geyer"
    unfl <- reach & d$p4 >= ALPHA & is.finite(d$zraw4) & is.finite(d$zraw3) & d$zraw4 <= 0.5 * d$zraw3
    ponly <- reach & d$p4 >= ALPHA
    data.frame(cell = cc, truth = RC$mechanism[RC$cell == cc], primary_floored = mean(prim), unfloored = mean(unfl),
               p_only_no_z = mean(ponly), reached_M4 = sum(reach))
  }))
  print(SE, digits = 3, row.names = FALSE)
  
  ## ---------------- 8. timing ----------------
  cat("\n== 8. Mean seconds per stage (workers, nsim 99) ==\n")
  cat(sprintf("M1 %.1f (n %d) | M3 %.1f (n %d) | M4 incl. profile %.1f (n %d) | replicate total %.1f\n",
              mean(A$sec1, na.rm = TRUE), sum(!is.na(A$sec1)), mean(A$sec3, na.rm = TRUE), sum(!is.na(A$sec3)),
              mean(A$sec4, na.rm = TRUE), sum(!is.na(A$sec4)), mean(A$sec)))
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "41_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")