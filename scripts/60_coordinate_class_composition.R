# =====================================================================
# 60_coordinate_class_composition.R
# Stage   : 13_sensitivity
# Folder  : 13_sensitivity/scripts/   (NEW)
# Purpose : EXPLORATORY (post hoc, not frozen; provenance-based) check of whether the pair excess is
#           carried by the low-precision GIS_MAP coordinates:
#           (A) class composition of close pairs (bands 0-2, 2-5, 5-10 km; 0-5 km descriptive)
#               against 9999 random permutations of the coordinate-class labels over the fixed
#               site locations (number of pairs with exactly one / both GIS_MAP sites);
#           (B) frozen M1 test (E3, nsim 999, generic code, 3 chunks of 333) on the 137 sites with
#               FIELD or AUTHOR_REPORTED coordinates (GIS_MAP excluded), with pair counts vs simulations.
# Outputs : outputs/60_report.txt, 60_pair_composition.csv, 60_gis_excluded_result.csv,
#           60_sessionInfo.txt, protocol_deviation_log_before_60.csv;
#           docs/protocol_deviation_log.csv (1 row)
# Notes   : RUN ONCE. Reading rules (set before running): GIS_MAP over-represented in close pairs if
#           p_up < 0.0083 (Bonferroni over 3 bands x 2 statistics) in a non-overlapping band; the
#           permutation null treats labels as exchangeable (labels may be spatially structured), so a
#           hit is a flag, not proof of an artefact. Subset M1: p < 0.04 retained; 0.04-0.06 precision
#           rule (nsim 10000, script 60b); > 0.06 not retained (inconclusive: Stage 1 power is modest).
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
GENERIC_FILE <- file.path(REPO_DIR, "R", "robust_generic.R"); stopifnot(file.exists(GENERIC_FILE))

OUT_DIR <- dir_of("13_sensitivity", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
f_res <- file.path(OUT_DIR, "60_gis_excluded_result.csv")
if (file.exists(f_res)) stop("60_gis_excluded_result.csv exists: this analysis runs once. Stop; do not rerun.")
NPERM <- 9999; NCH <- 3; CH <- 333; PAIR_R <- c(1, 2, 3, 5); ALPHA_B <- 0.05 / 6
SEED_G <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 150000

# ---- open inputs BEFORE sink ----
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
SITES_G <- read.csv(f_sites, stringsAsFactors = FALSE)
stopifnot(nrow(SITES_G) == 197, max(abs(X_obs_km$x - SITES_G$x_m / 1000), abs(X_obs_km$y - SITES_G$y_m / 1000)) < 1e-3)
f_r42 <- file.path(REPO_DIR, "11_generative_gof", "outputs", "42_m1_result.csv"); stopifnot(file.exists(f_r42))
R42 <- read.csv(f_r42, stringsAsFactors = FALSE)
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")
source(GENERIC_FILE)

ROW <- data.frame(
  date = "2026-11-14",
  item = "Exploratory (post hoc, not frozen): coordinate-class composition of close pairs and M1 test excluding GIS_MAP sites (script 60)",
  reason = "Coordinate classes FIELD (64), AUTHOR_REPORTED (73) and GIS_MAP (60) differ in precision and the GIS-map datum is unknown. A reviewer may ask whether the pair excess is produced by low-precision coordinates (for example village-level or snapped positions). The check uses a provenance variable fixed before the analysis, not a result; it is nevertheless post hoc and labelled exploratory.",
  effect_on_inference = "Reported as exploratory. It does not alter the primary result or any frozen rule, and no interaction claim is affected (Gate 6 FAIL stands).",
  approved_by = "author-delegated (assistant decision; exploratory, no frozen item changed)", stringsAsFactors = FALSE)
if (identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_60.csv"), overwrite = TRUE)
  if (!(ROW$item %in% LOG0$item)) write.csv(rbind(LOG0, ROW), f_log, row.names = FALSE)
  LOG_MSG <- "Log row appended (backup in outputs)."
} else { write.csv(ROW, file.path(OUT_DIR, "60_log_row_for_manual_append.csv"), row.names = FALSE); LOG_MSG <- "Log columns differ; row saved for manual append." }
pairs_at <- function(P, r) sum(closepaircounts(P, r)) / 2

REPORT <- file.path(OUT_DIR, "60_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("60_coordinate_class_composition (EXPLORATORY) | ", format(Sys.time()), "\n", sep = "")
  cat(LOG_MSG, "\n")
  cls <- SITES_G$coordinate_source
  cat("classes:", paste(sprintf("%s %d", names(table(cls)), as.integer(table(cls))), collapse = ", "), "\n")
  
  cat("\n== A. Class composition of close pairs vs label permutations (", NPERM, ") ==\n", sep = "")
  cp <- closepairs(X_obs_km, 10, twice = FALSE)
  BANDS <- list("0-2 km" = c(0, 2), "2-5 km" = c(2, 5), "5-10 km" = c(5, 10), "0-5 km" = c(0, 5))
  bm <- lapply(BANDS, function(b) cp$d >= b[1] & cp$d < b[2])
  gis <- as.integer(cls == "GIS_MAP")
  st <- function(g) { a <- g[cp$i] + g[cp$j]; unlist(lapply(bm, function(m) c(sum(a[m] == 1), sum(a[m] == 2)))) }
  o <- st(gis); set.seed(SEED_G); P <- vapply(seq_len(NPERM), function(k) st(sample(gis)), numeric(length(o)))
  bn <- rep(names(BANDS), each = 2); sn <- rep(c("one GIS_MAP site", "both GIS_MAP sites"), times = length(BANDS))
  COMP <- data.frame(band = bn, statistic = sn, pairs_in_band = rep(vapply(bm, sum, 1), each = 2), observed = o, perm_mean = rowMeans(P),
                     perm_lo = apply(P, 1, quantile, 0.025), perm_hi = apply(P, 1, quantile, 0.975),
                     p_up = vapply(seq_along(o), function(i) (1 + sum(P[i, ] >= o[i])) / (NPERM + 1), 1),
                     p_low = vapply(seq_along(o), function(i) (1 + sum(P[i, ] <= o[i])) / (NPERM + 1), 1), stringsAsFactors = FALSE)
  COMP$bonferroni_set <- COMP$band != "0-5 km"
  print(COMP, digits = 3, row.names = FALSE); write.csv(COMP, file.path(OUT_DIR, "60_pair_composition.csv"), row.names = FALSE)
  cat(sprintf("\nexpected share of pairs with at least one GIS_MAP site under random labels: %.3f | with both: %.3f\n",
              1 - (1 - 60 / 197) * (1 - 60 / 196) , (60 / 197) * (59 / 196)))
  cc <- paste(pmin(cls[cp$i], cls[cp$j]), pmax(cls[cp$i], cls[cp$j]), sep = " + ")
  cat("\nclass-pair counts for pairs within 5 km:\n"); print(table(cc[cp$d < 5]))
  hit <- COMP[COMP$bonferroni_set & COMP$p_up < ALPHA_B, ]
  cat("\ncomposition verdict (non-overlapping bands, p_up <", sprintf("%.4f", ALPHA_B), "):", if (nrow(hit)) paste(sprintf("%s (%s)", hit$band, hit$statistic), collapse = "; ") else "no over-representation of GIS_MAP sites in close pairs", "\n")
  
  cat("\n== B. Frozen M1 test excluding GIS_MAP sites ==\n")
  keep <- cls != "GIS_MAP"; S <- mk("R19_noGIS", X_obs_km[keep], W_km, covs, DUMMY_FIXED)
  cat(sprintf("sites %d | L range %.1f-%.1f km, g range %.1f-%.1f km | NA-covariate sites %d\n", S$n, min(S$rL), max(S$rL), min(S$rG), max(S$rG), sum(cov_na(S$X, S))))
  ft <- rob_fit(S$X, S); ob <- rob_curves(S$X, ft, S)
  sL <- NULL; sG <- NULL
  for (ch in seq_len(NCH)) { t0 <- unname(proc.time()["elapsed"]); r <- chunk_run(S, CH, SEED_G + ch); sL <- cbind(sL, r$sL); sG <- cbind(sG, r$sG)
  cat(sprintf("  chunk %d done in %.0f s [%s]\n", ch, unname(proc.time()["elapsed"]) - t0, format(Sys.time(), "%H:%M"))) }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = S$rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE])); csG <- GET::create_curve_set(list(r = S$rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  p <- erl_p(csL, csG, 2); N <- sum(okc); se <- sqrt(p * (1 - p) / N)
  lam <- exp(Reduce(`+`, Map(function(a, b) a * b, S$covW, ft$beta[-1])) + ft$beta[1]); set.seed(SEED_G + 99)
  psim <- t(vapply(seq_len(999), function(k) { X <- rob_poly(rpoispp(lam), S); vapply(PAIR_R, function(r) pairs_at(X, r), numeric(1)) }, numeric(length(PAIR_R))))
  pobs <- vapply(PAIR_R, function(r) pairs_at(S$X, r), numeric(1))
  flag <- if (is.na(p)) "NA" else if (p >= 0.04 && p <= 0.06) "PRECISION" else if (p < 0.04) "DEFICIENCY retained (GIS_MAP excluded)" else "NOT retained (GIS_MAP excluded; inconclusive: Stage 1 power is modest at n of about 137)"
  cat(sprintf("\nusable simulations %d of %d\nERL combined p (nstep 2, nsim %d): %.4f (SE %.4f) | primary (all 197 sites, script 42): %.4f\nclassification: %s\n", N, NCH * CH, N, p, se, R42$p_erl_nstep2, flag))
  PR <- data.frame(r_km = PAIR_R, observed = pobs, sim_mean = colMeans(psim), sim_lo = apply(psim, 2, quantile, 0.025), sim_hi = apply(psim, 2, quantile, 0.975), ratio = pobs / colMeans(psim))
  cat("\npair counts, GIS_MAP-excluded pattern vs its M1 simulations:\n"); print(PR, digits = 3, row.names = FALSE)
  write.csv(data.frame(setting = "R19_noGIS", n = S$n, nsim = N, p = p, mc_se = se, flag = flag, ratio_2 = PR$ratio[2], ratio_3 = PR$ratio[3], ratio_5 = PR$ratio[4],
                       p_script42 = R42$p_erl_nstep2, seed_base = SEED_G, date = format(Sys.time())), f_res, row.names = FALSE)
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "60_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")