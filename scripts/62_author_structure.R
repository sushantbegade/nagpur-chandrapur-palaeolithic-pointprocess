# =====================================================================
# 62_author_structure.R
# Stage   : 13_sensitivity
# Folder  : 13_sensitivity/scripts/   (NEW)
# Purpose : EXPLORATORY (post hoc, provenance-based) analysis of the source structure of close pairs
#           (follow-up to scripts 60-61):
#           (A) share of site pairs with the same lead author / same source reference by distance band,
#               and a label-permutation test: lead-author labels shuffled within 30 km and 60 km blocks
#               (keeps regional footprints, removes local author clustering); statistic = same-author pairs
#               in the 0-2, 2-5, 5-10 km bands; one-sided Monte Carlo p, Bonferroni 0.05/6;
#           (B) leave-one-author-out frozen M1 test (E3, nsim 999) for authors with >= 15 sites,
#               plus pair counts vs M1 simulations for each remaining pattern.
# Outputs : outputs/62_report.txt, 62_author_share_bands.csv, 62_block_perm.csv, 62_loao_result.csv,
#           62_sessionInfo.txt, protocol_deviation_log_before_62.csv; docs/protocol_deviation_log.csv (1 row)
# Notes   : RUN ONCE. Lead author = first author of the source reference (not necessarily the survey team).
#           No survey-effort layer is built from author or reference counts. Reading rules (set before
#           running): (A) local author clustering indicated if p_up < 0.0083 in a band for a block size;
#           (B) per author: p < 0.04 deficiency retained, 0.04-0.06 precision rule (script 62b), > 0.06 not
#           retained (the deficiency then depends on that author's sites, alone or through sample size).
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
f_res <- file.path(OUT_DIR, "62_loao_result.csv")
if (file.exists(f_res)) stop("62_loao_result.csv exists: this analysis runs once. Stop; do not rerun.")
NPERM <- 9999; NCH <- 3; CH <- 333; PAIR_R <- c(1, 2, 3, 5); BLK <- c(30, 60); ALPHA_B <- 0.05 / 6; MIN_AUTH <- 15
SEED_A <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 160000

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
  item = "Exploratory (post hoc): source structure of close pairs - same-author shares by distance, block-permutation test, leave-one-author-out M1 (script 62)",
  reason = "Script 61 showed that close pairs share a lead author far more often than site pairs overall (all 43 GIS_MAP-GIS_MAP pairs within 5 km; about 67% of all pairs within 5 km) and that the GIS_MAP class is mostly five authors' legacy sites. The source structure of the pair excess is a survey/discovery-structure question the design lists among competing explanations; the analyses use the lead_author provenance variable and are post hoc. No survey-effort layer is built from author or reference counts.",
  effect_on_inference = "Reported as exploratory and used for interpretation (observation and survey structure), not for site selection; no frozen rule or primary result changes; Gate 6 FAIL stands.",
  approved_by = "author-delegated (assistant decision; exploratory, no frozen item changed)", stringsAsFactors = FALSE)
if (identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_62.csv"), overwrite = TRUE)
  if (!(ROW$item %in% LOG0$item)) write.csv(rbind(LOG0, ROW), f_log, row.names = FALSE)
  LOG_MSG <- "Log row appended (backup in outputs)."
} else { write.csv(ROW, file.path(OUT_DIR, "62_log_row_for_manual_append.csv"), row.names = FALSE); LOG_MSG <- "Log columns differ; row saved for manual append." }
pairs_at <- function(P, r) sum(closepaircounts(P, r)) / 2

REPORT <- file.path(OUT_DIR, "62_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  options(width = 200)
  cat("62_author_structure (EXPLORATORY) | ", format(Sys.time()), "\n", sep = "")
  cat(LOG_MSG, "\n")
  x <- SITES_G$x_m / 1000; y <- SITES_G$y_m / 1000; au <- SITES_G$lead_author; rf <- SITES_G$source_reference
  D <- as.matrix(dist(cbind(x, y))); ut <- upper.tri(D); d <- D[ut]; ia <- row(D)[ut]; ja <- col(D)[ut]
  sa <- au[ia] == au[ja]; sr <- rf[ia] == rf[ja]
  tab <- sort(table(au), decreasing = TRUE); big <- names(tab)[tab >= MIN_AUTH]
  cat("lead authors:", length(tab), "| with >=", MIN_AUTH, "sites:", paste(sprintf("%s %d", big, as.integer(tab[big])), collapse = ", "), "\n")
  
  cat("\n== A1. Share of site pairs with the same lead author / source reference, by distance band ==\n")
  BD <- list("0-2" = c(0, 2), "2-5" = c(2, 5), "5-10" = c(5, 10), "10-20" = c(10, 20), "20-50" = c(20, 50), "50-100" = c(50, 100), "100+" = c(100, Inf))
  SH <- do.call(rbind, lapply(names(BD), function(b) { m <- d >= BD[[b]][1] & d < BD[[b]][2]
  data.frame(band_km = b, pairs = sum(m), same_author = sum(sa[m]), share_same_author = mean(sa[m]), same_reference = sum(sr[m]), share_same_reference = mean(sr[m])) }))
  print(SH, digits = 3, row.names = FALSE); write.csv(SH, file.path(OUT_DIR, "62_author_share_bands.csv"), row.names = FALSE)
  
  cat("\n== A2. Same-author pairs vs lead-author labels shuffled within blocks (", NPERM, " permutations) ==\n", sep = "")
  close <- d < 10; ia1 <- ia[close]; ja1 <- ja[close]; d1 <- d[close]
  bm1 <- list("0-2" = d1 < 2, "2-5" = d1 >= 2 & d1 < 5, "5-10" = d1 >= 5)
  stat <- function(a) { s <- a[ia1] == a[ja1]; vapply(bm1, function(m) sum(s & m), numeric(1)) }
  o <- stat(au); BP <- data.frame()
  for (blk in BLK) {
    cell <- paste(floor((x - min(x)) / blk), floor((y - min(y)) / blk)); idx <- split(seq_along(au), cell); set.seed(SEED_A + blk)
    P <- vapply(seq_len(NPERM), function(k) { a <- au; for (ix in idx) if (length(ix) > 1) a[ix] <- au[ix][sample.int(length(ix))]; stat(a) }, numeric(3))
    BP <- rbind(BP, data.frame(block_km = blk, cells = length(idx), band_km = names(bm1), pairs_in_band = vapply(bm1, sum, 1), observed_same_author = o,
                               perm_mean = rowMeans(P), perm_lo = apply(P, 1, quantile, 0.025), perm_hi = apply(P, 1, quantile, 0.975),
                               p_up = vapply(seq_along(o), function(i) (1 + sum(P[i, ] >= o[i])) / (NPERM + 1), 1), row.names = NULL))
  }
  print(BP, digits = 3, row.names = FALSE); write.csv(BP, file.path(OUT_DIR, "62_block_perm.csv"), row.names = FALSE)
  hit <- BP[BP$p_up < ALPHA_B, ]
  cat("local author clustering indicated (p_up <", sprintf("%.4f", ALPHA_B), "):", if (nrow(hit)) paste(sprintf("%g km blocks, %s km", hit$block_km, hit$band_km), collapse = "; ") else "none", "\n")
  
  cat("\n== B. Leave-one-author-out frozen M1 test (nsim 999) ==\n")
  OUT <- data.frame()
  for (a in big) {
    keep <- au != a; S <- mk(paste0("LOAO_", a), X_obs_km[keep], W_km, covs, DUMMY_FIXED)
    cat(sprintf("  %s removed (%d sites): n %d | NA-covariate sites %d\n", a, sum(!keep), S$n, sum(cov_na(S$X, S))))
    ft <- rob_fit(S$X, S); ob <- rob_curves(S$X, ft, S); sL <- NULL; sG <- NULL; t0 <- unname(proc.time()["elapsed"])
    for (ch in seq_len(NCH)) { r <- chunk_run(S, CH, SEED_A + 1000 * match(a, big) + ch); sL <- cbind(sL, r$sL); sG <- cbind(sG, r$sG) }
    okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
    csL <- GET::create_curve_set(list(r = S$rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE])); csG <- GET::create_curve_set(list(r = S$rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
    p <- erl_p(csL, csG, 2); N <- sum(okc); se <- sqrt(p * (1 - p) / N)
    lam <- exp(Reduce(`+`, Map(function(u, v) u * v, S$covW, ft$beta[-1])) + ft$beta[1]); set.seed(SEED_A + 500 + match(a, big))
    psim <- t(vapply(seq_len(999), function(k) { X <- rob_poly(rpoispp(lam), S); vapply(PAIR_R, function(r) pairs_at(X, r), numeric(1)) }, numeric(length(PAIR_R))))
    pobs <- vapply(PAIR_R, function(r) pairs_at(S$X, r), numeric(1)); rat <- pobs / colMeans(psim)
    flag <- if (is.na(p)) "NA" else if (p >= 0.04 && p <= 0.06) "PRECISION" else if (p < 0.04) "DEFICIENCY retained" else "NOT retained"
    cat(sprintf("     usable %d | p %.4f (SE %.4f) | %s | pair ratios 2/3/5 km %.2f/%.2f/%.2f | %.0f s\n", N, p, se, flag, rat[2], rat[3], rat[4], unname(proc.time()["elapsed"]) - t0))
    OUT <- rbind(OUT, data.frame(author_removed = a, sites_removed = sum(!keep), n = S$n, nsim = N, p = p, mc_se = se, flag = flag, ratio_2 = rat[2], ratio_3 = rat[3], ratio_5 = rat[4],
                                 p_script42 = R42$p_erl_nstep2, seed_base = SEED_A, date = format(Sys.time())))
  }
  write.csv(OUT, f_res, row.names = FALSE)
  cat("\n"); print(OUT[, c("author_removed", "sites_removed", "n", "p", "mc_se", "flag", "ratio_2", "ratio_3", "ratio_5")], digits = 3, row.names = FALSE)
  cat("\nprimary (all 197 sites, script 42): p", R42$p_erl_nstep2, "| pair ratios 2/3/5 km about 2.9 / 3.1 / 2.6\n")
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "62_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")