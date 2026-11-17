# =====================================================================
# 64_tables_and_register.R   (v2: T7b label fix; anisotropy split into T7d; T4a note once)
# Stage   : 15_tables
# Folder  : 15_tables/scripts/   (OVERWRITES v1)
# Purpose : assemble Tables 1-8 (Research_Design list) as CSV from EXISTING outputs only, and check
#           key values against those reported in the analysis conversation.
#           T1 provenance and database audit | T2 covariates | T3 models and mechanisms |
#           T4 parameter estimates | T5 generative goodness-of-fit | T6 Stage 1 recovery |
#           T7 robustness (M1/M3 settings, jitter, block validation, anisotropy) |
#           T8 sensitivity to observation, chronology and provenance assumptions.
# Outputs : outputs/64_report.txt, T1a ... T8 CSVs, 64_consistency_checks.csv, 64_sessionInfo.txt
# Notes   : No analysis is run. Missing inputs are listed, not guessed. v1 mislabelled the rows of
#           T7b (baseline and primary jitter swapped); v2 builds T7b by type and checks the table itself.
#           Figure 7 (Geyer diagnostics) is replaced by a Stage 1 stage-reach figure (M4 not fitted).
#           Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)
source(file.path(REPO_DIR, "R", "freeze_v1_params.R"))

D10 <- dir_of("10_simulation_recovery", "outputs"); D11 <- dir_of("11_generative_gof", "outputs")
D12 <- dir_of("12_validation", "outputs"); D13 <- dir_of("13_sensitivity", "outputs")
OUT <- dir_of("15_tables", "outputs"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

ENVM <- new.env(); ENVM$missing <- character(0)
rd <- function(dir, f) { p <- file.path(dir, f); if (!file.exists(p)) { ENVM$missing <- c(ENVM$missing, p); return(NULL) }; read.csv(p, stringsAsFactors = FALSE) }
g1 <- function(x) if (length(x) >= 1 && !is.na(x[1])) as.numeric(x[1]) else NA_real_
wilson <- function(k, n) as.numeric(suppressWarnings(prop.test(k, n, correct = FALSE)$conf.int))
mcse <- function(p, n) sqrt(p * (1 - p) / n)

# ---- open all inputs BEFORE sink ----
SITES <- read.csv(file.path(REPO_DIR, "data_derived", "sites_NA.csv"), stringsAsFactors = FALSE)
LOGS  <- read.csv(file.path(REPO_DIR, "docs", "site_inclusion_exclusion_log.csv"), stringsAsFactors = FALSE)
COVD  <- read.csv(file.path(REPO_DIR, "docs", "covariate_definitions.csv"), stringsAsFactors = FALSE)
R42 <- rd(D11, "42_m1_result.csv"); C42 <- rd(D11, "42_m1_coef.csv"); R44 <- rd(D11, "44_m2_result.csv"); R45 <- rd(D11, "45_m2_o2_precision_result.csv")
R46 <- rd(D11, "46_m3_result.csv"); B46 <- rd(D11, "46_m3_boot.csv"); R49 <- rd(D11, "49_robustness_summary.csv"); R49b <- rd(D11, "49b_precision_result.csv")
J50 <- rd(D11, "50_jitter_summary.csv"); R53 <- rd(D11, "53_m3_robustness_summary.csv"); R55 <- rd(D11, "55_exact_window_m1_result.csv")
A1 <- rd(D10, "40_stage1_all.csv"); V56 <- rd(D12, "56_blockcv.csv"); V57 <- rd(D12, "57_blockcv_exploratory.csv")
AN58 <- rd(D13, "58_anisotropy_bands.csv"); J59 <- rd(D13, "59_jitter_sens_summary.csv"); G60 <- rd(D13, "60_gis_excluded_result.csv")
L62 <- rd(D13, "62_loao_result.csv"); C63 <- rd(D13, "63_matched_control_summary.csv")

wr <- function(df, f) write.csv(df, file.path(OUT, f), row.names = FALSE)
safe <- function(label, expr) tryCatch(expr, error = function(e) cat(sprintf("[%s] FAILED: %s\n", label, conditionMessage(e))))
show <- function(df, n = 12) { print(head(df, n), digits = 4, row.names = FALSE); cat(sprintf("(%d rows)\n", nrow(df))) }

REPORT <- file.path(OUT, "64_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  options(width = 220)
  cat("64_tables_and_register v2 | ", format(Sys.time()), "\n", sep = "")
  
  cat("\n== T1a Provenance and database audit ==\n")
  safe("T1a", {
    act <- table(LOGS$action); cnt <- function(a) if (a %in% names(act)) as.integer(act[[a]]) else NA_integer_
    T1a <- data.frame(stage = c("Master inventory", "Identity / duplicate audit", "Identity / duplicate audit", "Provenance audit", "Provenance audit", "Provenance audit", "Provenance audit",
                                "Coordinate audit", "Coordinate audit", "Coordinate audit", "Inclusion criteria", "Spatial-window screening", "Analytical sample N_A"),
                      item = c("records in sites_master_decided.csv", "candidate pairs reviewed and confirmed distinct (no merge)", "records excluded as duplicates", "records with author-corrected fields, included",
                               "records flagged redeposited (included; sensitivity set)", "decisions stated by the author (AUTHOR_STATED)", "decisions by default inclusion (DEFAULT_UNQUEUED, not individually confirmed)",
                               "FIELD coordinates", "AUTHOR_REPORTED coordinates", "GIS_MAP coordinates (datum unknown)", "records excluded", "records inside the analytical window", "sites in primary dataset A"),
                      n = c(nrow(SITES), cnt("CONFIRMED_DISTINCT_NO_MERGE") / 2, 0, cnt("FIELD_CORRECTED_INCLUDE"), sum(SITES$sens_redeposited), sum(SITES$decision_basis == "AUTHOR_STATED"), sum(SITES$decision_basis == "DEFAULT_UNQUEUED"),
                            sum(SITES$coordinate_source == "FIELD"), sum(SITES$coordinate_source == "AUTHOR_REPORTED"), sum(SITES$coordinate_source == "GIS_MAP"), sum(SITES$inclusion_status != "INCLUDE"), sum(SITES$in_window), sum(SITES$dataset_A)),
                      stringsAsFactors = FALSE)
    wr(T1a, "T1a_audit.csv"); show(T1a, 20) })
  cat("\n== T1b Composition of the analytical sample ==\n")
  safe("T1b", {
    one <- function(v, lab) { t <- table(SITES[[v]]); data.frame(variable = lab, level = names(t), n = as.integer(t), stringsAsFactors = FALSE) }
    la <- sort(table(SITES$lead_author), decreasing = TRUE); top <- names(la)[1:8]
    T1b <- rbind(one("coordinate_source", "coordinate source"), one("site_type", "site type"), one("data_source", "data source"), one("phase_combo", "phase combination (pure and mixed)"),
                 data.frame(variable = "any-phase chronology", level = c("Lower Palaeolithic", "Middle Palaeolithic", "Upper Palaeolithic"), n = c(sum(SITES$LP), sum(SITES$MP), sum(SITES$UP))),
                 data.frame(variable = "analysis datasets", level = c("A (all)", "C (non-fluvial)", "redeposited (sensitivity)"), n = c(sum(SITES$dataset_A), sum(SITES$dataset_C), sum(SITES$sens_redeposited))),
                 data.frame(variable = "lead author (8 largest)", level = top, n = as.integer(la[top])),
                 data.frame(variable = "lead author (others)", level = paste(length(la) - 8, "authors"), n = sum(la) - sum(la[top])))
    wr(T1b, "T1b_composition.csv"); show(T1b, 40) })
  
  cat("\n== T2 Covariates ==\n")
  safe("T2", { T2 <- COVD[, intersect(c("covariate", "domain", "derivation", "transform", "df", "E1", "E2", "E3", "E4", "caveat"), names(COVD))]; wr(T2, "T2_covariates.csv"); show(T2) })
  cat("\n== T3 Models and mechanisms ==\n")
  safe("T3", {
    T3 <- data.frame(model = c("M0", "M1", "M2", "M3", "M4", "Strauss / Thomas"),
                     hypothesis = c("complete spatial randomness (exploratory reference)", "environmental inhomogeneous Poisson (E3, no observation offset)", "environment plus fixed observation offset (O2, O3) or free proxy coefficient (z_S2, z_S3)",
                                    "environment plus latent log-Gaussian field (exponential covariance)", "environment plus Geyer saturation interaction", "secondary / optional alternatives"),
                     estimator = c("-", "Poisson maximum pseudolikelihood (Berman-Turner)", "as M1 with fixed offset or free coefficient", "kppm composite likelihood (clik2); conditional simulation (rLGCP eps 0.5 + rpoint)", "maximum pseudolikelihood with profile search over r x saturation", "-"),
                     status_in_this_study = c("not run", "fitted (primary, O1); robustness in 17 settings", "fitted (O2, O3 primary robustness; z_S2, z_S3 in R15-R16)", "fitted (primary and 4 robustness settings)",
                                              "not fitted: hierarchy ended at M3 and Stage 1 showed Gibbs parameters not identifiable (Gate 6 FAIL)", "not run"), stringsAsFactors = FALSE)
    wr(T3, "T3_models.csv"); show(T3) })
  
  cat("\n== T4a E3 coefficients (Poisson standard errors assume independence; M1 is deficient, so they are anti-conservative; dual-role variables are not behavioural preferences) ==\n")
  safe("T4a", { T4a <- C42; T4a$rate_ratio <- exp(T4a$coef); wr(T4a, "T4a_E3_coefficients.csv"); show(T4a, 14) })
  cat("\n== T4b M3 parameters (estimator outputs; bootstrap at the point estimate) ==\n")
  safe("T4b", {
    q <- function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)
    T4b <- data.frame(parameter = c("sigma2 (latent variance)", "alpha (exponential scale, km)"), estimate = c(R46$sigma2_hat, R46$alpha_hat),
                      boot_lo = c(q(B46$sigma2)[1], q(B46$alpha)[1]), boot_median = c(q(B46$sigma2)[2], q(B46$alpha)[2]), boot_hi = c(q(B46$sigma2)[3], q(B46$alpha)[3]),
                      boot_share_at_bound = c(mean(B46$sigma2 < 0.01, na.rm = TRUE), mean(B46$alpha > R_MAX_DECISION_KM, na.rm = TRUE)),
                      note = c("estimator output; inflated in Stage 1 recovery", "estimator output; not a distance; alpha above the 10 km decision range in a share of refits"), stringsAsFactors = FALSE)
    wr(T4b, "T4b_M3_parameters.csv"); show(T4b) })
  cat("\n== T4c Stage 1 LGCP estimator reference (cells where M3 was reached) ==\n")
  safe("T4c", {
    RC <- RECOVERY_CELLS; lc <- RC$cell[RC$mechanism == "LGCP"]
    T4c <- do.call(rbind, lapply(lc, function(cc) { d <- A1[A1$cell == cc & !is.na(A1$sigma2), ]
    data.frame(cell = cc, nominal_sigma2 = RC$sigma2[RC$cell == cc], nominal_scale_km = RC$scale_km[RC$cell == cc], reached_M3 = nrow(d),
               sigma2_hat_median = median(d$sigma2), sigma2_hat_q25 = quantile(d$sigma2, .25), sigma2_hat_q75 = quantile(d$sigma2, .75),
               alpha_hat_median = median(d$alpha), alpha_hat_q25 = quantile(d$alpha, .25), alpha_hat_q75 = quantile(d$alpha, .75), share_alpha_gt_10 = mean(d$alpha > 10)) }))
    wr(T4c, "T4c_stage1_lgcp_estimators.csv"); show(T4c) })
  
  cat("\n== T5 Generative goodness-of-fit ==\n")
  safe("T5", {
    rows <- list()
    add <- function(model, setting, nsim, p, zr, zrr, note) rows[[length(rows) + 1]] <<- data.frame(model = model, setting = setting, nsim = nsim, p = p, mc_se = mcse(p, nsim), deficiency_at_0.05 = p < 0.05, raw_max_z = zr, at_r_km = zrr, note = note, stringsAsFactors = FALSE)
    add("M1 (O1)", "simplified window", R42$nsim, R42$p_erl_nstep2, R42$zraw_max, R42$zraw_r_km, "primary")
    if (!is.null(R55)) add("M1 (O1)", "exact window", R55$nsim, R55$p, NA, NA, "frozen primary window; full test")
    if (!is.null(R44)) for (k in seq_len(nrow(R44))) add(paste0("M2 (", R44$scenario[k], ")"), "fixed offset, simplified window", R44$nsim[k], R44$p_erl_nstep2[k], R44$zraw_max[k], R44$zraw_r_km[k], "assumed-strength offset")
    if (!is.null(R45)) add("M2 (O2)", "precision rerun", R45$nsim, R45$p_erl_nstep2, NA, NA, "nsim 10000; replaces the nsim-999 value")
    add("M3 (O1 trend + LGCP)", "simplified window", R46$nsim, R46$p_erl_nstep2, R46$zraw_max, R46$zraw_r_km, "adequate; hierarchy ends")
    T5 <- do.call(rbind, rows); wr(T5, "T5_generative_gof.csv"); show(T5) })
  
  cat("\n== T6 Stage 1 recovery (nsim 99, 100 replicates per cell) ==\n")
  safe("T6", {
    RC <- RECOVERY_CELLS
    T6 <- do.call(rbind, lapply(RC$cell, function(cc) { d <- A1[A1$cell == cc & !grepl("^failed", A1$class), ]; n <- nrow(d); tr <- RC$mechanism[RC$cell == cc]; r <- RC[RC$cell == cc, ]
    kp <- sum(d$class == "Poisson"); kl <- sum(d$class == "LGCP"); kg <- sum(d$class == "Geyer"); ki <- sum(d$class == "inadequate")
    main <- if (tr == "Geyer") kg else if (tr == "LGCP") kl else kp; ci <- wilson(main, n); cg <- wilson(kg, n)
    data.frame(cell = cc, truth = tr, n_target = r$n, sigma2 = r$sigma2, scale_km = r$scale_km, gamma = r$gamma, r_km = r$r_km, sat = r$sat, nuisance = r$nuisance, replicates = n, mean_n = mean(d$n),
               class_Poisson = kp, class_LGCP = kl, class_Geyer = kg, class_inadequate = ki, M1_rejection = mean(d$p1 < 0.05), reached_M3 = sum(!is.na(d$p3)), reached_M4 = sum(!is.na(d$p4)),
               P_correct_class = main / n, correct_lo = ci[1], correct_hi = ci[2], P_Geyer = kg / n, Geyer_lo = cg[1], Geyer_hi = cg[2],
               gamma_hat_median = suppressWarnings(median(d$gamma, na.rm = TRUE)), stringsAsFactors = FALSE) }))
    wr(T6, "T6_stage1_recovery.csv"); show(T6, 20)
    lg <- sum(T6$class_LGCP); cat(sprintf("\nLGCP-classed replicates: %d (from LGCP-truth cells %d, Geyer-truth cells %d, Poisson-truth cells %d)\n", lg,
                                          sum(T6$class_LGCP[T6$truth == "LGCP"]), sum(T6$class_LGCP[T6$truth == "Geyer"]), sum(T6$class_LGCP[T6$truth == "Poisson"]))) })
  
  cat("\n== T7a Robustness of M1 and M3 ==\n")
  safe("T7a", {
    a <- R49[, c("setting", "n", "nsim", "p", "mc_se", "flag", "ratio_2", "ratio_3", "ratio_5")]; a$test <- "M1"; a$note <- ""
    if (!is.null(R49b)) for (k in seq_len(nrow(R49b))) { i <- match(R49b$setting[k], a$setting); a$p[i] <- R49b$p_10000[k]; a$nsim[i] <- R49b$nsim[k]; a$mc_se[i] <- R49b$mc_se[k]; a$flag[i] <- R49b$decision[k]; a$note[i] <- paste("precision rerun; nsim 999 p", R49b$p_999[k]) }
    if (!is.null(R55)) a <- rbind(a, data.frame(setting = "R18_exact_window", n = R55$n, nsim = R55$nsim, p = R55$p, mc_se = R55$mc_se, flag = R55$flag, ratio_2 = NA, ratio_3 = NA, ratio_5 = NA, test = "M1", note = "frozen primary window", stringsAsFactors = FALSE))
    if (!is.null(R53)) a <- rbind(a, data.frame(setting = R53$setting, n = R53$n, nsim = R53$nsim, p = R53$p, mc_se = R53$mc_se, flag = R53$flag, ratio_2 = NA, ratio_3 = NA, ratio_5 = NA, test = "M3", note = paste0("sigma2 ", round(R53$sigma2, 2), ", alpha ", round(R53$alpha, 2), " km; M1 p ", R53$m1_p_script49), stringsAsFactors = FALSE))
    wr(a, "T7a_robustness_M1_M3.csv"); show(a, 30) })
  cat("\n== T7b Coordinate jitter (v2: rows looked up by type) ==\n")
  safe("T7b", {
    jb <- J50[match(c("base", "jitter"), J50$type), ]
    b1 <- data.frame(scenario = c("baseline (nominal coordinates)", "primary jitter (FIELD 10, AUTHOR 25, GIS 250 m)"), author_m = c(0, 25), gis_m = c(0, 250), reps = jb$reps, retained = jb$retained,
                     retention = jb$retention, ci_lo = jb$ci_lo, ci_hi = jb$ci_hi, p_median = jb$p_median, stringsAsFactors = FALSE)
    b2 <- data.frame(scenario = paste0("sensitivity ", J59$scenario, " (FIELD 10 m)"), author_m = J59$author_m, gis_m = J59$gis_m, reps = J59$reps, retained = J59$retained, retention = J59$retention,
                     ci_lo = J59$ci_lo, ci_hi = J59$ci_hi, p_median = J59$p_median, stringsAsFactors = FALSE)
    T7b <- rbind(b1, b2); T7b$criterion_0.80 <- T7b$retention >= 0.80; wr(T7b, "T7b_jitter.csv"); show(T7b) })
  cat("\n== T7c Block validation (first-order only) ==\n")
  safe("T7c", {
    pool <- function(d, v, w) sum(d[[v]] * d[[w]]) / sum(d[[w]])
    v1 <- do.call(rbind, lapply(split(V57, list(V57$block_km, V57$model)), function(d) { pg <- vapply(split(d, d$rep), pool, 1, v = "gain", w = "n_test"); pp <- vapply(split(d, d$rep), pool, 1, v = "pct", w = "n_test")
    data.frame(analysis = "script 57 draws (exploratory)", block_km = d$block_km[1], model = d$model[1], gain_median = median(pg), share_assignments_gain_gt0 = mean(pg > 0), percentile_median = median(pp), share_assignments_pct_gt_half = mean(pp > 0.5), stringsAsFactors = FALSE) }))
    v2 <- do.call(rbind, lapply(split(V56, V56$model), function(d) { pg <- vapply(split(d, d$rep), pool, 1, v = "gain_site_floor", w = "n_test_sites")
    data.frame(analysis = "script 56 (frozen design)", block_km = 60, model = d$model[1], gain_median = median(pg), share_assignments_gain_gt0 = mean(pg > 0), percentile_median = NA, share_assignments_pct_gt_half = NA, stringsAsFactors = FALSE) }))
    T7c <- rbind(v2, v1); wr(T7c, "T7c_validation.csv"); show(T7c, 30) })
  cat("\n== T7d Anisotropy diagnostic ==\n")
  safe("T7d", { T7d <- AN58; wr(T7d, "T7d_anisotropy.csv"); show(T7d) })
  
  cat("\n== T8 Sensitivity of the residual dependence to observation, chronology and provenance assumptions ==\n")
  safe("T8", {
    pick <- function(s, lab) { i <- match(s, R49$setting); data.frame(assumption = lab, setting = s, n = R49$n[i], nsim = R49$nsim[i], p = R49$p[i], ratio_2 = R49$ratio_2[i], ratio_3 = R49$ratio_3[i], ratio_5 = R49$ratio_5[i], note = R49$flag[i], stringsAsFactors = FALSE) }
    t8 <- rbind(pick("R00_baseline", "primary control (new seed)"), pick("R15_zS2", "observation proxy, free coefficient (access/open)"), pick("R16_zS3", "observation proxy, free coefficient (adds pattern-informed core term)"),
                pick("R01_redep", "redeposited sites excluded"), pick("R03_LP", "chronology: Lower Palaeolithic (any phase)"), pick("R04_MP", "chronology: Middle Palaeolithic (any phase)"), pick("R05_UP", "chronology: Upper Palaeolithic (any phase)"))
    t8 <- rbind(t8, data.frame(assumption = "Dataset C (non-fluvial)", setting = "R02_datasetC", n = 97, nsim = R49b$nsim[R49b$setting == "R02_datasetC"], p = R49b$p_10000[R49b$setting == "R02_datasetC"], ratio_2 = R49$ratio_2[R49$setting == "R02_datasetC"], ratio_3 = R49$ratio_3[R49$setting == "R02_datasetC"], ratio_5 = R49$ratio_5[R49$setting == "R02_datasetC"], note = "precision rerun", stringsAsFactors = FALSE))
    for (k in seq_len(nrow(R44))) t8 <- rbind(t8, data.frame(assumption = paste("fixed observation offset", R44$scenario[k]), setting = paste0("M2-", R44$scenario[k]), n = 197, nsim = R44$nsim[k], p = R44$p_erl_nstep2[k], ratio_2 = NA, ratio_3 = NA, ratio_5 = NA, note = "assumed strength", stringsAsFactors = FALSE))
    t8 <- rbind(t8, data.frame(assumption = "GIS_MAP-coordinate sites excluded (exploratory, post hoc)", setting = "R19_noGIS", n = G60$n, nsim = G60$nsim, p = G60$p, ratio_2 = G60$ratio_2, ratio_3 = G60$ratio_3, ratio_5 = G60$ratio_5, note = G60$flag, stringsAsFactors = FALSE))
    for (k in seq_len(nrow(L62))) t8 <- rbind(t8, data.frame(assumption = paste("lead author removed (exploratory):", L62$author_removed[k]), setting = paste0("LOAO_", L62$author_removed[k]), n = L62$n[k], nsim = L62$nsim[k], p = L62$p[k], ratio_2 = L62$ratio_2[k], ratio_3 = L62$ratio_3[k], ratio_5 = L62$ratio_5[k], note = L62$flag[k], stringsAsFactors = FALSE))
    for (k in seq_len(nrow(C63))) t8 <- rbind(t8, data.frame(assumption = paste("matched control (exploratory):", C63$description[k]), setting = paste0("control_", C63$group[k]), n = C63$n_median[k], nsim = 199, p = C63$p_median[k], ratio_2 = C63$pair_ratio_2[k], ratio_3 = C63$pair_ratio_3[k], ratio_5 = C63$pair_ratio_5[k], note = sprintf("retention %.2f (%.2f-%.2f); p is the median over replicates", C63$retention[k], C63$ci_lo[k], C63$ci_hi[k]), stringsAsFactors = FALSE))
    wr(t8, "T8_assumption_sensitivity.csv"); show(t8, 30) })
  
  cat("\n== Consistency checks against values reported in the analysis conversation (v2 adds table-level T7b checks) ==\n")
  safe("checks", {
    ck <- function(name, exp, found, tol = 1e-3) data.frame(check = name, expected = exp, found = found, ok = isTRUE(abs(exp - found) <= tol), stringsAsFactors = FALSE)
    bl <- T7b$scenario == "baseline (nominal coordinates)"; pj <- T7b$scenario == "primary jitter (FIELD 10, AUTHOR 25, GIS 250 m)"
    CK <- rbind(ck("M1 p (script 42)", 0.0310, g1(R42$p_erl_nstep2)), ck("exact-window M1 p", 0.0180, g1(R55$p)), ck("M2-O2 p at nsim 10000", 0.0387, g1(R45$p_erl_nstep2), 5e-4),
                ck("M2-O3 p", 0.0240, g1(R44$p_erl_nstep2[R44$scenario == "O3"])), ck("M3 p", 0.1920, g1(R46$p_erl_nstep2)), ck("M3 sigma2_hat", 1.2298, g1(R46$sigma2_hat), 1e-3), ck("M3 alpha_hat", 6.048, g1(R46$alpha_hat), 1e-3),
                ck("jitter retention (primary), source file", 0.94, g1(J50$retention[J50$type == "jitter"])), ck("jitter retention (baseline), source file", 0.92, g1(J50$retention[J50$type == "base"])), ck("jitter S3 retention", 0.84, g1(J59$retention[J59$scenario == "S3"])),
                ck("T7b table: baseline reps", 100, g1(T7b$reps[bl]), 0), ck("T7b table: baseline retention", 0.92, g1(T7b$retention[bl])), ck("T7b table: primary jitter reps", 200, g1(T7b$reps[pj]), 0), ck("T7b table: primary jitter retention", 0.94, g1(T7b$retention[pj])),
                ck("R16 free z_S3 p", 0.096, g1(R49$p[R49$setting == "R16_zS3"])), ck("R04 MP p", 0.397, g1(R49$p[R49$setting == "R04_MP"])), ck("R12 precision p", 0.0290, g1(R49b$p_10000[R49b$setting == "R12_quad05"]), 5e-4),
                ck("R02 precision p", 0.0461, g1(R49b$p_10000[R49b$setting == "R02_datasetC"]), 5e-4), ck("Stage 1 rows", 1900, nrow(A1), 0), ck("Stage 1 Geyer-classed replicates in Geyer-truth cells", 3, sum(A1$class == "Geyer" & A1$cell %in% c("G1", "G2", "G3", "G4", "C1", "C3")), 0),
                ck("GIS-excluded M1 p", 0.1270, g1(G60$p)), ck("matched control retention G", 0.45, g1(C63$retention[C63$group == "G"])), ck("matched control retention R", 0.75, g1(C63$retention[C63$group == "R"])),
                ck("leave-one-author-out authors", 5, nrow(L62), 0), ck("sites in dataset A", 197, nrow(SITES), 0))
    print(CK, row.names = FALSE); cat(sprintf("\nchecks passed %d of %d\n", sum(CK$ok), nrow(CK))); wr(CK, "64_consistency_checks.csv") })
  
  cat("\n== Missing inputs ==\n"); cat(if (length(ENVM$missing)) paste(ENVM$missing, collapse = "\n") else "none", "\n")
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT, "64_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")