# =============================================================================
# 15_freeze_v1.R   (FINAL: author confirmed 2026-11-07)
# Stage: pre-analysis freeze | no models fitted
# Purpose: write freeze v1 (parameter file, human-readable freeze document, hash
#   manifest of frozen inputs, derived layers and scripts).
# Safeguard: refuses to overwrite an existing FROZEN manifest.
# Outputs: R/freeze_v1_params.R, docs/freeze_v1.md, docs/freeze_v1_manifest.csv
# Not executed by assistant. Run top to bottom in RStudio. Hashing takes a few minutes.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
AUTHOR_CONFIRMED_DATE <- "2026-11-07"

while (sink.number() > 0) sink()
source(REPO_CONFIG)
dir.create(file.path(REPO_DIR, "docs"), showWarnings = FALSE)
params_f <- file.path(REPO_DIR, "R", "freeze_v1_params.R")
doc_f    <- file.path(REPO_DIR, "docs", "freeze_v1.md")
man_f    <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv")

if (file.exists(man_f)) {
  old <- tryCatch(read.csv(man_f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(old) && "freeze_status" %in% names(old) && any(grepl("^FROZEN", old$freeze_status)))
    stop("A FROZEN freeze_v1 manifest already exists (", unique(old$freeze_status)[1], "). ",
         "Not overwriting. Log a deviation in docs/protocol_deviation_log.csv; to redo on purpose, delete the manifest by hand.")
}
status <- if (is.na(AUTHOR_CONFIRMED_DATE)) "PROPOSED (not frozen until the author confirms)" else
  paste("FROZEN on", AUTHOR_CONFIRMED_DATE, "(author-confirmed)")
DELEGATION <- "Author delegated routine methodological and design decisions to the assistant on 2026-11-07. Each decision is logged in the decision register as APPROVED (author-delegated). Integrity rules unchanged: frozen items change only through docs/protocol_deviation_log.csv; no site is added, dropped or altered to obtain a result; post-hoc choices are labelled exploratory."

# ---- PARAMETERS ---------------------------------------------------------------
FREEZE_VERSION <- "v1"; FREEZE_STATUS <- status; FREEZE_DELEGATION <- DELEGATION
N_ANALYTICAL <- 197L
DATASETS <- list(primary = "A (all 197)", sensitivity = c("C non-fluvial (n=97)", "A minus redeposited IDs 48,49 (n=195)"))
REDEPOSITED_IDS <- c(48L, 49L)
CHRONOLOGY_SETS <- c(anyLP = 75L, anyMP = 114L, anyUP = 79L)
WINDOW_AREA_KM2 <- 21262.3
WINDOW <- list(primary = "exact boundary polygon, EPSG:32644, 21262.3 km2", fallback_simplify_m = 100,
               erode_sensitivity_km = 5, erode_n_sites = 169L)
PRIMARY_COORDINATES <- "nominal (best-supported) coordinates; coordinate jitter is robustness only"
GRIDS_M <- list(primary = 500, sensitivity = c(250, 1000), nested_origin_m = 1000)
QUADRATURE_KM <- list(primary = 1, sensitivity = 0.5)
E1 <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5")
COVARIATE_SETS <- list(E1 = E1, E2 = c(E1, "lith5"), E3 = c(E1, "lith5", "geomorph3"),
                       E4 = c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "lith5"))
COVARIATE_DF <- c(E1 = 6, E2 = 10, E3 = 12, E4 = 8)
PRIMARY_SET <- "E3"
REFERENCE_LEVELS <- c(lith5 = "crystalline_felsic", geomorph3 = "pediplain")
OBS <- list(primary = "O1", sensitivity = c("O2", "O3"), free_variants = c("z_S2", "z_S3"),
            kappa = c(O2 = 1.1192, O3 = 2.5366),
            rule = "sufficiency under O2 stops the hierarchy; sufficiency only under O3 = observation-sensitive, hierarchy continues")
R_MIN_KM <- c(A = 1, C = 2); R_MAX_DECISION_KM <- 10; R_MAX_EXPLORE_KM <- 20; R_STEP_KM <- 0.5
PCF_STOYAN_C <- 0.15
N_BY_SET <- c(A = 197, C = 97)
LAMBDA_PER_KM2 <- N_BY_SET / WINDOW_AREA_KM2
PCF_BW_KM <- PCF_STOYAN_C / sqrt(LAMBDA_PER_KM2)
GINHOM_LOWER_KM <- pmax(R_MIN_KM[c("A", "C")], ceiling(PCF_BW_KM / 0.5) * 0.5)
GEYER_R_KM <- c(2, 4, 7, 10); GEYER_SAT <- c(2, 4)
SUMMARIES <- list(primary = c("Linhom", "ginhom"), secondary = c("G", "F", "J"), joint_test = "ERL combined global envelope",
                  edge_correction = "isotropic", edge_sensitivity = c("translate", "border"), pcf_stoyan_c = PCF_STOYAN_C,
                  ginhom_range_rule = "lower bound = max(R_MIN, kernel bandwidth rounded up to 0.5 km); upper = R_MAX_DECISION",
                  intensity_for_inhom = "trend-only E3 Poisson fit, refit to each pattern")
TESTS <- list(alpha = 0.05, nsim_empirical = 999L, nsim_empirical_mcmc = 499L, nsim_recovery = 199L,
              nsim_recovery_floor = 99L, recovery_reps = 200L, recovery_reps_floor = 100L, jitter_reps = 200L)
JITTER_SIGMA_M <- c(FIELD = 10, AUTHOR_REPORTED = 25, GIS_MAP = 250)
JITTER_SENS_M  <- list(AUTHOR_REPORTED = c(10, 100, 250), GIS_MAP = c(100, 500, 1000))
ESTIMATORS <- c(M1 = "Poisson maximum pseudolikelihood / Berman-Turner (spatstat ppm)", M2 = "as M1 with fixed log E offset or free z_S",
                M3 = "LGCP, exponential covariance, composite likelihood (spatstat kppm, method clik2)",
                M4 = "Geyer saturation, maximum pseudolikelihood, profile search over GEYER_R_KM x GEYER_SAT",
                note = "model classes frozen; argument lists verified against documentation before coding")
DECISION <- list(deficiency = "combined ERL p < 0.05", reduction_min = 0.50, power_min = 0.80, fpr_max = 0.10,
                 mechanism_recovery_min = 0.70, jitter_retention_min = 0.80,
                 required_settings = c("O2", "E2", "E4", "grid250", "grid1000", "eroded_window", "redeposited_excluded"),
                 untestable_ok = c("Dataset C", "chronology subsets"))
VALIDATION <- list(block_km = 60, folds = 5L, buffer_km = 10)
SEEDS <- list(master = MASTER_SEED, offsets = c(recovery = 10000, empirical = 20000, jitter = 30000, validation = 40000, bootstrap = 50000),
              env_truth_offset = 5000)
RECOVERY_ENV <- list(coef_sd = 0.3, levels = c(E_weak = 0.5, E_mod = 1), seed_offset = 5000, expected_n = 197)
cell <- function(cell, mechanism, env = "E_mod", n = 197, sigma2 = NA, scale_km = NA, gamma = NA, r_km = NA, sat = NA, nuisance = "none")
  data.frame(cell, mechanism, env, n, sigma2, scale_km, gamma, r_km, sat, nuisance, stringsAsFactors = FALSE)
RECOVERY_CELLS <- rbind(
  cell("P1", "Poisson", env = "E_weak"), cell("P2", "Poisson"),
  cell("L1", "LGCP", sigma2 = 0.5, scale_km = 3), cell("L2", "LGCP", sigma2 = 0.5, scale_km = 8),
  cell("L3", "LGCP", sigma2 = 1.5, scale_km = 3), cell("L4", "LGCP", sigma2 = 1.5, scale_km = 8),
  cell("G1", "Geyer", gamma = 1.5, r_km = 3, sat = 3), cell("G2", "Geyer", gamma = 1.5, r_km = 7, sat = 3),
  cell("G3", "Geyer", gamma = 2.5, r_km = 3, sat = 3), cell("G4", "Geyer", gamma = 2.5, r_km = 7, sat = 3),
  cell("N1", "Poisson", nuisance = "true O3 offset, fitted without it"),
  cell("N2", "Poisson", nuisance = "superposition of two independent patterns (100+97)"),
  cell("N3", "Poisson", nuisance = "coordinate jitter at primary sigma"),
  cell("N4", "Poisson", nuisance = "misspecified truth: quadratic slope + slope x lithology interaction"),
  cell("N5", "LGCP", sigma2 = 0.5, scale_km = 8, nuisance = "LGCP + true O2 offset + jitter, no Gibbs"),
  cell("C1", "Geyer", n = 75, gamma = 2.5, r_km = 3, sat = 3), cell("C2", "LGCP", n = 75, sigma2 = 1.5, scale_km = 8),
  cell("C3", "Geyer", n = 114, gamma = 2.5, r_km = 3, sat = 3), cell("C4", "LGCP", n = 114, sigma2 = 1.5, scale_km = 8))

PARAM_NAMES <- c("FREEZE_VERSION", "FREEZE_STATUS", "FREEZE_DELEGATION", "N_ANALYTICAL", "DATASETS", "REDEPOSITED_IDS",
                 "CHRONOLOGY_SETS", "WINDOW_AREA_KM2", "WINDOW", "PRIMARY_COORDINATES", "GRIDS_M", "QUADRATURE_KM", "COVARIATE_SETS",
                 "COVARIATE_DF", "PRIMARY_SET", "REFERENCE_LEVELS", "OBS", "R_MIN_KM", "R_MAX_DECISION_KM", "R_MAX_EXPLORE_KM", "R_STEP_KM",
                 "PCF_STOYAN_C", "N_BY_SET", "LAMBDA_PER_KM2", "PCF_BW_KM", "GINHOM_LOWER_KM", "GEYER_R_KM", "GEYER_SAT", "SUMMARIES", "TESTS",
                 "JITTER_SIGMA_M", "JITTER_SENS_M", "ESTIMATORS", "DECISION", "VALIDATION", "SEEDS", "RECOVERY_ENV", "RECOVERY_CELLS")
writeLines(c(sprintf("# freeze_v1_params.R : %s | written %s", status, format(Sys.time())),
             "# Source this file in confirmatory scripts. Do not edit; changes go through docs/protocol_deviation_log.csv."), params_f)
dump(PARAM_NAMES, file = params_f, append = TRUE, control = c("keepNA", "keepInteger", "niceNames", "showAttributes"))

# ---- HUMAN-READABLE DOCUMENT -------------------------------------------------
doc <- c(
  "# Freeze v1: Beyond Environmental Suitability (Nagpur-Chandrapur)", "", paste("Status:", status), "",
  "Protocol of record: Research_Design.docx. Freeze v1 fixes the choices below BEFORE any archaeological point-process model is fitted. Any later change is logged in docs/protocol_deviation_log.csv.", "",
  paste("Delegation:", DELEGATION), "",
  "## 1. Data and window",
  "- NA = 197 (Gate 1). Primary dataset A. Sensitivity: C non-fluvial (97); A without redeposited IDs 48, 49 (195). Dataset B retired.",
  "- Primary analysis uses nominal (best-supported) coordinates. Coordinate jitter is robustness only.",
  "- Chronology: pooled primary; sensitivity on overlapping any-phase sets LP 75, MP 114, UP 79; pure sets (48/62/27) descriptive only.",
  "- Window: exact boundary polygon; fallback 100 m simplified polygon; sensitivity eroded by 5 km (169 sites).", "",
  "## 2. Covariates and grids",
  "- Primary baseline E3 (12 df). Robustness E1 (6), E2 (10), E4 (8). Reference levels: crystalline felsic, pediplain.",
  "- Grids 500 m primary, 250 and 1000 m sensitivity. Quadrature spacing 1 km, sensitivity 0.5 km. Unwinsorised-TPI sensitivity (E3).", "",
  "## 3. Observation",
  "- Primary O1. O2, O3 fixed offsets (kappa 1.1192 and 2.5366, calibrated at 500 m to 5th-95th percentile contrast 2 and 4). Free-coefficient variants z_S2, z_S3.",
  "- Sufficiency under O2 stops the hierarchy; sufficiency only under O3 is labelled observation-sensitive and the hierarchy continues.", "",
  "## 4. Scales, summaries, tests",
  sprintf("- L_inhom decision range: %.0f-%.0f km (A), %.0f-%.0f km (C). g_inhom lower bound = max(r_min, kernel bandwidth rounded up to 0.5 km): %.1f km (A, bw %.2f km), %.1f km (C, bw %.2f km); upper %.0f km. Exploratory plots to %.0f km; step %.1f km.",
          R_MIN_KM["A"], R_MAX_DECISION_KM, R_MIN_KM["C"], R_MAX_DECISION_KM, GINHOM_LOWER_KM["A"], PCF_BW_KM["A"], GINHOM_LOWER_KM["C"], PCF_BW_KM["C"], R_MAX_DECISION_KM, R_MAX_EXPLORE_KM, R_STEP_KM),
  "- Primary summaries L_inhom and g_inhom tested jointly (ERL combined global envelope, alpha 0.05). Secondary G, F, J. Edge correction isotropic; sensitivity translation, border.",
  "- Inhomogeneous intensity = trend-only E3 Poisson fit refit to each pattern. pcf bandwidth: Stoyan rule c = 0.15 (default value verified against documentation at coding).",
  "- Interaction grid: r in {2,4,7,10} km x saturation in {2,4}, profile pseudolikelihood; the recovery study repeats the grid search.",
  "- Estimators: M1/M2 Poisson MPL; M3 LGCP (exponential, kppm method clik2); M4 Geyer MPL with profile search. Classes frozen; argument lists verified before coding.",
  "- Simulations: empirical 999 (M1, M3), 499 (M4); recovery 199 per test, 200 replicates per cell, floors 99 and 100 if a timing pilot requires.", "",
  "## 5. Recovery design (19 cells)",
  "- Environment truth: fixed draw of 12 coefficients N(0, 0.3^2) per SD, seed MASTER_SEED + 5000; E_weak = 0.5x, E_mod = 1x; expected n 197.",
  "- Poisson (2), LGCP sigma2 {0.5,1.5} x scale {3,8} km (4), Geyer gamma {1.5,2.5} x r {3,7} km with sat 3 (4), nuisance (5), n = 75 and 114 cells (4).",
  "- Nuisance cells: O3 truth fitted without offset; superposition of two patterns; primary jitter; misspecified truth; LGCP + O2 + jitter.", "",
  "## 6. Decision rules",
  "- Deficiency: combined ERL p < 0.05. M4 substantially reduces a deficiency: p >= 0.05 and maximum standardised deviation reduced by at least 50% relative to M3.",
  "- Gate 6: power >= 0.80, false-positive rate <= 0.10 in every non-Gibbs cell, mechanism recovery >= 0.70 (LGCP and Geyer), in cells bracketing the empirical estimates. Extra bracketing cells allowed if estimates fall outside the grid (logged as extension).",
  "- Interaction claim also needs retention in >= 80% of jitter replicates (200) and in O2, E2, E4, 250 m, 1000 m, eroded window, redeposited-excluded. Dataset C and chronology subsets: untestable is not failure.",
  "- Coordinate jitter sigma (m): FIELD 10, AUTHOR_REPORTED 25, GIS_MAP 250; sensitivity AUTHOR 10/100/250, GIS 100/500/1000.",
  "- Validation (not interaction evidence): spatial blocks 60 km, 5 folds, buffer 10 km.",
  "- Seeds: MASTER_SEED + offsets recovery 10000, empirical 20000, jitter 30000, validation 40000, bootstrap 50000, plus index.", "",
  "## 7. Order of execution",
  "1 timing pilot; 2 recovery study (Gate 6); 3 empirical M1..M4 with hierarchy; 4 robustness; 5 interpretation. Empirical g/K/L plots are exploratory until step 2 is complete.", "",
  "## 8. Disclosed non-blind elements",
  "- O3 core-position term chosen after the descriptive edge deficit was seen (7 sites observed vs 20 expected within 2 km of the edge).",
  "- Saturation grid {2,4} brackets crude pair counts from script 07 (no model fits). Radius set follows nearest-neighbour distances and edge erosion.",
  "- Geometry and descriptive counts only; no site-by-covariate tables and no model fits informed any frozen choice.")
writeLines(doc, doc_f)

# ---- HASH MANIFEST ---------------------------------------------------------------
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
keep <- !(man$domain %in% c("sentinel2_band", "sentinel2_index", "soil")) & man$layer != "MineralDeposit"
raw_paths <- file.path(RAW_DIR, man$path[keep])
shp  <- raw_paths[grepl("\\.shp$", raw_paths, ignore.case = TRUE)]
side <- unlist(lapply(shp, function(p) paste0(sub("\\.shp$", "", p, ignore.case = TRUE), c(".shp", ".dbf", ".shx", ".prj", ".cpg"))))
raw_all <- unique(c(raw_paths, side))
env_layers <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5", "lith5", "geomorph3",
                "logE_O1", "logE_O2", "logE_O3", "z_S2", "z_S3", "inwin", "incenter")
env_files <- unlist(lapply(c(250, 500, 1000), function(r) file.path(REPO_DIR, "data_derived", "env", paste0(r, "m"), paste0(env_layers, ".tif"))))
der_files <- file.path(REPO_DIR, c("data_derived/sites_NA.csv", "data_derived/sites_master_decided.csv", "data_derived/ppp_NA_window.rds",
                                   "docs/site_review_DECIDED.csv", "docs/site_inclusion_exclusion_log.csv", "docs/covariate_definitions.csv",
                                   "docs/observation_scenarios.csv", "docs/input_manifest.csv", "docs/external_download_log.csv", "docs/protocol_deviation_log.csv",
                                   "docs/freeze_v1.md"))
scr <- list.files(REPO_DIR, pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
scr <- scr[!grepl("/renv/", scr)]
groups <- list(raw_input = raw_all, derived_layer = env_files, derived_table = der_files, script = scr)
hash <- do.call(rbind, lapply(names(groups), function(g) {
  f <- groups[[g]]; ex <- file.exists(f)
  if (any(!ex)) cat("MISSING (", g, "):", paste(basename(f[!ex]), collapse = ", "), "\n")
  f <- f[ex]; if (!length(f)) return(NULL)
  cat("hashing", length(f), g, "files ...\n"); flush.console()
  data.frame(group = g, path = f, bytes = file.size(f), md5 = unname(tools::md5sum(f)), stringsAsFactors = FALSE)
}))
hash$freeze_status <- status
write.csv(hash, man_f, row.names = FALSE)

cat("\n==== FREEZE v1 ====\nStatus:", status, "\n")
cat("Files written:\n ", params_f, "\n ", doc_f, "\n ", man_f, "\n")
print(aggregate(bytes ~ group, data = hash, FUN = function(x) c(files = length(x), MB = round(sum(x) / 1e6, 1))))
cat("Recovery cells:", nrow(RECOVERY_CELLS), "| Geyer fit grid:", length(GEYER_R_KM) * length(GEYER_SAT), "combinations\n")
cat("df by set:", paste(names(COVARIATE_DF), COVARIATE_DF, collapse = ", "), "\n")
cat(sprintf("Mean intensity per km2: A %.5f, C %.5f | pcf bandwidth (km): A %.2f, C %.2f | g_inhom lower bound (km): A %.1f, C %.1f\n",
            LAMBDA_PER_KM2["A"], LAMBDA_PER_KM2["C"], PCF_BW_KM["A"], PCF_BW_KM["C"], GINHOM_LOWER_KM["A"], GINHOM_LOWER_KM["C"]))
cat(if (is.na(AUTHOR_CONFIRMED_DATE)) "\nNEXT: confirm in chat, set AUTHOR_CONFIRMED_DATE, re-run.\n" else "\nFROZEN. Later changes go to docs/protocol_deviation_log.csv.\n")
writeLines(capture.output(sessionInfo()), file.path(REPO_DIR, "docs", "sessionInfo_15.txt"))