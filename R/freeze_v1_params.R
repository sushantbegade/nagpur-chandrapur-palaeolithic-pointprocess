# freeze_v1_params.R : FROZEN on 2026-11-07 (author-confirmed) | written 2026-11-07 12:28:34
# Source this file in confirmatory scripts. Do not edit; changes go through docs/protocol_deviation_log.csv.
FREEZE_VERSION <-
"v1"
FREEZE_STATUS <-
"FROZEN on 2026-11-07 (author-confirmed)"
FREEZE_DELEGATION <-
"Author delegated routine methodological and design decisions to the assistant on 2026-11-07. Each decision is logged in the decision register as APPROVED (author-delegated). Integrity rules unchanged: frozen items change only through docs/protocol_deviation_log.csv; no site is added, dropped or altered to obtain a result; post-hoc choices are labelled exploratory."
N_ANALYTICAL <-
197L
DATASETS <-
list(primary = "A (all 197)", sensitivity = c("C non-fluvial (n=97)", 
"A minus redeposited IDs 48,49 (n=195)"))
REDEPOSITED_IDS <-
48:49
CHRONOLOGY_SETS <-
c(anyLP = 75L, anyMP = 114L, anyUP = 79L)
WINDOW_AREA_KM2 <-
21262.3
WINDOW <-
list(primary = "exact boundary polygon, EPSG:32644, 21262.3 km2", 
    fallback_simplify_m = 100, erode_sensitivity_km = 5, erode_n_sites = 169L)
PRIMARY_COORDINATES <-
"nominal (best-supported) coordinates; coordinate jitter is robustness only"
GRIDS_M <-
list(primary = 500, sensitivity = c(250, 1000), nested_origin_m = 1000)
QUADRATURE_KM <-
list(primary = 1, sensitivity = 0.5)
COVARIATE_SETS <-
list(E1 = c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", 
"z_ldist_riv3", "z_ldist_riv5"), E2 = c("z_elev_p1", "z_elev_p2", 
"z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5", "lith5"
), E3 = c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", 
"z_ldist_riv5", "lith5", "geomorph3"), E4 = c("z_elev_p1", "z_elev_p2", 
"z_lslope", "z_tpi2km", "lith5"))
COVARIATE_DF <-
c(E1 = 6, E2 = 10, E3 = 12, E4 = 8)
PRIMARY_SET <-
"E3"
REFERENCE_LEVELS <-
c(lith5 = "crystalline_felsic", geomorph3 = "pediplain")
OBS <-
list(primary = "O1", sensitivity = c("O2", "O3"), free_variants = c("z_S2", 
"z_S3"), kappa = c(O2 = 1.1192, O3 = 2.5366), rule = "sufficiency under O2 stops the hierarchy; sufficiency only under O3 = observation-sensitive, hierarchy continues")
R_MIN_KM <-
c(A = 1, C = 2)
R_MAX_DECISION_KM <-
10
R_MAX_EXPLORE_KM <-
20
R_STEP_KM <-
0.5
PCF_STOYAN_C <-
0.15
N_BY_SET <-
c(A = 197, C = 97)
LAMBDA_PER_KM2 <-
c(A = 0.00926522530488235, C = 0.00456206525164258)
PCF_BW_KM <-
c(A = 1.55834376155835, C = 2.22080544040303)
GINHOM_LOWER_KM <-
c(A = 2, C = 2.5)
GEYER_R_KM <-
c(2, 4, 7, 10)
GEYER_SAT <-
c(2, 4)
SUMMARIES <-
list(primary = c("Linhom", "ginhom"), secondary = c("G", "F", 
"J"), joint_test = "ERL combined global envelope", edge_correction = "isotropic", 
    edge_sensitivity = c("translate", "border"), pcf_stoyan_c = 0.15, 
    ginhom_range_rule = "lower bound = max(R_MIN, kernel bandwidth rounded up to 0.5 km); upper = R_MAX_DECISION", 
    intensity_for_inhom = "trend-only E3 Poisson fit, refit to each pattern")
TESTS <-
list(alpha = 0.05, nsim_empirical = 999L, nsim_empirical_mcmc = 499L, 
    nsim_recovery = 199L, nsim_recovery_floor = 99L, recovery_reps = 200L, 
    recovery_reps_floor = 100L, jitter_reps = 200L)
JITTER_SIGMA_M <-
c(FIELD = 10, AUTHOR_REPORTED = 25, GIS_MAP = 250)
JITTER_SENS_M <-
list(AUTHOR_REPORTED = c(10, 100, 250), GIS_MAP = c(100, 500, 
1000))
ESTIMATORS <-
c(M1 = "Poisson maximum pseudolikelihood / Berman-Turner (spatstat ppm)", 
M2 = "as M1 with fixed log E offset or free z_S", M3 = "LGCP, exponential covariance, composite likelihood (spatstat kppm, method clik2)", 
M4 = "Geyer saturation, maximum pseudolikelihood, profile search over GEYER_R_KM x GEYER_SAT", 
note = "model classes frozen; argument lists verified against documentation before coding"
)
DECISION <-
list(deficiency = "combined ERL p < 0.05", reduction_min = 0.5, 
    power_min = 0.8, fpr_max = 0.1, mechanism_recovery_min = 0.7, 
    jitter_retention_min = 0.8, required_settings = c("O2", "E2", 
    "E4", "grid250", "grid1000", "eroded_window", "redeposited_excluded"
    ), untestable_ok = c("Dataset C", "chronology subsets"))
VALIDATION <-
list(block_km = 60, folds = 5L, buffer_km = 10)
SEEDS <-
list(master = 20260101, offsets = c(recovery = 10000, empirical = 20000, 
jitter = 30000, validation = 40000, bootstrap = 50000), env_truth_offset = 5000)
RECOVERY_ENV <-
list(coef_sd = 0.3, levels = c(E_weak = 0.5, E_mod = 1), seed_offset = 5000, 
    expected_n = 197)
RECOVERY_CELLS <-
structure(list(cell = c("P1", "P2", "L1", "L2", "L3", "L4", "G1", 
"G2", "G3", "G4", "N1", "N2", "N3", "N4", "N5", "C1", "C2", "C3", 
"C4"), mechanism = c("Poisson", "Poisson", "LGCP", "LGCP", "LGCP", 
"LGCP", "Geyer", "Geyer", "Geyer", "Geyer", "Poisson", "Poisson", 
"Poisson", "Poisson", "LGCP", "Geyer", "LGCP", "Geyer", "LGCP"
), env = c("E_weak", "E_mod", "E_mod", "E_mod", "E_mod", "E_mod", 
"E_mod", "E_mod", "E_mod", "E_mod", "E_mod", "E_mod", "E_mod", 
"E_mod", "E_mod", "E_mod", "E_mod", "E_mod", "E_mod"), n = c(197, 
197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 197, 
197, 75, 75, 114, 114), sigma2 = c(NA, NA, 0.5, 0.5, 1.5, 1.5, 
NA, NA, NA, NA, NA, NA, NA, NA, 0.5, NA, 1.5, NA, 1.5), scale_km = c(NA, 
NA, 3, 8, 3, 8, NA, NA, NA, NA, NA, NA, NA, NA, 8, NA, 8, NA, 
8), gamma = c(NA, NA, NA, NA, NA, NA, 1.5, 1.5, 2.5, 2.5, NA, 
NA, NA, NA, NA, 2.5, NA, 2.5, NA), r_km = c(NA, NA, NA, NA, NA, 
NA, 3, 7, 3, 7, NA, NA, NA, NA, NA, 3, NA, 3, NA), sat = c(NA, 
NA, NA, NA, NA, NA, 3, 3, 3, 3, NA, NA, NA, NA, NA, 3, NA, 3, 
NA), nuisance = c("none", "none", "none", "none", "none", "none", 
"none", "none", "none", "none", "true O3 offset, fitted without it", 
"superposition of two independent patterns (100+97)", "coordinate jitter at primary sigma", 
"misspecified truth: quadratic slope + slope x lithology interaction", 
"LGCP + true O2 offset + jitter, no Gibbs", "none", "none", "none", 
"none")), row.names = c(NA, -19L), class = "data.frame")
