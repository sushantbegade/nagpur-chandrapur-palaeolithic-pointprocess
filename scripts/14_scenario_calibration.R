# =============================================================================
# 14_scenario_calibration.R
# Stage: 04_observation | Gate 5 closure
# Purpose: calibrate O2/O3 offset strength to a 5th-95th percentile contrast of E
#   over window-centre cells at the primary grid, apply the same kappa at all grids,
#   rewrite logE_O2 and logE_O3, record parameters.
# Does NOT read site data. No models.
# Requires outputs of 13_observation_proxies.R.
# Outputs: data_derived/env/<res>m/{logE_O2,logE_O3}.tif (overwritten),
#   docs/observation_scenarios.csv, 04_observation/outputs/scenario_calibration_report.txt,
#   scenario_range_v2.csv, sessionInfo_14.txt
# Not executed by assistant. Run top to bottom in RStudio.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
RES_SET <- c(250, 500, 1000); PRIMARY <- 500
CONTRAST_O2 <- 2; CONTRAST_O3 <- 4; QLO <- 0.05; QHI <- 0.95

while (sink.number() > 0) sink()
source(REPO_CONFIG)
suppressPackageStartupMessages(library(terra))
out <- dir_of("04_observation"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
env_dir <- file.path(REPO_DIR, "data_derived", "env")
for (r in RES_SET) stopifnot(all(file.exists(file.path(env_dir, paste0(r, "m"),
                                                       paste0(c("S2", "S3", "inwin", "incenter"), ".tif")))))
rdl <- function(res, n) rast(file.path(env_dir, paste0(res, "m"), paste0(n, ".tif")))

sink(file.path(out, "scenario_calibration_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== CALIBRATION AT", PRIMARY, "m ====\n")
  ic0 <- !is.na(values(rdl(PRIMARY, "incenter"), mat = FALSE))
  S2 <- values(rdl(PRIMARY, "S2"), mat = FALSE)[ic0]; S3 <- values(rdl(PRIMARY, "S3"), mat = FALSE)[ic0]
  q2 <- quantile(S2, c(QLO, QHI), na.rm = TRUE); q3 <- quantile(S3, c(QLO, QHI), na.rm = TRUE)
  kap2 <- log(CONTRAST_O2) / (q2[2] - q2[1]); kap3 <- log(CONTRAST_O3) / (q3[2] - q3[1])
  cat(sprintf("S2 p5 %.4f p95 %.4f | S3 p5 %.4f p95 %.4f\n", q2[1], q2[2], q3[1], q3[2]))
  cat(sprintf("kappa_O2 = %.4f (target p95/p5 = %g) | kappa_O3 = %.4f (target p95/p5 = %g)\n",
              kap2, CONTRAST_O2, kap3, CONTRAST_O3))
  cat(sprintf("Implied factor over full proxy range 0 -> 1: O2 %.2f | O3 %.2f\n", exp(kap2), exp(kap3)))
  
  rows <- list()
  for (res in RES_SET) {
    d <- file.path(env_dir, paste0(res, "m"))
    inw <- rdl(res, "inwin"); ic <- !is.na(values(rdl(res, "incenter"), mat = FALSE))
    for (sc in c("O2", "O3")) {
      S <- rdl(res, if (sc == "O2") "S2" else "S3")
      k <- if (sc == "O2") kap2 else kap3
      lE <- mask(k * S, inw); names(lE) <- paste0("logE_", sc)
      writeRaster(lE, file.path(d, paste0("logE_", sc, ".tif")), overwrite = TRUE, datatype = "FLT4S")
      v <- values(lE, mat = FALSE)[ic]; e <- exp(v)
      rows[[length(rows) + 1]] <- data.frame(res_m = res, scenario = sc, kappa = round(k, 4),
                                             E_p5 = round(quantile(e, QLO, na.rm = TRUE), 3), E_p95 = round(quantile(e, QHI, na.rm = TRUE), 3),
                                             ratio_p95_p5 = round(quantile(e, QHI, na.rm = TRUE) / quantile(e, QLO, na.rm = TRUE), 3),
                                             E_max_over_min = round(exp(max(v, na.rm = TRUE) - min(v, na.rm = TRUE)), 3), NA_centre = sum(is.na(v)), row.names = NULL)
    }
    O1 <- mask(rdl(res, "S2") * 0, inw); names(O1) <- "logE_O1"
    writeRaster(O1, file.path(d, "logE_O1.tif"), overwrite = TRUE, datatype = "FLT4S")
  }
  tab <- do.call(rbind, rows)
  cat("\nScenario offset range after recalibration:\n"); options(width = 200); print(tab, row.names = FALSE)
  write.csv(tab, file.path(out, "scenario_range_v2.csv"), row.names = FALSE)
  
  sc_def <- data.frame(
    scenario = c("O1", "O2", "O3"),
    definition = c("log E = 0 (no observation offset)",
                   "log E = kappa_O2 * mean(S_access, S_open)",
                   "log E = kappa_O3 * mean(S_access, S_open, S_core)"),
    kappa = c(0, round(kap2, 4), round(kap3, 4)),
    p95_over_p5_target_500m = c(1, CONTRAST_O2, CONTRAST_O3),
    proxies = c("none", "OSM roads (motorway-tertiary), WorldCover 2021 built-up/water/tree",
                "O2 proxies + distance to window edge (pattern-informed, not blind)"),
    note = "Fixed offsets are assumed strengths, not estimates. Free-coefficient variant uses z_S2 / z_S3. Calibrated at 500 m, applied unchanged at 250 and 1000 m.",
    status = "APPROVED (author-delegated 2026-11-06); not FROZEN", stringsAsFactors = FALSE)
  write.csv(sc_def, file.path(REPO_DIR, "docs", "observation_scenarios.csv"), row.names = FALSE)
  cat("\nWrote docs/observation_scenarios.csv\n")
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_14.txt"))
cat("\nDONE. Outputs in:", out, "\n")