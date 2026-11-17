# =====================================================================
# 37_inspect_observation_layers.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : locate and describe the observation-offset layers (O2, O3, log E)
#           needed for recovery cells N1 and N5. Read-only. Opens only files
#           listed in docs/freeze_v1_manifest.csv (no folder scanning).
# Outputs : outputs/37_report.txt, 37_obs_layers.csv, 37_sessionInfo.txt
# Notes   : Manifest 'path' may be relative to REPO_DIR, to RAW_DIR, or absolute;
#           all three are tried. Rows are selected by a name pattern (INFERENCE:
#           layer names are not known to me); the full group table is printed.
#           Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)
suppressPackageStartupMessages(library(terra))

OUT_DIR <- dir_of("10_simulation_recovery", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
PAT <- "obs|offset|kappa|scenario|logE|log_E|worldcover|road|core|edge|proxy|z_S|z_elev_p1"
MAX_CELLS <- 5e7

# ---- open all inputs BEFORE sink ----
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
resolve <- function(p) {
  cand <- c(p, file.path(REPO_DIR, p), file.path(RAW_DIR, p))
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else NA_character_
}
MAN$resolved <- vapply(MAN$path, resolve, character(1))
sel <- grepl(PAT, MAN$path, ignore.case = TRUE) | grepl(PAT, MAN$group, ignore.case = TRUE)
S <- MAN[sel, , drop = FALSE]
rinfo <- function(p) {
  if (is.na(p) || !grepl("\\.(tif|tiff)$", p, ignore.case = TRUE)) return(NULL)
  tryCatch({
    r <- terra::rast(p)
    e <- as.vector(terra::ext(r))
    big <- terra::ncell(r) > MAX_CELLS
    st <- if (big) c(NA_real_, NA_real_, NA_real_, NA_real_) else {
      g <- terra::global(r[[1]], c("min", "max", "mean"), na.rm = TRUE)
      c(g$min, g$max, g$mean, terra::global(is.na(r[[1]]), "mean")[1, 1])
    }
    data.frame(path = p, nlyr = terra::nlyr(r), ncol = ncol(r), nrow = nrow(r), xres = terra::res(r)[1],
               xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4], vmin = st[1], vmax = st[2], vmean = st[3],
               na_frac = st[4], stringsAsFactors = FALSE)
  }, error = function(e) data.frame(path = p, nlyr = NA, ncol = NA, nrow = NA, xres = NA, xmin = NA, xmax = NA,
                                    ymin = NA, ymax = NA, vmin = NA, vmax = NA, vmean = NA, na_frac = NA,
                                    stringsAsFactors = FALSE))
}
RI <- do.call(rbind, lapply(S$resolved, rinfo))
CSVS <- S$resolved[!is.na(S$resolved) & grepl("\\.csv$", S$resolved, ignore.case = TRUE)]
CSVD <- lapply(CSVS, function(p) tryCatch(read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL))
names(CSVD) <- CSVS

REPORT <- file.path(OUT_DIR, "37_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("37_inspect_observation_layers | ", format(Sys.time()), "\n", sep = "")
  cat(sprintf("Manifest rows %d | resolvable %d | selected by pattern %d\n", nrow(MAN), sum(!is.na(MAN$resolved)), nrow(S)))
  cat("\n== Manifest groups ==\n"); print(table(MAN$group, useNA = "ifany"))
  cat("\n== Example manifest paths (first 3) ==\n"); print(head(MAN$path, 3))
  cat("\n== Selected rows (group, freeze_status, resolved, path) ==\n")
  print(data.frame(group = S$group, status = S$freeze_status, found = !is.na(S$resolved), path = S$path),
        row.names = FALSE, right = FALSE)
  cat("\n== Raster info for selected .tif rows ==\n")
  if (!is.null(RI) && nrow(RI)) print(RI, digits = 5, row.names = FALSE) else cat("no .tif rows resolved\n")
  if (!is.null(RI) && nrow(RI)) write.csv(RI, file.path(OUT_DIR, "37_obs_layers.csv"), row.names = FALSE)
  cat("\n== Selected .csv tables (full print if <= 25 rows) ==\n")
  for (p in names(CSVD)) {
    d <- CSVD[[p]]
    cat("\n--", basename(p), "--\n")
    if (is.null(d)) cat("could not read\n") else print(if (nrow(d) <= 25) d else head(d, 10), row.names = FALSE)
  }
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "37_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")