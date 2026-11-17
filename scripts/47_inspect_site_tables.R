# =====================================================================
# 47_inspect_site_tables.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW; read-only)
# Purpose : locate the per-site information needed for the robustness stage:
#           coordinate class (FIELD / AUTHOR_REPORTED / GIS_MAP), chronology
#           (LP / MP / UP any-phase sets), site type (fluvial / non-fluvial, Dataset C),
#           redeposited flags (IDs 48, 49), and the contents of ppp_NA_window.rds.
#           Opens only data_derived/ppp_NA_window.rds and the 'derived_table' rows of
#           docs/freeze_v1_manifest.csv. No folder scanning.
# Outputs : outputs/47_report.txt, 47_sessionInfo.txt
# Notes   : Prints column names, classes, unique-value counts, tabulations of columns with
#           up to 15 unique values, and the first 3 rows of each table. No analysis.
#           Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)
suppressPackageStartupMessages(library(spatstat.geom))

OUT_DIR <- dir_of("11_generative_gof", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
f_ppp <- file.path(REPO_DIR, "data_derived", "ppp_NA_window.rds"); stopifnot(file.exists(f_ppp))
MAX_UNI <- 15

# ---- open all inputs BEFORE sink ----
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
OBJ <- readRDS(f_ppp)
TB  <- MAN[MAN$group == "derived_table", , drop = FALSE]
TB$exists <- file.exists(TB$path)
CSV <- TB$path[TB$exists & grepl("\\.csv$", TB$path, ignore.case = TRUE)]
TD  <- lapply(CSV, function(p) tryCatch(read.csv(p, stringsAsFactors = FALSE), error = function(e) NULL))
names(TD) <- CSV
LAY <- MAN[MAN$group == "derived_layer", , drop = FALSE]

REPORT <- file.path(OUT_DIR, "47_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("47_inspect_site_tables | ", format(Sys.time()), "\n", sep = "")
  
  cat("\n== A. ppp_NA_window.rds ==\n")
  cat("class:", paste(class(OBJ), collapse = ","), "| length", length(OBJ), "\n")
  if (is.list(OBJ)) for (nm in names(OBJ)) {
    x <- OBJ[[nm]]
    cat(sprintf("  $%-14s class %-22s %s\n", nm, paste(class(x), collapse = ","),
                if (is.ppp(x)) sprintf("npoints %d | window %s | marks: %s", npoints(x), Window(x)$type,
                                       if (is.null(marks(x))) "none" else paste(class(marks(x)), collapse = ","))
                else if (is.data.frame(x)) sprintf("data.frame %d x %d: %s", nrow(x), ncol(x), paste(names(x), collapse = ","))
                else if (is.atomic(x)) sprintf("length %d", length(x)) else ""))
  }
  X <- OBJ$X_km
  if (!is.null(X) && !is.null(marks(X))) { cat("marks of X_km:\n"); print(head(as.data.frame(marks(X)), 3)) }
  
  cat("\n== B. Manifest: derived_table rows (", nrow(TB), ") ==\n", sep = "")
  print(data.frame(exists = TB$exists, path = TB$path), row.names = FALSE, right = FALSE)
  
  cat("\n== C. Derived layers by directory ==\n")
  print(table(dirname(LAY$path)))
  
  cat("\n== D. Per-table structure ==\n")
  for (p in names(TD)) {
    d <- TD[[p]]
    cat("\n-----", basename(p), "-----\n")
    if (is.null(d)) { cat("could not read\n"); next }
    cat("rows", nrow(d), "| cols", ncol(d), "\n")
    info <- data.frame(column = names(d), class = vapply(d, function(x) class(x)[1], ""),
                       n_unique = vapply(d, function(x) length(unique(x)), 1L),
                       n_na = vapply(d, function(x) sum(is.na(x) | x == ""), 1L), stringsAsFactors = FALSE)
    print(info, row.names = FALSE)
    for (cn in names(d)) {
      u <- unique(d[[cn]])
      if (length(u) <= MAX_UNI && length(u) >= 1 && nrow(d) > 20) { cat("  [", cn, "]\n", sep = ""); print(table(d[[cn]], useNA = "ifany")) }
    }
    cat("first 3 rows:\n"); print(head(d, 3), row.names = FALSE)
  }
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "47_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")