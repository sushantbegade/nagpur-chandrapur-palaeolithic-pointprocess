# =============================================================================
# 04_site_review_sheet.R
# Stage: 01_data_audit | Gate 1 manual-review support
# Purpose: (a) verify name-duplicate flags, (b) verify site-type mapping,
#   (c) list pairs within 2 km with name similarity, (d) descriptive edge check,
#   (e) write review sheet for AUTHOR decisions.
# NOTHING is dropped or edited. Decisions are the author's, recorded in the sheet.
# Outputs (01_data_audit/outputs/): site_review_report.txt, site_review_sheet.csv,
#   pairs_le2000m_name_check.csv, sessionInfo_04.txt
# Requires master_schema_draft.csv from 02_site_audit.R (v2).
# Not executed by assistant. Run top to bottom in RStudio.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"

while (sink.number() > 0) sink()
source(REPO_CONFIG)
suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)

out <- dir_of("01_data_audit"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
getp <- function(l) file.path(RAW_DIR, man$path[man$layer == l])

# ---- READ (before sink) ------------------------------------------------------
f <- getp("SiteTable_csv")
txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
sites <- read.csv(text = txt, stringsAsFactors = FALSE, check.names = FALSE)
cn <- c("site_id","site_name","latitude","longitude","district","state","river",
        "basin","coord_precision_raw","lithic_assemblage","site_type_raw",
        "raw_material","surface_context","cultural_attribution","chronology_raw",
        "pub_year_raw","primary_reference","data_source_raw")
names(sites) <- cn
n <- nrow(sites)
sch <- read.csv(file.path(out, "master_schema_draft.csv"), stringsAsFactors = FALSE)
stopifnot(nrow(sch) == n, all(sch$site_id == sites$site_id))

pts <- st_transform(st_as_sf(sites, coords = c("longitude", "latitude"),
                             crs = 4326, remove = FALSE), CRS_PROJ)
D <- matrix(as.numeric(st_distance(pts)), nrow = n); diag(D) <- Inf
bnd <- st_union(st_transform(st_make_valid(st_read(getp("Boundary"), quiet = TRUE)), CRS_PROJ))

sink(file.path(out, "site_review_report.txt"), split = TRUE)
tryCatch({
  
  # ---- A. NAME-DUPLICATE VERIFICATION -----------------------------------------
  cat("==== A. NAME-DUPLICATE VERIFICATION ====\n")
  nm <- trimws(sites$site_name)
  non_ascii <- grepl("[^ -~]", nm, perl = TRUE)
  cat("Names with non-ASCII characters:", sum(non_ascii), "\n")
  if (any(non_ascii)) print(data.frame(site_id = sites$site_id[non_ascii], name = nm[non_ascii]))
  tr <- iconv(nm, to = "ASCII//TRANSLIT", sub = "")
  cat("Names where iconv gave NA:", sum(is.na(tr)), "\n")
  tr[is.na(tr)] <- nm[is.na(tr)]
  nm_norm <- tolower(gsub("[^[:alnum:]]", "", tr))
  name_dup <- duplicated(nm_norm) | duplicated(nm_norm, fromLast = TRUE)
  cat("Duplicate names after NA-safe normalisation:", sum(name_dup), "records\n")
  if (any(name_dup)) print(sites[name_dup, c("site_id", "site_name", "river", "primary_reference")])
  chk <- c(2, 35, 80, 84, 101, 187)
  cat("\nBytes of names for IDs flagged in script 02 (hidden characters show here):\n")
  for (i in chk) cat(sprintf("ID %3d [%s]: %s\n", i, nm[i], paste(charToRaw(nm[i]), collapse = " ")))
  cat("\nKargaon pair (80, 84) distance (m):", round(D[80, 84]), "\n")
  
  # ---- B. SITE-TYPE MAPPING CHECK ---------------------------------------------
  cat("\n==== B. SITE TYPE: RAW vs MAPPED ====\n")
  print(table(raw = sprintf("[%s]", sites$site_type_raw), mapped = sch$site_type_norm))
  cat("\nRecords with raw site type blank / 'Not Reported' / trench:\n")
  sel <- grepl("^\\s*$|not reported|trench|factory", sites$site_type_raw, ignore.case = TRUE)
  print(sites[sel, c("site_id", "site_name", "site_type_raw", "surface_context", "data_source_raw")])
  
  # ---- C. MEGALITH / REDEPOSITED CONTEXT --------------------------------------
  cat("\n==== C. MEGALITH / CAIRN CONTEXT ====\n")
  meg <- grepl("megalith|cairn", paste(sites$surface_context, sites$site_name,
                                       sites$cultural_attribution, sites$lithic_assemblage),
               ignore.case = TRUE)
  cat("Records:", sum(meg), "\n")
  if (any(meg)) print(sites[meg, c("site_id", "site_name", "surface_context", "cultural_attribution",
                                   "chronology_raw", "primary_reference")])
  
  # ---- D. PAIRS <= 2 km WITH NAME SIMILARITY ----------------------------------
  cat("\n==== D. PAIRS <= 2000 m (names, references) ====\n")
  ij <- which(D <= 2000 & upper.tri(D), arr.ind = TRUE)
  pr <- data.frame(
    id_a = sites$site_id[ij[, 1]], name_a = nm[ij[, 1]],
    id_b = sites$site_id[ij[, 2]], name_b = nm[ij[, 2]],
    dist_m = round(D[ij]),
    prec_a = sch$coordinate_source[ij[, 1]], prec_b = sch$coordinate_source[ij[, 2]],
    ref_a = sites$primary_reference[ij[, 1]], ref_b = sites$primary_reference[ij[, 2]],
    stringsAsFactors = FALSE)
  pr$name_dist <- round(mapply(function(x, y)
    adist(tolower(x), tolower(y))[1, 1] / max(nchar(x), nchar(y)),
    pr$name_a, pr$name_b, USE.NAMES = FALSE), 2)
  pr$same_ref <- pr$ref_a == pr$ref_b
  pr <- pr[order(pr$dist_m), ]
  write.csv(pr, file.path(out, "pairs_le2000m_name_check.csv"), row.names = FALSE,
            fileEncoding = "UTF-8")
  cat("Pairs <=2000 m:", nrow(pr), "(name_dist 0 = identical, 1 = fully different)\n")
  options(width = 220); print(pr, row.names = FALSE)
  
  # ---- E. EDGE PROXIMITY (descriptive, no inference) --------------------------
  cat("\n==== E. EDGE PROXIMITY: observed vs area-expected ====\n")
  A <- sum(as.numeric(st_area(bnd)))
  d_in_m <- sch$dist_edge_km_signed * 1000
  for (d in c(2000, 5000, 10000)) {
    inner <- st_buffer(bnd, -d)
    frac <- 1 - sum(as.numeric(st_area(inner))) / A
    cat(sprintf("<%2d km of edge: window area share %.1f%% -> expected %.0f sites | observed %d\n",
                d / 1000, 100 * frac, n * frac, sum(d_in_m < d)))
  }
  cat("Descriptive only. Sites are inhomogeneous and not independent.\n")
  
  # ---- F. REVIEW SHEET ---------------------------------------------------------
  sim_partner <- sapply(sites$site_id, function(id) {
    k <- which((pr$id_a == id | pr$id_b == id) & pr$name_dist <= 0.35 & pr$dist_m <= 2000)
    if (!length(k)) return("")
    paste(ifelse(pr$id_a[k] == id, pr$id_b[k], pr$id_a[k]), collapse = ";")
  })
  partners <- sapply(seq_len(n), function(i) {
    j <- which(D[i, ] <= 1000)
    if (length(j)) paste(sprintf("%d:%dm", sites$site_id[j], round(D[i, j])), collapse = "; ") else ""
  })
  extra <- cbind(MEGALITH_CONTEXT = meg, NAME_DUP_CONFIRMED = name_dup,
                 NAME_SIMILAR_NEAR = sim_partner != "")
  extra_txt <- apply(extra, 1, function(r) paste(names(r)[which(r)], collapse = ";"))
  base_flags <- sch$data_quality_flag; base_flags[is.na(base_flags)] <- ""
  all_flags <- ifelse(base_flags == "", extra_txt,
                      ifelse(extra_txt == "", base_flags, paste(base_flags, extra_txt, sep = ";")))
  strip <- function(s, drop) vapply(strsplit(s, ";", fixed = TRUE),
                                    function(x) paste(setdiff(x, drop), collapse = ";"), "")
  needs_review <- strip(all_flags, c("GIS_MAP_COORD", "EXCAVATION")) != "" | partners != ""
  
  review <- data.frame(
    site_id = sites$site_id, site_name = nm,
    latitude = sites$latitude, longitude = sites$longitude,
    coordinate_source = sch$coordinate_source,
    site_type_raw = sites$site_type_raw, site_type_norm = sch$site_type_norm,
    surface_context = sites$surface_context, chronology_raw = sites$chronology_raw,
    source_reference = sites$primary_reference, data_source_raw = sites$data_source_raw,
    auto_flags = all_flags, partners_le1000m = partners, similar_name_partners_le2km = sim_partner,
    needs_review = needs_review,
    decision_inclusion = "PENDING",   # author: INCLUDE / EXCLUDE
    decision_site_unit = "PENDING",   # author: site / locality / sub-locality / findspot / duplicate
    duplicate_of_id = NA, exclusion_reason = NA, decision_note = NA,
    stringsAsFactors = FALSE)
  write.csv(review, file.path(out, "site_review_sheet.csv"), row.names = FALSE,
            fileEncoding = "UTF-8")
  cat("\n==== F. REVIEW SHEET ====\n")
  cat("Rows:", nrow(review), "| needs_review (ignoring GIS_MAP/EXCAVATION flags alone):",
      sum(needs_review), "\n")
  cat("Saved: site_review_sheet.csv\n")
  
}, finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_04.txt"))
cat("\nDONE. Outputs in:", out, "\n")