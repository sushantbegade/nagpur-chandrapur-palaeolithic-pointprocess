# =============================================================================
# 02_site_audit.R   (v2: CSV read before sink; tapply bug fixed)
# Stage: 01_data_audit | Gate 1 (data integrity) + Gate 2 (spatial precision)
# Purpose: automatic audit of the 197-record master inventory. NOTHING is dropped,
#          edited or re-coordinated. Flags only. Inclusion decisions are manual.
# Outputs (01_data_audit/outputs/): site_audit_report.txt, master_schema_draft.csv,
#   near_pairs_le1000m.csv, reconciliation_auto.csv, sessionInfo_02.txt
# Not executed by assistant. Run top to bottom in RStudio.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"

while (sink.number() > 0) sink()          # close any stale sink from earlier failed runs
source(REPO_CONFIG)
suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)

out <- dir_of("01_data_audit"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
getp <- function(l) file.path(RAW_DIR, man$path[man$layer == l])

# PROVISIONAL sigma (m), used ONLY for the Gate 2 feasibility diagnostic below.
# NOT a protocol decision. Author must define and justify per class at freeze.
PROVISIONAL_SIGMA_M <- c(FIELD = 10, AUTHOR_REPORTED = 10, GIS_MAP = 250, UNKNOWN = 250)

# ---- 1. READ CSV (before any sink) ------------------------------------------
f <- getp("SiteTable_csv")
cat("Site file:", f, "| exists:", file.exists(f), "\n")
stopifnot(length(f) == 1, file.exists(f))
lines <- readLines(f, warn = FALSE)
txt   <- paste(lines, collapse = "\n")     # rejoin so quoted embedded newlines survive
sites   <- read.csv(text = txt, stringsAsFactors = FALSE, check.names = FALSE)
raw_chr <- read.csv(text = txt, stringsAsFactors = FALSE, check.names = FALSE,
                    colClasses = "character")
stopifnot(ncol(sites) == 18)
cn <- c("site_id","site_name","latitude","longitude","district","state","river",
        "basin","coord_precision_raw","lithic_assemblage","site_type_raw",
        "raw_material","surface_context","cultural_attribution","chronology_raw",
        "pub_year_raw","primary_reference","data_source_raw")
names(sites) <- cn; names(raw_chr) <- cn
n <- nrow(sites)

sink(file.path(out, "site_audit_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== 1. INTEGRITY ====\n")
  cat("Rows:", n, "\n")
  cat("site_id unique:", !anyDuplicated(sites$site_id),
      "| consecutive 1..n:", all(sites$site_id == seq_len(n)), "\n")
  cat("Lat NA:", sum(is.na(sites$latitude)), "| Lon NA:", sum(is.na(sites$longitude)), "\n")
  cat("Lat range:", range(sites$latitude, na.rm = TRUE), "\n")
  cat("Lon range:", range(sites$longitude, na.rm = TRUE), "\n")
  
  # ---- 2. COORDINATE PRECISION CLASS + DECIMALS --------------------------------
  pr <- sites$coord_precision_raw
  prec <- ifelse(grepl("Field Verified", pr, ignore.case = TRUE), "FIELD",
                 ifelse(grepl("Reported By", pr, ignore.case = TRUE), "AUTHOR_REPORTED",
                        ifelse(grepl("GIS", pr, ignore.case = TRUE), "GIS_MAP", "UNKNOWN")))
  ndec <- function(s) { s <- trimws(s); ifelse(grepl("\\.", s), nchar(sub("^[^.]*\\.", "", s)), 0L) }
  dec_lat <- ndec(raw_chr$latitude); dec_lon <- ndec(raw_chr$longitude)
  
  cat("\n==== 2. COORDINATE PRECISION CLASSES ====\n")
  print(table(prec, useNA = "ifany"))
  cat("\nDecimal places of latitude by class (as stored in CSV):\n")
  print(table(prec, dec_lat))
  cat("\nDecimal places of longitude by class:\n")
  print(table(prec, dec_lon))
  cat("(1e-4 deg ~ 11 m; 1e-3 deg ~ 111 m; 1e-5 ~ 1.1 m)\n")
  
  # ---- 3. SPATIAL: window, district, edge distance ----------------------------
  pts_ll <- st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)
  pts    <- st_transform(pts_ll, CRS_PROJ)
  bnd    <- st_transform(st_make_valid(st_read(getp("Boundary"), quiet = TRUE)), CRS_PROJ)
  win    <- st_union(bnd)
  win_area_m2 <- sum(as.numeric(st_area(win)))
  
  in_window <- lengths(st_intersects(pts, win)) > 0
  d_edge    <- as.numeric(drop(st_distance(pts, st_boundary(win))))
  d_edge_signed_km <- round(ifelse(in_window, d_edge, -d_edge) / 1000, 2)
  
  idx <- sapply(st_intersects(pts, bnd), function(i) if (length(i)) i[1] else NA_integer_)
  poly_district <- toupper(bnd$District[idx])
  district_mismatch <- toupper(trimws(sites$district)) != poly_district
  district_mismatch[is.na(district_mismatch)] <- FALSE
  
  cat("\n==== 3. WINDOW / DISTRICT ====\n")
  cat("Window area (km2):", round(win_area_m2 / 1e6, 1), "\n")
  cat("Sites inside window:", sum(in_window), " outside:", sum(!in_window), "\n")
  if (any(!in_window)) { cat("Outside-window sites:\n")
    print(data.frame(site_id = sites$site_id[!in_window], name = sites$site_name[!in_window],
                     dist_outside_km = -d_edge_signed_km[!in_window])) }
  cat("Reported district != polygon district:", sum(district_mismatch), "\n")
  if (any(district_mismatch)) print(data.frame(site_id = sites$site_id[district_mismatch],
                                               name = sites$site_name[district_mismatch], reported = sites$district[district_mismatch],
                                               polygon = poly_district[district_mismatch]))
  cat("\nDistance to window edge (km), inside sites:\n")
  print(round(quantile(d_edge[in_window] / 1000, c(0, .05, .1, .25, .5, .75, 1)), 2))
  cat("Sites within 5 km of edge:", sum(in_window & d_edge < 5000),
      "| within 10 km:", sum(in_window & d_edge < 10000), "\n")
  
  # ---- 4. DUPLICATES / NEAR-DUPLICATES / NN DISTANCES -------------------------
  nm_norm <- tolower(gsub("[^a-z0-9]", "", iconv(trimws(sites$site_name), to = "ASCII//TRANSLIT")))
  name_dup <- duplicated(nm_norm) | duplicated(nm_norm, fromLast = TRUE)
  exact_xy_dup <- duplicated(raw_chr[, c("latitude", "longitude")]) |
    duplicated(raw_chr[, c("latitude", "longitude")], fromLast = TRUE)
  
  D <- matrix(as.numeric(st_distance(pts)), nrow = n); diag(D) <- Inf
  nn <- apply(D, 1, min)
  near100 <- rowSums(D <= 100) > 0; near500 <- rowSums(D <= 500) > 0; near1000 <- rowSums(D <= 1000) > 0
  
  cat("\n==== 4. DUPLICATES / NEAREST NEIGHBOUR ====\n")
  cat("Duplicate (normalised) names:", sum(name_dup), "records\n")
  if (any(name_dup)) print(sites[name_dup, c("site_id", "site_name", "river", "primary_reference")])
  cat("Identical coordinate strings:", sum(exact_xy_dup), "records\n")
  cat("Sites with a neighbour <=100 m:", sum(near100), "| <=500 m:", sum(near500),
      "| <=1000 m:", sum(near1000), "\n")
  cat("\nNearest-neighbour distance (m), all sites:\n")
  print(round(quantile(nn, c(0, .05, .1, .25, .5, .75, .9, 1))))
  cat("\nNN distance (m) by coordinate class (median / min):\n")
  print(round(t(sapply(split(nn, prec), function(x) c(median = median(x), min = min(x))))))
  
  ij <- which(D <= 1000 & upper.tri(D), arr.ind = TRUE)
  pairs <- if (nrow(ij)) data.frame(
    id_a = sites$site_id[ij[, 1]], name_a = sites$site_name[ij[, 1]],
    id_b = sites$site_id[ij[, 2]], name_b = sites$site_name[ij[, 2]],
    dist_m = round(D[ij]), prec_a = prec[ij[, 1]], prec_b = prec[ij[, 2]],
    ref_a = sites$primary_reference[ij[, 1]], ref_b = sites$primary_reference[ij[, 2]],
    stringsAsFactors = FALSE) else data.frame()
  if (nrow(pairs)) pairs <- pairs[order(pairs$dist_m), ]
  write.csv(pairs, file.path(out, "near_pairs_le1000m.csv"), row.names = FALSE)
  cat("\nPairs <=1000 m:", nrow(pairs), "(saved). First 25:\n"); print(head(pairs, 25))
  
  # crude CSR reference, descriptive only (no edge correction, no inference)
  exp_nn <- 0.5 / sqrt(n / win_area_m2)
  cat("\nCrude Clark-Evans ratio (uncorrected, descriptive only): mean NN =",
      round(mean(nn)), "m; CSR expectation =", round(exp_nn), "m; R =", round(mean(nn) / exp_nn, 2), "\n")
  
  # ---- 5. GATE 2 FEASIBILITY (provisional sigma) ------------------------------
  sig <- PROVISIONAL_SIGMA_M[prec]
  cat("\n==== 5. GATE 2 DIAGNOSTIC (PROVISIONAL sigma, author must set) ====\n")
  print(PROVISIONAL_SIGMA_M)
  cat("Sites whose NN distance < 2*sigma_i:", sum(nn < 2 * sig), "of", n, "\n")
  cat("Median sigma / median NN:", round(median(sig) / median(nn), 4), "\n")
  
  # ---- 6. SITE TYPE, SOURCE, CHRONOLOGY, REFERENCES ---------------------------
  st_raw <- tolower(trimws(sites$site_type_raw))
  site_type <- ifelse(grepl("riverbed", st_raw), "riverbed",
                      ifelse(grepl("lake", st_raw), "lake_shore",
                             ifelse(grepl("open", st_raw), "open_air",
                                    ifelse(grepl("factory|trench", st_raw), "excavated_workshop", "not_reported"))))
  ds <- trimws(sites$data_source_raw)
  blank_ds <- is.na(ds) | ds == ""
  has_exc  <- grepl("excavation", ds, ignore.case = TRUE)
  
  chr <- tolower(sites$chronology_raw)
  LP <- grepl("lower", chr); MP <- grepl("middle", chr); UP <- grepl("upper", chr)
  n_phase <- LP + MP + UP
  phase <- apply(cbind(ifelse(LP, "LP", ""), ifelse(MP, "MP", ""), ifelse(UP, "UP", "")), 1,
                 function(r) paste(r[r != ""], collapse = "+"))
  
  yrs <- regmatches(sites$pub_year_raw, gregexpr("[0-9]{4}", sites$pub_year_raw))
  yr_min <- sapply(yrs, function(x) if (length(x)) min(as.integer(x)) else NA_integer_)
  yr_max <- sapply(yrs, function(x) if (length(x)) max(as.integer(x)) else NA_integer_)
  n_refs <- sapply(strsplit(sites$primary_reference, ";"), length)
  
  cat("\n==== 6. SITE TYPE / SOURCE / CHRONOLOGY / REFERENCES ====\n")
  cat("Site type (normalised):\n"); print(table(site_type))
  cat("\nSite type x data source (blank -> '(blank)'):\n")
  print(table(site_type, ifelse(blank_ds, "(blank)", ds)))
  cat("\nPhase combinations:\n"); print(table(phase, useNA = "ifany"))
  cat("\nAny-phase membership: LP =", sum(LP), " MP =", sum(MP), " UP =", sum(UP),
      "| unparsed:", sum(n_phase == 0), "\n")
  cat("\nPublication year (min per record): "); print(summary(yr_min))
  cat("Records with >1 reference:", sum(n_refs > 1), "\n")
  cat("\nDistinct primary-reference strings:", length(unique(sites$primary_reference)), "\nTop 25:\n")
  print(head(sort(table(sites$primary_reference), decreasing = TRUE), 25))
  cat("\nSurface-context values (top 15):\n")
  print(head(sort(table(sites$surface_context), decreasing = TRUE), 15))
  cat("\nRaw-material values (top 15):\n")
  print(head(sort(table(sites$raw_material), decreasing = TRUE), 15))
  
  # ---- 7. FLAGS, SCHEMA DRAFT, RECONCILIATION ---------------------------------
  flags <- cbind(
    OUTSIDE_WINDOW = !in_window, DISTRICT_MISMATCH = district_mismatch,
    NAME_DUP = name_dup, EXACT_XY_DUP = exact_xy_dup, NEAR_LE500M = near500,
    BLANK_DATA_SOURCE = blank_ds, GIS_MAP_COORD = prec == "GIS_MAP",
    NO_PHASE_PARSED = n_phase == 0, SITE_TYPE_NOT_REPORTED = site_type == "not_reported",
    EXCAVATION = has_exc)
  dq <- apply(flags, 1, function(r) paste(names(r)[which(r)], collapse = ";"))
  
  schema <- data.frame(
    site_id = sites$site_id, canonical_name = trimws(sites$site_name),
    source_reference = sites$primary_reference,
    latitude = sites$latitude, longitude = sites$longitude,
    coordinate_precision = sites$coord_precision_raw, coordinate_source = prec,
    coordinate_confidence = NA_character_,
    archaeological_context = paste(site_type, "|", sites$surface_context),
    Palaeolithic_phase = phase, chronological_confidence = NA_character_,
    site_unit_status = "PENDING", survey_source = ifelse(blank_ds, NA, ds),
    survey_project = NA_character_, inclusion_status = "PENDING_REVIEW",
    exclusion_reason = NA_character_, data_quality_flag = dq,
    in_window = in_window, dist_edge_km_signed = d_edge_signed_km,
    nn_dist_m = round(nn), site_type_norm = site_type,
    pub_year_min = yr_min, pub_year_max = yr_max, n_references = n_refs,
    stringsAsFactors = FALSE)
  write.csv(schema, file.path(out, "master_schema_draft.csv"), row.names = FALSE)
  
  recon <- data.frame(
    stage = c("Master records", "Exact duplicate rows", "Identical coordinate strings (records)",
              "Duplicate names (records)", "Outside window", "District mismatch",
              "Neighbour <=500 m (records, nested-locality candidates)",
              "GIS-map coordinates (records)", "Blank data source", "No phase parsed",
              "Site type not reported", "NOTE"),
    n = c(n, sum(duplicated(sites)), sum(exact_xy_dup), sum(name_dup), sum(!in_window),
          sum(district_mismatch), sum(near500), sum(prec == "GIS_MAP"), sum(blank_ds),
          sum(n_phase == 0), sum(site_type == "not_reported"), NA),
    status = c(rep("FLAG_ONLY - no record removed", 11), "Final NA requires manual review"))
  write.csv(recon, file.path(out, "reconciliation_auto.csv"), row.names = FALSE)
  cat("\n==== 7. AUTOMATIC RECONCILIATION (flags only) ====\n"); print(recon, right = FALSE)
  cat("\nRecords with >=1 flag:", sum(dq != ""), "| flag counts:\n"); print(colSums(flags))
  
}, finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_02.txt"))
cat("\nDONE. Outputs in:", out, "\n")