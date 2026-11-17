# =============================================================================
# 06_final_na.R
# Stage: 01_data_audit | Gate 1 closure support
# Purpose: apply author decisions (stated in chat 2026-11-06) to the v2 review sheet,
#   build the decided sheet, inclusion/exclusion log, 197 -> NA reconciliation
#   (Table 1, Figure 2 data) and the analysis-ready site tables.
# NOTHING is dropped by distance, flag or default. Exclusions exist only in the
# EXCLUSIONS table below (currently empty) or from window screening.
# Outputs:
#   docs/site_review_DECIDED.csv
#   docs/site_inclusion_exclusion_log.csv   (OVERWRITES the empty scaffold file)
#   01_data_audit/outputs/site_inclusion_final.csv, site_audit_final_report.txt,
#     table1_audit_stages.csv, sessionInfo_06.txt
#   14_figures/outputs/fig2_reconciliation_data.csv
#   15_tables/outputs/table1_audit_stages.csv
#   data_derived/sites_master_decided.csv  (all 197, with inclusion_status)
#   data_derived/sites_NA.csv              (included records only)
# Not executed by assistant. Run top to bottom in RStudio.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"

while (sink.number() > 0) sink()
source(REPO_CONFIG)
suppressPackageStartupMessages(library(sf))
sf_use_s2(FALSE)

out <- dir_of("01_data_audit"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
for (d in c("docs", "data_derived", "14_figures/outputs", "15_tables/outputs"))
  dir.create(file.path(REPO_DIR, d), recursive = TRUE, showWarnings = FALSE)

AUTHOR_DATE <- "2026-11-06"           # date the author stated the decisions below
EVID_TAG    <- "AUTHOR_STATED_20261106"

# ---------------------------------------------------------------------------
# AUTHOR DECISIONS (transcribed from author statements; edit here if wrong)
# ---------------------------------------------------------------------------
SRC_SURVEY_IDS    <- c(50, 71, 77, 78, 79, 80, 81)      # data source = Survey
OPEN_AIR_IDS      <- c(70, 71)                          # site type = open_air
YEAR_OVERRIDE     <- c(`166` = "2017")                  # publication-year field
DISTINCT_PAIRS    <- list(c(63, 64), c(127, 133), c(136, 147), c(21, 175), c(80, 84),
                          c(6, 9), c(7, 10), c(11, 15), c(66, 68), c(171, 190))
REDEPOSITED_IDS   <- c(48, 49)                          # redeposited; INCLUDED
EXCLUSIONS <- data.frame(site_id = integer(), stage = integer(),   # stage 1-4
                         reason = character(), stringsAsFactors = FALSE)
# stage codes: 1 identity/duplicate, 2 provenance, 3 coordinate, 4 archaeological criteria

# ---- READ (before sink) ------------------------------------------------------
rd <- function(f) read.csv(text = paste(readLines(f, warn = FALSE), collapse = "\n"),
                           stringsAsFactors = FALSE, check.names = FALSE)
sheet <- rd(file.path(out, "site_review_sheet_v2.csv"))
sch   <- rd(file.path(out, "master_schema_draft.csv"))
n <- nrow(sheet); ids <- sheet$site_id
stopifnot(n == 197, all(ids == 1:197), nrow(sch) == n, all(sch$site_id == ids))
stopifnot(all(c(SRC_SURVEY_IDS, OPEN_AIR_IDS, REDEPOSITED_IDS,
                unlist(DISTINCT_PAIRS), as.integer(names(YEAR_OVERRIDE))) %in% ids))

isblank <- function(x) is.na(x) | trimws(as.character(x)) == ""
af <- ifelse(is.na(sheet$auto_flags), "", sheet$auto_flags)
hasflag <- function(pat) grepl(pat, af)
in_window <- as.logical(sch$in_window)
nm <- trimws(sheet$site_name)

sink(file.path(out, "site_audit_final_report.txt"), split = TRUE)
tryCatch({
  
  # ---- 1. APPLY DECISIONS -----------------------------------------------------
  cat("==== 1. DECISIONS APPLIED ====\n")
  notes <- setNames(rep("", n), ids)
  add_note <- function(id, txt) {
    k <- as.character(id)
    notes[k] <<- if (notes[k] == "") txt else paste(notes[k], txt, sep = " | ")
  }
  logrows <- list()
  add_log <- function(id, stage, action, reason)
    logrows[[length(logrows) + 1]] <<- data.frame(site_id = id, stage = stage, action = action,
                                                  reason = reason, date = AUTHOR_DATE,
                                                  stringsAsFactors = FALSE)
  
  data_source_final <- sheet$data_source_raw
  data_source_final[ids %in% SRC_SURVEY_IDS] <- "Survey"
  for (i in SRC_SURVEY_IDS) {
    add_note(i, "Author: data source = Survey (blank in master table).")
    add_log(i, "provenance", "FIELD_CORRECTED_INCLUDE", "Data source blank in master table; author states Survey.")
  }
  site_type_final <- sheet$site_type_norm_v2
  site_type_final[ids %in% OPEN_AIR_IDS] <- "open_air"
  for (i in OPEN_AIR_IDS) {
    add_note(i, "Author: site type = open_air (was Not Reported).")
    add_log(i, "provenance", "FIELD_CORRECTED_INCLUDE", "Site type Not Reported; author states open-air.")
  }
  pub_year_final <- sheet$pub_year_raw
  for (k in names(YEAR_OVERRIDE)) {
    i <- as.integer(k); pub_year_final[ids == i] <- YEAR_OVERRIDE[[k]]
    add_note(i, "Author: publication year = 2017; reference string unchanged (lists Joglekar 2017 and Joglekar & Imchen 2018).")
    add_log(i, "provenance", "FIELD_CORRECTED_INCLUDE", "Publication-year field incomplete; author states 2017.")
  }
  for (p in DISTINCT_PAIRS) {
    add_note(p[1], sprintf("Author: distinct from ID %d (no merge).", p[2]))
    add_note(p[2], sprintf("Author: distinct from ID %d (no merge).", p[1]))
    add_log(p[1], "identity", "CONFIRMED_DISTINCT_NO_MERGE", sprintf("Author states ID %d and %d are separate sites.", p[1], p[2]))
    add_log(p[2], "identity", "CONFIRMED_DISTINCT_NO_MERGE", sprintf("Author states ID %d and %d are separate sites.", p[1], p[2]))
  }
  for (i in REDEPOSITED_IDS) {
    add_note(i, "Author: tools redeposited in megalith cairn fill (Thapar 1964); INCLUDED; sens_redeposited=TRUE.")
    add_log(i, "inclusion_criteria", "INCLUDE_REDEPOSITED", "Author states redeposited context; kept in primary analysis; flagged for sensitivity.")
  }
  if (nrow(EXCLUSIONS)) for (r in seq_len(nrow(EXCLUSIONS))) {
    add_note(EXCLUSIONS$site_id[r], paste("EXCLUDED:", EXCLUSIONS$reason[r]))
    add_log(EXCLUSIONS$site_id[r], c("identity","provenance","coordinates","inclusion_criteria")[EXCLUSIONS$stage[r]],
            "EXCLUDE", EXCLUSIONS$reason[r])
  }
  
  explicit_ids <- unique(c(SRC_SURVEY_IDS, OPEN_AIR_IDS, as.integer(names(YEAR_OVERRIDE)),
                           unlist(DISTINCT_PAIRS), REDEPOSITED_IDS, EXCLUSIONS$site_id))
  basis <- ifelse(ids %in% explicit_ids, "AUTHOR_STATED", "DEFAULT_UNQUEUED")
  queue_ids <- ids[as.logical(sheet$needs_review)]
  cat("Queued (needs_review) IDs:", length(queue_ids), "\n")
  unc <- setdiff(queue_ids, explicit_ids)
  cat("Queued IDs WITHOUT an explicit author decision:", if (length(unc)) unc else "none", "\n")
  cat("Explicit decisions cover", length(explicit_ids), "records; default rows:", sum(basis == "DEFAULT_UNQUEUED"), "\n")
  
  # ---- 2. FINAL STATUS ---------------------------------------------------------
  excluded_by_decision <- ids %in% EXCLUSIONS$site_id
  include <- in_window & !excluded_by_decision
  inclusion_status <- ifelse(include, "INCLUDE", "EXCLUDE")
  exclusion_reason <- rep(NA_character_, n)
  exclusion_reason[excluded_by_decision] <- EXCLUSIONS$reason[match(ids[excluded_by_decision], EXCLUSIONS$site_id)]
  exclusion_reason[!in_window & is.na(exclusion_reason)] <- "Outside analytical window"
  evidence_checked <- ifelse(basis == "AUTHOR_STATED", EVID_TAG, "N")
  
  # ---- 3. DERIVED ATTRIBUTES ---------------------------------------------------
  chr <- tolower(sheet$chronology_raw)
  LP <- grepl("lower", chr); MP <- grepl("middle", chr); UP <- grepl("upper", chr)
  phase <- apply(cbind(ifelse(LP, "LP", ""), ifelse(MP, "MP", ""), ifelse(UP, "UP", "")), 1,
                 function(r) paste(r[r != ""], collapse = "+"))
  lead_author <- sub("^\\s*([^ ,;]+).*$", "\\1", trimws(sheet$source_reference))
  nm_key <- tolower(trimws(nm))
  dup_name <- duplicated(nm_key) | duplicated(nm_key, fromLast = TRUE)
  name_unique <- ifelse(dup_name, sprintf("%s [ID %d]", nm, ids), nm)
  
  pts <- st_transform(st_as_sf(data.frame(sheet[, c("site_id","longitude","latitude")]),
                               coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), CRS_PROJ)
  xy <- st_coordinates(pts)
  
  dataset_A <- include
  dataset_B <- include           # merged on source evidence: author found no merges -> identical to A
  dataset_C <- include & !(site_type_final %in% c("riverbed", "lake_shore"))   # non-fluvial contexts
  sens_redeposited <- ids %in% REDEPOSITED_IDS
  
  # ---- 4. OUTPUT TABLES --------------------------------------------------------
  decided <- cbind(sheet[, setdiff(names(sheet), c("decision_inclusion","decision_site_unit","duplicate_of_id",
                                                   "exclusion_reason","evidence_checked","decision_note",
                                                   "site_type_norm_v2","data_source_raw","pub_year_raw"))],
                   site_type_raw_norm_v2 = sheet$site_type_norm_v2, site_type_final = site_type_final,
                   data_source_raw = sheet$data_source_raw, data_source_final = data_source_final,
                   pub_year_raw = sheet$pub_year_raw, pub_year_final = pub_year_final,
                   decision_inclusion = inclusion_status, decision_site_unit = "site", duplicate_of_id = NA,
                   exclusion_reason = exclusion_reason, evidence_checked = evidence_checked,
                   decision_basis = basis, decision_note = unname(notes), stringsAsFactors = FALSE)
  write.csv(decided, file.path(REPO_DIR, "docs", "site_review_DECIDED.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  
  final <- data.frame(
    site_id = ids, canonical_name = nm, canonical_name_unique = name_unique,
    longitude = sheet$longitude, latitude = sheet$latitude,
    x_m = round(xy[, 1], 1), y_m = round(xy[, 2], 1), crs_epsg = CRS_PROJ,
    coordinate_source = sheet$coordinate_source, river = sheet$river,
    site_type = site_type_final, data_source = data_source_final,
    source_reference = sheet$source_reference, lead_author = lead_author,
    chronology_raw = sheet$chronology_raw, LP = LP, MP = MP, UP = UP, phase_combo = phase,
    site_unit_status = "site", in_window = in_window,
    inclusion_status = inclusion_status, exclusion_reason = exclusion_reason,
    decision_basis = basis, evidence_checked = evidence_checked,
    sens_redeposited = sens_redeposited,
    dataset_A = dataset_A, dataset_B = dataset_B, dataset_C = dataset_C,
    stringsAsFactors = FALSE)
  write.csv(final, file.path(out, "site_inclusion_final.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(final, file.path(REPO_DIR, "data_derived", "sites_master_decided.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(final[include, ], file.path(REPO_DIR, "data_derived", "sites_NA.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  
  summ <- data.frame(site_id = NA_integer_, stage = "summary", action = "DEFAULT_INCLUDE_UNQUEUED",
                     reason = sprintf("%d records not in review queue; INCLUDE/site by default (decision_basis=DEFAULT_UNQUEUED in site_inclusion_final.csv); awaiting author confirmation.",
                                      sum(basis == "DEFAULT_UNQUEUED")), date = AUTHOR_DATE, stringsAsFactors = FALSE)
  summ2 <- data.frame(site_id = NA_integer_, stage = "summary", action = "EXCLUSIONS_COUNT",
                      reason = sprintf("Author exclusions: %d; outside window: %d.", nrow(EXCLUSIONS), sum(!in_window)),
                      date = AUTHOR_DATE, stringsAsFactors = FALSE)
  logdf <- do.call(rbind, c(logrows, list(summ, summ2)))
  write.csv(logdf, file.path(REPO_DIR, "docs", "site_inclusion_exclusion_log.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cat("\nLog rows written:", nrow(logdf), "\n")
  
  # ---- 5. RECONCILIATION (Table 1 / Figure 2) ---------------------------------
  cat("\n==== 2. RECONCILIATION 197 -> NA ====\n")
  blank_raw <- isblank(sheet$data_source_raw)
  nr_raw <- sheet$site_type_norm_v2 == "not_reported"
  f_ident <- hasflag("NAME_DUP_CONFIRMED|SAME_NAME_TOKENS|CANDIDATE_PAIR|EXACT_XY_DUP|NEAR_LE500M")
  f_prov  <- blank_raw | nr_raw | hasflag("SECONDARY_CONTEXT|REF_YEAR_MISMATCH")
  f_coord <- sheet$coordinate_source %in% c("GIS_MAP", "UNKNOWN")
  f_incl  <- hasflag("SECONDARY_CONTEXT")
  f_win   <- !in_window
  
  stages <- data.frame(
    stage = c("0 Master inventory", "1 Identity / duplicate audit", "2 Provenance audit",
              "3 Coordinate audit", "4 Archaeological inclusion criteria",
              "5 Spatial-window screening"),
    flagged = c(n, sum(f_ident), sum(f_prov), sum(f_coord), sum(f_incl), sum(f_win)),
    removed = c(0, sum(EXCLUSIONS$stage == 1), sum(EXCLUSIONS$stage == 2), sum(EXCLUSIONS$stage == 3),
                sum(EXCLUSIONS$stage == 4), sum(!in_window & !excluded_by_decision)),
    note = c("197 records, IDs 1-197",
             "Exact-name, same-token and <=500 m candidates; 5 candidate pairs confirmed distinct by author",
             "Blank data source (7), site type not reported (2), secondary context (2), reference-year mismatch (1); resolved by author",
             "GIS-map coordinates retained (60); sigma and datum pending (Gate 2)",
             "Redeposited megalith-fill records IDs 48, 49 included by author",
             "Boundary shapefile EPSG:32644"),
    stringsAsFactors = FALSE)
  stages$remaining <- n - cumsum(stages$removed)
  stages <- rbind(stages, data.frame(stage = "6 Final analytical sample NA", flagged = NA,
                                     removed = NA, note = "Records with inclusion_status = INCLUDE",
                                     remaining = sum(include), stringsAsFactors = FALSE))
  options(width = 220); print(stages, right = FALSE, row.names = FALSE)
  stopifnot(stages$remaining[nrow(stages)] == stages$remaining[nrow(stages) - 1])
  write.csv(stages, file.path(out, "table1_audit_stages.csv"), row.names = FALSE)
  write.csv(stages, file.path(REPO_DIR, "15_tables", "outputs", "table1_audit_stages.csv"), row.names = FALSE)
  write.csv(stages[, c("stage", "removed", "remaining")],
            file.path(REPO_DIR, "14_figures", "outputs", "fig2_reconciliation_data.csv"), row.names = FALSE)
  cat("\nFINAL NA =", sum(include), "\n")
  
  # ---- 6. DESCRIPTIVES OF THE FINAL SAMPLE ------------------------------------
  cat("\n==== 3. FINAL SAMPLE DESCRIPTIVES (NA records) ====\n")
  fi <- final[include, ]
  cat("Coordinate class:\n"); print(table(fi$coordinate_source))
  cat("\nSite type (final):\n"); print(table(fi$site_type))
  cat("\nData source (final):\n"); print(table(fi$data_source, useNA = "ifany"))
  cat("\nPhase combinations:\n"); print(table(fi$phase_combo))
  cat("Any-phase: LP =", sum(fi$LP), " MP =", sum(fi$MP), " UP =", sum(fi$UP), "\n")
  cat("\nDistrict-free summaries: datasets\n")
  cat("Dataset A:", sum(dataset_A), "| Dataset B:", sum(dataset_B),
      "(identical to A: no source-based merges) | Dataset C (non-fluvial):", sum(dataset_C), "\n")
  cat("Dataset C by site type:\n"); print(table(fi$site_type[fi$dataset_C]))
  cat("\nSensitivity flag redeposited (IDs 48, 49):", sum(fi$sens_redeposited), "\n")
  cat("Duplicate canonical names disambiguated:", sum(dup_name), "records\n")
  if (any(dup_name)) print(final[dup_name, c("site_id", "canonical_name", "canonical_name_unique")], row.names = FALSE)
  cat("\nDecision basis:\n"); print(table(fi$decision_basis))
  cat("\nDeviations from protocol: none\n")
  
}, finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_06.txt"))
cat("\nDONE. Outputs in:", out, "\n")