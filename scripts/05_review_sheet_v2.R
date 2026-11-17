# =============================================================================
# 05_review_sheet_v2.R
# Stage: 01_data_audit | Gate 1 manual-review support (v2)
# Fixes vs scripts 02/04: (1) false NAME_DUP removed, (2) site-type mapping v2
#   (factory_site separate from excavated_in_situ), (3) candidate same-locality
#   pairs by token-set name match, (4) reference-string checks, (5) source
#   concentration table. NOTHING dropped/edited. Decisions = author.
# Outputs (01_data_audit/outputs/): site_review_report_v2.txt,
#   site_review_sheet_v2.csv, candidate_pairs.csv, source_concentration.csv,
#   sessionInfo_05.txt
# Requires master_schema_draft.csv from 02_site_audit.R (v2).
# Not executed by assistant. Run top to bottom in RStudio.
#
# DECISION VOCABULARY (author fills in v2 sheet):
#   decision_inclusion : INCLUDE / EXCLUDE
#   decision_site_unit : site / locality / sub-locality / findspot / duplicate
#   duplicate_of_id    : site_id of the record kept (only if duplicate)
#   exclusion_reason   : text, mandatory if EXCLUDE
#   evidence_checked   : Y / N  (did you consult the source paper?)
# Decide on provenance/archaeology only. Never on effect on a statistic.
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
nm  <- trimws(sites$site_name)
ref <- trimws(sites$primary_reference)

sink(file.path(out, "site_review_report_v2.txt"), split = TRUE)
tryCatch({
  
  # ---- A. SITE TYPE v2 --------------------------------------------------------
  cat("==== A. SITE-TYPE MAPPING v2 ====\n")
  st_raw <- tolower(trimws(sites$site_type_raw))
  site_type2 <- ifelse(grepl("riverbed", st_raw), "riverbed",
                       ifelse(grepl("lake", st_raw), "lake_shore",
                              ifelse(grepl("insitu|trench", st_raw), "excavated_in_situ",
                                     ifelse(grepl("factory", st_raw), "factory_site",
                                            ifelse(grepl("^open", st_raw), "open_air", "not_reported")))))
  print(table(raw = sprintf("[%s]", sites$site_type_raw), mapped = site_type2))
  cat("\nMapped totals:\n"); print(table(site_type2))
  
  # ---- B. CANDIDATE SAME-LOCALITY PAIRS ---------------------------------------
  cat("\n==== B. CANDIDATE SAME-LOCALITY PAIRS ====\n")
  nm_norm <- tolower(gsub("[^[:alnum:]]", "", nm))
  tok_key <- vapply(strsplit(tolower(gsub("[^[:alnum:]]+", " ", nm)), "\\s+"),
                    function(x) paste(sort(unique(x[x != ""])), collapse = " "), "")
  name_dup_exact <- duplicated(nm_norm) | duplicated(nm_norm, fromLast = TRUE)
  tok_dup <- duplicated(tok_key) | duplicated(tok_key, fromLast = TRUE)
  cat("Exact normalised-name duplicates (records):", sum(name_dup_exact), "\n")
  cat("Same-name-token duplicates, any order (records):", sum(tok_dup), "\n")
  
  ND <- adist(tolower(nm)) / outer(nchar(nm), nchar(nm), pmax)
  c1 <- outer(tok_key, tok_key, "==")
  c2 <- (ND <= 0.35) & (D <= 2000)
  c3 <- D <= 500
  ij <- which(upper.tri(D) & (c1 | c2 | c3), arr.ind = TRUE)
  
  if (nrow(ij)) {
    reason <- apply(cbind(c1[ij], c2[ij], c3[ij]), 1, function(r)
      paste(c("SAME_NAME_TOKENS", "SIMILAR_NAME_LE2KM", "LE500M")[r], collapse = "+"))
    cand <- data.frame(
      id_a = sites$site_id[ij[, 1]], name_a = nm[ij[, 1]],
      id_b = sites$site_id[ij[, 2]], name_b = nm[ij[, 2]],
      dist_m = round(D[ij]), reason = reason,
      name_dist = round(ND[ij], 2),
      prec_a = sch$coordinate_source[ij[, 1]], prec_b = sch$coordinate_source[ij[, 2]],
      river_a = sites$river[ij[, 1]], river_b = sites$river[ij[, 2]],
      ref_a = ref[ij[, 1]], ref_b = ref[ij[, 2]], stringsAsFactors = FALSE)
    cand <- cand[order(cand$dist_m), ]
  } else cand <- data.frame()
  write.csv(cand, file.path(out, "candidate_pairs.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  cat("Candidate pairs:", nrow(cand), "\n")
  options(width = 230); if (nrow(cand)) print(cand, row.names = FALSE)
  
  cand_partner <- vapply(seq_len(n), function(i) {
    k <- which(ij[, 1] == i | ij[, 2] == i)
    if (!length(k)) "" else paste(ifelse(ij[k, 1] == i, sites$site_id[ij[k, 2]],
                                         sites$site_id[ij[k, 1]]), collapse = ";")
  }, "")
  
  # ---- C. SECONDARY CONTEXT ----------------------------------------------------
  cat("\n==== C. MEGALITH / CAIRN CONTEXT ====\n")
  meg <- grepl("megalith|cairn", paste(sites$surface_context, sites$site_name,
                                       sites$cultural_attribution, sites$lithic_assemblage),
               ignore.case = TRUE)
  cat("Records:", sum(meg), "\n")
  if (any(meg)) print(sites[meg, c("site_id", "site_name", "surface_context",
                                   "chronology_raw", "primary_reference")])
  
  # ---- D. REFERENCE CHECKS -----------------------------------------------------
  cat("\n==== D. REFERENCE CHECKS ====\n")
  yrs_ref <- regmatches(ref, gregexpr("[0-9]{4}", ref))
  yrs_pub <- regmatches(sites$pub_year_raw, gregexpr("[0-9]{4}", sites$pub_year_raw))
  yr_mismatch <- mapply(function(a, b) !setequal(a, b), yrs_ref, yrs_pub)
  cat("Records where years in reference string != publication-year field:", sum(yr_mismatch), "\n")
  if (any(yr_mismatch)) print(data.frame(site_id = sites$site_id[yr_mismatch],
                                         name = nm[yr_mismatch], reference = ref[yr_mismatch],
                                         pub_year_field = sites$pub_year_raw[yr_mismatch]), row.names = FALSE)
  
  pat <- "[A-Z][A-Za-z'\\-]+(?=(?: et al\\.)?(?: and [A-Z][A-Za-z]+)?,? [0-9]{4})"
  auth <- unlist(regmatches(ref, gregexpr(pat, ref, perl = TRUE)))
  at <- sort(table(auth), decreasing = TRUE)
  cat("\nAuthor tokens extracted:", length(at), "(count of mentions)\n")
  an <- names(at); AD <- adist(tolower(an))
  w <- which(AD <= 2 & upper.tri(AD), arr.ind = TRUE)
  cat("\nSimilar author spellings (possible typos), edit distance <= 2:\n")
  if (nrow(w)) print(data.frame(a = an[w[, 1]], n_a = as.integer(at[w[, 1]]),
                                b = an[w[, 2]], n_b = as.integer(at[w[, 2]]), dist = AD[w]), row.names = FALSE) else cat("none\n")
  
  cat("\nAll author tokens with count (singletons often typos):\n"); print(at)
  
  # ---- E. SOURCE CONCENTRATION ------------------------------------------------
  cat("\n==== E. SOURCE CONCENTRATION (lead author of reference string) ====\n")
  lead <- sub("^\\s*([^ ,;]+).*$", "\\1", ref)
  has2k <- rowSums(D <= 2000) > 0
  sc <- data.frame(lead_author = names(table(lead)), n_sites = as.integer(table(lead)),
                   stringsAsFactors = FALSE)
  sc$n_with_neighbour_le2km <- as.integer(tapply(has2k, lead, sum)[sc$lead_author])
  sc <- sc[order(-sc$n_sites), ]
  sc$pct_of_197 <- round(100 * sc$n_sites / n, 1)
  sc$cum_pct <- round(cumsum(sc$pct_of_197), 1)
  write.csv(sc, file.path(out, "source_concentration.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  print(head(sc, 15), row.names = FALSE)
  cat("Lead authors:", nrow(sc), "| top-3 share:", sc$cum_pct[3], "% | top-5 share:", sc$cum_pct[5], "%\n")
  
  # ---- F. REVIEW SHEET v2 ------------------------------------------------------
  base <- sch$data_quality_flag; base[is.na(base)] <- ""
  drop_tok <- function(s, drop) vapply(strsplit(s, ";", fixed = TRUE),
                                       function(x) paste(setdiff(x, drop), collapse = ";"), "")
  base2 <- drop_tok(base, "NAME_DUP")                      # remove false flags from script 02
  extra <- cbind(NAME_DUP_CONFIRMED = name_dup_exact,
                 SAME_NAME_TOKENS = tok_dup,
                 CANDIDATE_PAIR = cand_partner != "",
                 SECONDARY_CONTEXT = meg,
                 REF_YEAR_MISMATCH = yr_mismatch)
  extra_txt <- apply(extra, 1, function(r) paste(names(r)[which(r)], collapse = ";"))
  all_flags <- ifelse(base2 == "", extra_txt,
                      ifelse(extra_txt == "", base2, paste(base2, extra_txt, sep = ";")))
  partners <- vapply(seq_len(n), function(i) {
    j <- which(D[i, ] <= 1000)
    if (length(j)) paste(sprintf("%d:%dm", sites$site_id[j], round(D[i, j])), collapse = "; ") else ""
  }, "")
  needs_review <- drop_tok(all_flags, c("GIS_MAP_COORD", "EXCAVATION")) != "" | partners != ""
  
  review <- data.frame(
    site_id = sites$site_id, site_name = nm,
    latitude = sites$latitude, longitude = sites$longitude,
    coordinate_source = sch$coordinate_source, river = sites$river,
    site_type_raw = sites$site_type_raw, site_type_norm_v2 = site_type2,
    surface_context = sites$surface_context, chronology_raw = sites$chronology_raw,
    source_reference = ref, pub_year_raw = sites$pub_year_raw,
    data_source_raw = sites$data_source_raw,
    auto_flags = all_flags, partners_le1000m = partners,
    candidate_same_locality_ids = cand_partner, needs_review = needs_review,
    decision_inclusion = "PENDING", decision_site_unit = "PENDING",
    duplicate_of_id = NA, exclusion_reason = NA, evidence_checked = NA, decision_note = NA,
    stringsAsFactors = FALSE)
  write.csv(review, file.path(out, "site_review_sheet_v2.csv"), row.names = FALSE,
            fileEncoding = "UTF-8")
  
  cat("\n==== F. REVIEW SHEET v2 ====\n")
  cat("Rows:", nrow(review), "| needs_review:", sum(needs_review), "\n")
  cat("Flag counts (v2):\n")
  ft <- unlist(strsplit(all_flags[all_flags != ""], ";", fixed = TRUE))
  print(sort(table(ft), decreasing = TRUE))
  cat("\nneeds_review IDs:\n"); print(sites$site_id[needs_review])
  cat("Saved: site_review_sheet_v2.csv (use THIS, not v1)\n")
  
}, finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_05.txt"))
cat("\nDONE. Outputs in:", out, "\n")