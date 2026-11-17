# =====================================================================
# 61_gis_class_provenance.R
# Stage   : 13_sensitivity
# Folder  : 13_sensitivity/scripts/   (NEW; read-only, descriptive)
# Purpose : describe what the GIS_MAP coordinate class consists of, to support the provenance
#           question raised by script 60: (A) class by lead author, river, site type, data source,
#           chronology; (B) spatial summary by class and nearest-neighbour class mixing;
#           (C) the GIS-GIS pairs within 5 km (names, distance, rivers, authors, references);
#           (D) shared-author / shared-reference shares of close pairs by class pair.
# Outputs : outputs/61_report.txt, 61_gis_gis_pairs.csv, 61_sessionInfo.txt
# Notes   : No inference, no test, no change to any result. Opens only sites_NA.csv (manifest).
#           Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

OUT_DIR <- dir_of("13_sensitivity", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
S <- read.csv(f_sites, stringsAsFactors = FALSE)
stopifnot(nrow(S) == 197)

REPORT <- file.path(OUT_DIR, "61_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  options(width = 220)
  cat("61_gis_class_provenance (DESCRIPTIVE) | ", format(Sys.time()), "\n", sep = "")
  cls <- S$coordinate_source; gis <- cls == "GIS_MAP"
  
  cat("\n== A. Class composition ==\n")
  cat("coordinate class x lead author (ordered by GIS_MAP count):\n")
  ta <- table(S$lead_author, cls); ta <- ta[order(-ta[, "GIS_MAP"], -rowSums(ta)), , drop = FALSE]; print(ta)
  cat("\nGIS_MAP sites per river (all classes shown):\n"); tr <- table(S$river, cls); print(tr[order(-tr[, "GIS_MAP"], -rowSums(tr)), , drop = FALSE])
  cat("\ncoordinate class x site type:\n"); print(table(S$site_type, cls))
  cat("\ncoordinate class x data source:\n"); print(table(S$data_source, cls))
  cat("\ncoordinate class x phase combination:\n"); print(table(S$phase_combo, cls))
  cat("\ndistinct source references per class:", paste(tapply(S$source_reference, cls, function(x) length(unique(x))), collapse = " / "), "(", paste(names(table(cls)), collapse = " / "), ")\n")
  
  cat("\n== B. Space ==\n")
  x <- S$x_m / 1000; y <- S$y_m / 1000
  SP <- do.call(rbind, lapply(sort(unique(cls)), function(k) data.frame(class = k, n = sum(cls == k), cx_km = mean(x[cls == k]), cy_km = mean(y[cls == k]),
                                                                        sd_x = sd(x[cls == k]), sd_y = sd(y[cls == k]), mean_dist_to_class_centroid = mean(sqrt((x[cls == k] - mean(x[cls == k]))^2 + (y[cls == k] - mean(y[cls == k]))^2)))))
  print(SP, digits = 4, row.names = FALSE)
  D <- as.matrix(dist(cbind(x, y))); diag(D) <- Inf
  nn <- apply(D, 1, which.min)
  cat(sprintf("\nnearest neighbour in the same class: GIS_MAP sites %.3f (random expectation %.3f) | FIELD %.3f (%.3f) | AUTHOR_REPORTED %.3f (%.3f)\n",
              mean(cls[nn][gis] == "GIS_MAP"), 59 / 196, mean(cls[nn][cls == "FIELD"] == "FIELD"), 63 / 196,
              mean(cls[nn][cls == "AUTHOR_REPORTED"] == "AUTHOR_REPORTED"), 72 / 196))
  cat("nearest-neighbour distance (km) by class: ", paste(sprintf("%s median %.2f", sort(unique(cls)), tapply(apply(D, 1, min), cls, median)), collapse = " | "), "\n")
  
  cat("\n== C. GIS_MAP - GIS_MAP pairs within 5 km ==\n")
  ij <- which(D < 5 & upper.tri(D), arr.ind = TRUE); ij <- ij[gis[ij[, 1]] & gis[ij[, 2]], , drop = FALSE]
  P <- data.frame(id_a = S$site_id[ij[, 1]], id_b = S$site_id[ij[, 2]], name_a = S$canonical_name[ij[, 1]], name_b = S$canonical_name[ij[, 2]], dist_km = D[ij],
                  river_a = S$river[ij[, 1]], river_b = S$river[ij[, 2]], author_a = S$lead_author[ij[, 1]], author_b = S$lead_author[ij[, 2]],
                  same_author = S$lead_author[ij[, 1]] == S$lead_author[ij[, 2]], same_reference = S$source_reference[ij[, 1]] == S$source_reference[ij[, 2]],
                  type_a = S$site_type[ij[, 1]], type_b = S$site_type[ij[, 2]], stringsAsFactors = FALSE)
  P <- P[order(P$dist_km), ]
  cat(sprintf("pairs: %d | same lead author %.2f | same source reference %.2f | same river %.2f\n", nrow(P), mean(P$same_author), mean(P$same_reference), mean(P$river_a == P$river_b)))
  print(P, digits = 3, row.names = FALSE); write.csv(P, file.path(OUT_DIR, "61_gis_gis_pairs.csv"), row.names = FALSE)
  
  cat("\n== D. Shared author / reference among close pairs (within 5 km), by class pair ==\n")
  ij2 <- which(D < 5 & upper.tri(D), arr.ind = TRUE)
  cc <- paste(pmin(cls[ij2[, 1]], cls[ij2[, 2]]), pmax(cls[ij2[, 1]], cls[ij2[, 2]]), sep = " + ")
  sa <- S$lead_author[ij2[, 1]] == S$lead_author[ij2[, 2]]; sr <- S$source_reference[ij2[, 1]] == S$source_reference[ij2[, 2]]
  print(data.frame(class_pair = names(table(cc)), pairs = as.integer(table(cc)), share_same_author = as.numeric(tapply(sa, cc, mean)), share_same_reference = as.numeric(tapply(sr, cc, mean))), digits = 3, row.names = FALSE)
  cat(sprintf("\nall sites, any distance: share of site pairs with the same lead author %.3f | same source reference %.3f\n",
              mean(outer(S$lead_author, S$lead_author, "==")[upper.tri(D)]), mean(outer(S$source_reference, S$source_reference, "==")[upper.tri(D)])))
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "61_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")