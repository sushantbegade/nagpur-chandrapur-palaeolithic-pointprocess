# =============================================================================
# 03_layer_window_check.R
# Stage: 01_data_audit | Gate 4 prep (environment) + Gate 3 prep (window)
# Purpose: layer vs window diagnostics. Uses ONLY window + layers. No site data.
# Outputs: raster_window_check.csv, polygon_coverage_check.csv,
#   polygon_class_areas.csv, rivers_in_window.txt, waterbody_in_window.txt,
#   layer_window_report.txt, sessionInfo_03.txt
# Not executed by assistant. Run top to bottom in RStudio. May take minutes.
# If error mid-run, type sink() in console.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"

source(REPO_CONFIG)
suppressPackageStartupMessages({ library(terra); library(sf) })
sf_use_s2(FALSE)

out <- dir_of("01_data_audit"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
getp <- function(l) file.path(RAW_DIR, man$path[man$layer == l])

sink(file.path(out, "layer_window_report.txt"), split = TRUE)
tryCatch({
  
  win_sf <- st_transform(st_make_valid(st_read(getp("Boundary"), quiet = TRUE)), CRS_PROJ)
  win <- st_union(win_sf)
  win_km2 <- sum(as.numeric(st_area(win))) / 1e6
  bb <- st_bbox(win)
  cat("Window area (km2):", round(win_km2, 1), "\n")
  cat("Window bbox, EPSG", CRS_PROJ, ":", round(bb), "\n")
  bb_ll <- st_bbox(st_transform(win, 4326))
  cat("Window bbox, EPSG 4326 (lon/lat):", round(bb_ll, 4), "\n")
  
  set.seed(MASTER_SEED)
  smp <- st_sf(geometry = st_sample(win, 5000, type = "regular"))
  cat("Regular sample points inside window:", nrow(smp), "\n")
  
  # ---- A. RASTERS: CRS, extent, within-window NA ------------------------------
  cat("\n==== A. RASTERS vs WINDOW ====\n")
  rr <- man[man$type == "raster", ]
  vals_store <- list()
  ras <- do.call(rbind, lapply(seq_len(nrow(rr)), function(i) {
    m <- rr[i, ]
    r <- rast(file.path(RAW_DIR, m$path))
    crs_txt <- crs(r)
    has_crs <- nzchar(crs_txt)
    nm <- if (has_crs) tryCatch(crs(r, describe = TRUE)$name, error = function(e) NA) else NA
    e <- as.vector(ext(r))
    row <- data.frame(layer = m$layer, has_crs = has_crs, crs_name = nm,
                      xmin = e[1], xmax = e[2], ymin = e[3], ymax = e[4], res = res(r)[1],
                      extent_looks_like = if (all(abs(e) <= 360)) "degrees?" else "projected units?",
                      pct_NA_in_window = NA_real_, v_min = NA_real_, v_max = NA_real_,
                      n_unique = NA_integer_, stringsAsFactors = FALSE)
    if (has_crs) tryCatch({
      p <- st_transform(smp, st_crs(crs_txt))
      v <- terra::extract(r, vect(p))[, 2]
      row$pct_NA_in_window <- round(100 * mean(is.na(v)), 2)
      row$v_min <- suppressWarnings(min(v, na.rm = TRUE))
      row$v_max <- suppressWarnings(max(v, na.rm = TRUE))
      row$n_unique <- length(unique(v[!is.na(v)]))
      vals_store[[m$layer]] <<- v
    }, error = function(e) cat("  extract failed:", m$layer, conditionMessage(e), "\n"))
    row
  }))
  write.csv(ras, file.path(out, "raster_window_check.csv"), row.names = FALSE)
  options(width = 200)
  print(ras[, c("layer", "has_crs", "crs_name", "res", "extent_looks_like",
                "pct_NA_in_window", "v_min", "v_max", "n_unique")], right = FALSE)
  
  cat("\nFull extents (for CRS-less soil rasters and bbox comparison):\n")
  print(ras[, c("layer", "xmin", "xmax", "ymin", "ymax")], right = FALSE)
  
  cat("\nSoil raster CRS WKT (first 300 chars):\n")
  for (l in grep("^Soil_", rr$layer, value = TRUE)) {
    r <- rast(file.path(RAW_DIR, rr$path[rr$layer == l]))
    cat(l, ":", if (nzchar(crs(r))) substr(gsub("\\s+", " ", crs(r)), 1, 300) else "NO CRS", "\n")
  }
  cat("\nSoil raster value tables inside window (top 12 values) where CRS exists:\n")
  for (l in intersect(names(vals_store), grep("^Soil_", names(vals_store), value = TRUE))) {
    cat("--", l, "\n"); print(head(sort(table(vals_store[[l]]), decreasing = TRUE), 12))
  }
  if ("DEM" %in% names(vals_store)) {
    cat("\nDEM elevation (m) inside window, quantiles:\n")
    print(round(quantile(vals_store$DEM, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE)))
  }
  
  # ---- B. POLYGON COVERAGE + CLASS AREAS --------------------------------------
  cat("\n==== B. POLYGON LAYERS: COVERAGE OF WINDOW + CLASS AREAS ====\n")
  poly <- list(Geology = "index_", Geomorphology = "legend_sho", Lithology = "lithologic")
  cov_rows <- list(); cls_rows <- list()
  for (l in names(poly)) tryCatch({
    v <- st_transform(st_make_valid(st_read(getp(l), quiet = TRUE)), CRS_PROJ)
    v <- v[, poly[[l]]]; names(v)[1] <- "class"
    u_in <- suppressWarnings(st_intersection(win, st_union(v)))
    cov_km2 <- sum(as.numeric(st_area(u_in))) / 1e6
    v_in <- suppressWarnings(st_intersection(v, win))
    v_in$a_km2 <- as.numeric(st_area(v_in)) / 1e6
    ca <- aggregate(a_km2 ~ class, data = st_drop_geometry(v_in), FUN = sum)
    ca$pct_window <- round(100 * ca$a_km2 / win_km2, 2); ca$a_km2 <- round(ca$a_km2, 1)
    ca <- ca[order(-ca$a_km2), ]; ca$layer <- l; ca$class_field <- poly[[l]]
    cov_rows[[l]] <- data.frame(layer = l, window_km2 = round(win_km2, 1),
                                union_covered_km2 = round(cov_km2, 1), pct_window_covered = round(100 * cov_km2 / win_km2, 2),
                                sum_class_area_km2 = round(sum(ca$a_km2), 1),
                                n_classes_in_window = nrow(ca), n_poly_in_window = nrow(v_in))
    cls_rows[[l]] <- ca
    cat("\n--", l, "(field:", poly[[l]], ") --\n")
    print(cov_rows[[l]], row.names = FALSE); print(head(ca[, c("class", "a_km2", "pct_window")], 25), row.names = FALSE)
  }, error = function(e) cat("FAILED", l, ":", conditionMessage(e), "\n"))
  write.csv(do.call(rbind, cov_rows), file.path(out, "polygon_coverage_check.csv"), row.names = FALSE)
  write.csv(do.call(rbind, cls_rows), file.path(out, "polygon_class_areas.csv"), row.names = FALSE)
  cat("\n(sum_class_area > union_covered means overlapping polygons; investigate if large)\n")
  
  # ---- C. RIVERS --------------------------------------------------------------
  cat("\n==== C. RIVERS IN WINDOW ====\n")
  tryCatch({
    rv <- st_transform(st_read(getp("Rivers"), quiet = TRUE), CRS_PROJ)
    rv <- rv[lengths(st_intersects(rv, win)) > 0, ]
    rc <- suppressWarnings(st_intersection(rv, win))
    rc$len_km <- as.numeric(st_length(rc)) / 1000
    cat("Segments in window:", nrow(rv), " total length (km):", round(sum(rc$len_km)), "\n")
    cat("\nLength (km) by stream order (ordsh):\n"); print(round(tapply(rc$len_km, rc$ordsh, sum), 1))
    cat("\nSegments by stream order:\n"); print(table(rc$ordsh, useNA = "ifany"))
    cat("\nLength (km) by rilcode:\n"); print(round(tapply(rc$len_km, rc$rilcode, sum), 1))
    cat("\nNamed rivers, length (km), top 15:\n")
    print(head(sort(round(tapply(rc$len_km, rc$rivname, sum), 1), decreasing = TRUE), 15))
    cat("\nDrainage density (km/km2):", round(sum(rc$len_km) / win_km2, 3), "\n")
  }, error = function(e) cat("RIVERS FAILED:", conditionMessage(e), "\n"))
  
  # ---- D. WATERBODY -----------------------------------------------------------
  cat("\n==== D. WATERBODY IN WINDOW ====\n")
  tryCatch({
    wb <- st_transform(st_read(getp("Waterbody"), quiet = TRUE), CRS_PROJ)
    inside <- lengths(st_intersects(wb, win)) > 0
    cat("Total features:", nrow(wb), " inside window:", sum(inside), " outside:", sum(!inside), "\n")
    cat("\nInside window by district field:\n"); print(table(wb$district[inside]))
    cat("\nOutside window by district field:\n"); print(table(wb$district[!inside]))
    cat("\nInside window by wbcode:\n"); print(table(wb$wbcode[inside]))
    cat("\nArea (ha) inside window, quantiles:\n"); print(round(quantile(wb$area_ha[inside], c(0, .25, .5, .75, .95, 1), na.rm = TRUE), 2))
  }, error = function(e) cat("WATERBODY FAILED:", conditionMessage(e), "\n"))
  
}, finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_03.txt"))
cat("\nDONE. Outputs in:", out, "\n")