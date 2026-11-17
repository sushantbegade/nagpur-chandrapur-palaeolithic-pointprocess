# =============================================================================
# 07_window_and_ppp.R   (v2: cross-check column guarded; sections isolated)
# Stage: 02_spatial_window | Gate 3 (analytical window) support
# Purpose: build the analytical window (spatstat owin) from the boundary shapefile,
#   build the NA=197 point pattern (metres and km) with marks table, window/edge
#   diagnostics, pair counts by distance, DEM metadata dump.
# Uses ONLY geometry. No covariates, no models, no site-response relationships.
# Outputs:
#   02_spatial_window/outputs/: window_report.txt, window_summary.csv,
#     eroded_area_table.csv, pair_counts.csv, dem_gdal_info.txt, sessionInfo_07.txt
#   data_derived/ppp_NA_window.rds  (list: W_m, X_m, X_km, XC_m, marks, notes)
# Not executed by assistant. Run top to bottom in RStudio. Takes a few minutes
# (edge erosion of a ~19,000-vertex polygon).
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"

while (sink.number() > 0) sink()
source(REPO_CONFIG)

need <- c("spatstat.geom", "sf", "terra")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: install.packages(c(",
                       paste(sprintf('"%s"', miss), collapse = ", "), "))")
suppressPackageStartupMessages({ library(sf); library(terra); library(spatstat.geom) })
sf_use_s2(FALSE)

out <- dir_of("02_spatial_window"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(REPO_DIR, "data_derived"), showWarnings = FALSE)
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
getp <- function(l) file.path(RAW_DIR, man$path[man$layer == l])
rd <- function(f) read.csv(text = paste(readLines(f, warn = FALSE), collapse = "\n"),
                           stringsAsFactors = FALSE, check.names = FALSE)

# ---- READ (before sink) ------------------------------------------------------
bnd <- st_transform(st_make_valid(st_read(getp("Boundary"), quiet = TRUE)), CRS_PROJ)
sites_f <- file.path(REPO_DIR, "data_derived", "sites_NA.csv")
stopifnot(file.exists(sites_f))
S <- rd(sites_f)
sch_f <- file.path(dir_of("01_data_audit"), "master_schema_draft.csv")
sch <- if (file.exists(sch_f)) rd(sch_f) else NULL
dem_info <- tryCatch(terra::describe(getp("DEM")),
                     error = function(e) paste("describe failed:", conditionMessage(e)))
writeLines(dem_info, file.path(out, "dem_gdal_info.txt"))

# ---- helpers -----------------------------------------------------------------
ring_area <- function(x, y) { k <- c(2:length(x), 1); 0.5 * sum(x * y[k] - x[k] * y) }
manual_owin <- function(g) {
  mp <- st_cast(st_union(g), "MULTIPOLYGON")
  parts <- st_cast(st_sfc(mp, crs = st_crs(g)), "POLYGON")
  rings <- list()
  for (i in seq_along(parts)) {
    co <- st_coordinates(parts[i])
    for (r in unique(co[, "L1"])) {
      m <- co[co[, "L1"] == r, c("X", "Y"), drop = FALSE]
      m <- m[-nrow(m), , drop = FALSE]
      a <- ring_area(m[, 1], m[, 2])
      if ((a > 0) != (r == 1)) m <- m[nrow(m):1, , drop = FALSE]
      rings[[length(rings) + 1]] <- list(x = m[, 1], y = m[, 2])
    }
  }
  owin(poly = rings)
}
safe <- function(label, expr) {
  tryCatch(expr, error = function(e)
    cat(sprintf("\n[SECTION FAILED: %s] %s\n", label, conditionMessage(e))))
}
buffer_area <- function(g, r_km) {
  e <- tryCatch(st_buffer(g, -r_km * 1000), error = function(e) NULL)
  if (is.null(e) || length(e) == 0 || all(st_is_empty(e))) 0
  else sum(as.numeric(st_area(e))) / 1e6
}

sink(file.path(out, "window_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== 0. VERSIONS ====\n")
  for (p in c("spatstat.geom", "sf", "terra")) cat(sprintf("%-14s %s\n", p, as.character(packageVersion(p))))
  cat("CRS_PROJ: EPSG", CRS_PROJ, "\n")
  
  # ---- 1. WINDOW GEOMETRY (sf) ------------------------------------------------
  cat("\n==== 1. WINDOW GEOMETRY (sf) ====\n")
  win_sfc <- st_union(st_geometry(bnd))
  parts <- st_cast(st_sfc(st_cast(win_sfc, "MULTIPOLYGON"), crs = CRS_PROJ), "POLYGON")
  nrings <- vapply(seq_along(parts), function(i) length(unique(st_coordinates(parts[i])[, "L1"])), 1L)
  nvert  <- nrow(st_coordinates(win_sfc))
  A_sf   <- sum(as.numeric(st_area(win_sfc))) / 1e6
  perim  <- sum(as.numeric(st_length(st_boundary(win_sfc)))) / 1000
  cat("Dissolved window: parts =", length(parts), "| holes =", sum(nrings - 1),
      "| vertices =", nvert, "\n")
  cat("Area (km2):", round(A_sf, 1), "| perimeter (km):", round(perim, 1),
      "| compactness 4*pi*A/P^2:", round(4 * pi * A_sf / perim^2, 3), "\n")
  
  # ---- 2. spatstat OWIN --------------------------------------------------------
  cat("\n==== 2. spatstat owin ====\n")
  W <- tryCatch(as.owin(win_sfc), error = function(e) { cat("sf-route FAILED:", conditionMessage(e), "\n"); NULL })
  route <- "sf as.owin"
  if (is.null(W) || abs(area.owin(W) / 1e6 - A_sf) / A_sf > 1e-3) {
    W <- manual_owin(win_sfc); route <- "manual polygon builder"
  }
  unitname(W) <- c("metre", "metres")
  cat("Construction route:", route, "| type:", W$type, "| area (km2):", round(area.owin(W) / 1e6, 1), "\n")
  
  # ---- 3. POINT PATTERN --------------------------------------------------------
  cat("\n==== 3. POINT PATTERN (NA records) ====\n")
  n <- nrow(S)
  stopifnot(n == 197, all(S$inclusion_status == "INCLUDE"), !anyDuplicated(S$site_id))
  pts <- st_transform(st_as_sf(S, coords = c("longitude", "latitude"), crs = 4326, remove = FALSE), CRS_PROJ)
  xy <- st_coordinates(pts)
  rej <- 0
  X <- withCallingHandlers(
    ppp(x = xy[, 1], y = xy[, 2], window = W, check = TRUE),
    warning = function(w) { rej <<- rej + 1; cat("ppp warning:", conditionMessage(w), "\n"); invokeRestart("muffleWarning") })
  cat("Points in ppp:", npoints(X), "of", n, "| warnings:", rej,
      "| duplicated:", sum(duplicated(X)), "\n")
  stopifnot(npoints(X) == n)
  unitname(X) <- c("metre", "metres")
  marks_tab <- S[, c("site_id", "coordinate_source", "site_type", "LP", "MP", "UP", "phase_combo",
                     "dataset_C", "sens_redeposited", "lead_author", "decision_basis")]
  X_km  <- rescale(X, 1000, c("kilometre", "kilometres"))
  keepC <- as.logical(S$dataset_C)
  XC    <- X[keepC]
  cat("Dataset C ppp points:", npoints(XC), "(expected 97)\n")
  lam <- n / (area.owin(W) / 1e6)
  cat(sprintf("Intensity (sites/km2): %.5f | CSR mean NN distance: %.2f km\n", lam, 0.5 / sqrt(lam)))
  
  # ---- 4. DISTANCE TO EDGE (cross-check) --------------------------------------
  safe("4 edge distance", {
    cat("\n==== 4. DISTANCE TO EDGE (spatstat bdist.points) ====\n")
    bd <- bdist.points(X) / 1000
    print(round(quantile(bd, c(0, .05, .1, .25, .5, .75, 1)), 2))
    if (!is.null(sch) && "dist_edge_km_signed" %in% names(sch)) {
      ref <- abs(as.numeric(sch$dist_edge_km_signed[match(S$site_id, sch$site_id)]))
      cat("Max abs diff vs script 02 edge distance (km):", round(max(abs(bd - ref), na.rm = TRUE), 3), "\n")
    } else cat("Script-02 edge-distance column not available: cross-check skipped.\n")
  })
  
  # ---- 5. ERODED AREA ----------------------------------------------------------
  fine_km <- NULL; fa <- NULL; r_half <- NA; r_quart <- NA
  safe("5 eroded area", {
    cat("\n==== 5. ERODED WINDOW AREA (sf buffer; spatstat cross-check) ====\n")
    rr_km <- c(1, 2, 3, 5, 7.5, 10, 15, 20, 25, 30, 40)
    er <- do.call(rbind, lapply(rr_km, function(r) {
      a <- buffer_area(win_sfc, r)
      data.frame(r_km = r, eroded_area_km2 = round(a, 1), fraction_of_window = round(a / A_sf, 3))
    }))
    er$spatstat_eroded_km2 <- tryCatch(round(eroded.areas(W, rr_km * 1000) / 1e6, 1), error = function(e) NA)
    print(er, row.names = FALSE)
    write.csv(er, file.path(out, "eroded_area_table.csv"), row.names = FALSE)
    fine_km <- seq(1, 60, by = 1)
    fa <- vapply(fine_km, function(r) buffer_area(win_sfc, r) / A_sf, 0)
    r_half  <<- fine_km[which(fa < 0.5)[1]]
    r_quart <<- fine_km[which(fa < 0.25)[1]]
    cat("\nr (km) at which eroded window first falls below 50% of area:", r_half,
        "| below 25%:", r_quart, "\n")
    cat("Descriptive geometry only. r_max is chosen at freeze (PROPOSED), not here.\n")
  })
  
  # ---- 6. PAIR COUNTS ----------------------------------------------------------
  safe("6 pair counts", {
    cat("\n==== 6. OBSERVED PAIRS WITHIN r (descriptive) ====\n")
    pd  <- pairdist(X) / 1000; pdC <- pairdist(XC) / 1000
    r_set <- c(0.5, 1, 1.5, 2, 3, 5, 7.5, 10, 15, 20, 30)
    cntp <- function(P, r) sum(P[upper.tri(P)] <= r)
    pc <- data.frame(
      r_km = r_set,
      pairs_A = sapply(r_set, function(r) cntp(pd, r)),
      pairs_C = sapply(r_set, function(r) cntp(pdC, r)),
      expected_neighbours_per_site_CSR_noedge = round(lam * pi * r_set^2, 2))
    pc$mean_observed_neighbours_per_site <- round(2 * pc$pairs_A / n, 2)
    print(pc, row.names = FALSE)
    write.csv(pc, file.path(out, "pair_counts.csv"), row.names = FALSE)
    cat("Observed vs CSR expectation here has no edge correction and no inhomogeneity.\n")
  })
  
  # ---- 7. DEM METADATA --------------------------------------------------------
  cat("\n==== 7. DEM GDAL METADATA (first 70 lines; full in dem_gdal_info.txt) ====\n")
  cat(head(dem_info, 70), sep = "\n")
  
  # ---- 8. SAVE ----------------------------------------------------------------
  safe("8 save", {
    saveRDS(list(W_m = W, X_m = X, X_km = X_km, XC_m = XC, marks = marks_tab,
                 crs_epsg = CRS_PROJ, window_route = route,
                 notes = "Window = dissolved boundary shapefile; n=197 NA; Dataset B retired; C = non-fluvial."),
            file.path(REPO_DIR, "data_derived", "ppp_NA_window.rds"))
    summ <- data.frame(item = c("window_area_km2", "window_parts", "window_holes", "window_vertices_sf",
                                "spatstat_route", "n_points_A", "n_points_C", "intensity_per_km2",
                                "r_eroded_below_50pct_km", "r_eroded_below_25pct_km"),
                       value = c(round(A_sf, 1), length(parts), sum(nrings - 1), nvert, route,
                                 npoints(X), npoints(XC), signif(lam, 4), r_half, r_quart))
    write.csv(summ, file.path(out, "window_summary.csv"), row.names = FALSE)
    cat("\n==== SAVED ====\ndata_derived/ppp_NA_window.rds and window_summary.csv\n")
  })
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_07.txt"))
cat("\nDONE. Outputs in:", out, "\n")