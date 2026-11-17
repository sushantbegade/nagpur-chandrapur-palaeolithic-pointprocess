# =============================================================================
# 01_inspect_inputs.R
# Stage: 01_data_audit | Purpose: metadata inventory of the LISTED inputs only
# No recursive scanning. Only files in MANIFEST are opened.
# Outputs (01_data_audit/outputs/): input_manifest_check.csv, input_inventory.csv,
#   vector_attributes.txt, boundary_check.txt, site_table_profile.txt,
#   sessionInfo_01.txt ; also docs/input_manifest.csv (reused by later scripts)
# Not executed by assistant. Run top to bottom in RStudio.
# =============================================================================

# ---- EDIT ONLY THIS LINE if repo path differs --------------------------------
REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
# ------------------------------------------------------------------------------

source(REPO_CONFIG)   # defines RAW_DIR (= Nagpur-Chandrapur Enhanced Dataset), CRS_PROJ, dir_of()
suppressPackageStartupMessages({ library(terra); library(sf) })

out <- dir_of("01_data_audit")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- MANIFEST: paths relative to RAW_DIR -------------------------------------
sent <- c("B01","B02","B03","B04","B05","B06","B07","B08","B09","B11","B12","B8A")
manifest <- rbind(
  data.frame(type = "raster", layer = "DEM", domain = "terrain",
             path = "DEM/DEM.tif"),
  data.frame(type = "raster", layer = paste0("S2_", sent), domain = "sentinel2_band",
             path = sprintf("Sentinel 2A Multispectral Data/%s_Nagpur_Chandrapur.tif", sent)),
  data.frame(type = "raster",
             layer = paste0("S2_", c("BSI","MNDWI","MSAVI","NDBI","NDVI","NDWI","SAVI")),
             domain = "sentinel2_index",
             path = sprintf("Sentinel 2A Multispectral Data/%s.tif",
                            c("BSI","MNDWI","MSAVI","NDBI","NDVI","NDWI","SAVI"))),
  data.frame(type = "raster",
             layer = paste0("Soil_", c("Depth","Erosion","Productivity","Slope","Texture")),
             domain = "soil",
             path = sprintf("Soil Data/Soil %s.tif",
                            c("Depth","Erosion","Productivity","Slope","Texture"))),
  data.frame(type = "vector", layer = "Boundary", domain = "window",
             path = "District Boundary/Nagpur_and_Chandrapur_Boundry.shp"),
  data.frame(type = "vector", layer = "Geology", domain = "geology",
             path = "Geology/Geology.shp"),
  data.frame(type = "vector", layer = "Geomorphology", domain = "geomorphology",
             path = "Geomorphology/Geomorphology.shp"),
  data.frame(type = "vector", layer = "Lithology", domain = "lithology",
             path = "Lithology/Lithology.shp"),
  data.frame(type = "vector", layer = "MineralDeposit", domain = "raw_material_candidate",
             path = "Mineral Deposit/Mineral Deposit.shp"),
  data.frame(type = "vector", layer = "Rivers", domain = "hydrology",
             path = "Rivers/Rivers.shp"),
  data.frame(type = "vector", layer = "Waterbody", domain = "hydrology",
             path = "Waterbody/Waterbody.shp"),
  data.frame(type = "table", layer = "SiteTable_csv", domain = "sites",
             path = "Nagpur Chandrapur Palaeolithic Site Raw Data.csv")
)
manifest$full_path <- file.path(RAW_DIR, manifest$path)
manifest$exists    <- file.exists(manifest$full_path)

dir.create(file.path(REPO_DIR, "docs"), showWarnings = FALSE)
write.csv(manifest[, c("type","layer","domain","path")],
          file.path(REPO_DIR, "docs", "input_manifest.csv"), row.names = FALSE)
write.csv(manifest, file.path(out, "input_manifest_check.csv"), row.names = FALSE)

cat("==== MANIFEST CHECK:", sum(manifest$exists), "of", nrow(manifest), "found ====\n")
if (any(!manifest$exists)) { cat("MISSING:\n"); print(manifest$path[!manifest$exists]) }

# ---- 1. RASTERS (manifest only) ----------------------------------------------
na_row <- function(m, status) data.frame(
  type = m$type, layer = m$layer, file = m$path, status = status, rotated = NA,
  crs_epsg = NA, crs_name = NA, res_x = NA, res_y = NA, nlyr = NA, ncol = NA,
  nrow = NA, xmin = NA, xmax = NA, ymin = NA, ymax = NA, pct_NA_sample = NA,
  min = NA, max = NA, geom_type = NA, n_invalid_geom = NA, fields = NA,
  n_features = NA, stringsAsFactors = FALSE)

inspect_raster <- function(m) {
  if (!m$exists) return(na_row(m, "MISSING"))
  cat("Reading raster:", m$path, "\n")
  rotated <- FALSE
  r <- withCallingHandlers(
    tryCatch(rast(m$full_path), error = function(e) NULL),
    warning = function(w) {
      if (grepl("rotated", conditionMessage(w), ignore.case = TRUE)) rotated <<- TRUE
      invokeRestart("muffleWarning")
    })
  if (is.null(r)) return(na_row(m, "READ_ERROR"))
  
  st <- tryCatch({
    s <- spatSample(r[[1]], size = 20000, method = "regular", warn = FALSE)[, 1]
    list(na = round(100 * mean(is.na(s)), 2),
         mn = suppressWarnings(min(s, na.rm = TRUE)),
         mx = suppressWarnings(max(s, na.rm = TRUE)))
  }, error = function(e) list(na = NA, mn = NA, mx = NA))
  cd <- tryCatch(crs(r, describe = TRUE),
                 error = function(e) data.frame(code = NA, name = NA))
  
  data.frame(
    type = m$type, layer = m$layer, file = m$path,
    status = if (rotated) "ROTATED" else "ok", rotated = rotated,
    crs_epsg = cd$code, crs_name = cd$name,
    res_x = res(r)[1], res_y = res(r)[2], nlyr = nlyr(r),
    ncol = ncol(r), nrow = nrow(r),
    xmin = xmin(r), xmax = xmax(r), ymin = ymin(r), ymax = ymax(r),
    pct_NA_sample = st$na, min = st$mn, max = st$mx,
    geom_type = NA, n_invalid_geom = NA, fields = NA, n_features = NA,
    stringsAsFactors = FALSE)
}
mr <- manifest[manifest$type == "raster", ]
rast_inv <- do.call(rbind, lapply(seq_len(nrow(mr)), function(i) inspect_raster(mr[i, ])))

# ---- 2. VECTORS (manifest only) ----------------------------------------------
inspect_vector <- function(m) {
  if (!m$exists) return(na_row(m, "MISSING"))
  cat("Reading vector:", m$path, "\n")
  v <- tryCatch(st_read(m$full_path, quiet = TRUE), error = function(e) NULL)
  if (is.null(v)) return(na_row(m, "READ_ERROR"))
  bb <- as.numeric(st_bbox(v))
  inval <- tryCatch(sum(!st_is_valid(v), na.rm = TRUE), error = function(e) NA)
  data.frame(
    type = m$type, layer = m$layer, file = m$path, status = "ok", rotated = NA,
    crs_epsg = st_crs(v)$epsg, crs_name = st_crs(v)$Name,
    res_x = NA, res_y = NA, nlyr = NA, ncol = NA, nrow = NA,
    xmin = bb[1], ymin = bb[2], xmax = bb[3], ymax = bb[4],
    pct_NA_sample = NA, min = NA, max = NA,
    geom_type = paste(unique(as.character(st_geometry_type(v))), collapse = "/"),
    n_invalid_geom = inval,
    fields = paste(setdiff(names(v), attr(v, "sf_column")), collapse = "; "),
    n_features = nrow(v), stringsAsFactors = FALSE)
}
mv <- manifest[manifest$type == "vector", ]
vec_inv <- do.call(rbind, lapply(seq_len(nrow(mv)), function(i) inspect_vector(mv[i, ])))

# ---- 3. MERGE + REPORT -------------------------------------------------------
inv <- rbind(rast_inv, vec_inv)
write.csv(inv, file.path(out, "input_inventory.csv"), row.names = FALSE)

cat("\n==== INPUT INVENTORY ====\n")
print(inv[, c("type","layer","status","crs_epsg","res_x","res_y",
              "n_features","pct_NA_sample")], right = FALSE)
cat("\n==== ROTATED RASTERS ====\n");   print(inv$layer[which(inv$status == "ROTATED")])
cat("\n==== READ ERRORS ====\n");       print(inv$layer[which(inv$status == "READ_ERROR")])
cat("\n==== MISSING ====\n");           print(inv$layer[which(inv$status == "MISSING")])

# ---- 4. VECTOR ATTRIBUTES + CLASS COUNTS -------------------------------------
sink(file.path(out, "vector_attributes.txt"))
tryCatch({
  for (i in seq_len(nrow(mv))) {
    m <- mv[i, ]; if (!m$exists) next
    v <- tryCatch(st_read(m$full_path, quiet = TRUE), error = function(e) NULL)
    if (is.null(v)) next
    d <- st_drop_geometry(v)
    cat("\n====", m$layer, "(", m$path, ") ====\n")
    cat("Features:", nrow(d), " Fields:", ncol(d), "\n")
    if (ncol(d) == 0) { cat("(no attribute fields)\n"); next }
    str(d, vec.len = 3, give.attr = FALSE)
    for (cn in names(d)[sapply(d, function(x) is.character(x) || is.factor(x))]) {
      tb <- sort(table(d[[cn]], useNA = "ifany"), decreasing = TRUE)
      cat("\nField '", cn, "': ", length(tb), " unique values; top 15:\n", sep = "")
      print(head(tb, 15))
    }
  }
}, finally = sink())
cat("\n==== VECTOR ATTRIBUTES ====\n")
cat(readLines(file.path(out, "vector_attributes.txt")), sep = "\n")

# ---- 5. BOUNDARY CHECK (area vs ~21,335 km2 in design) -----------------------
bm <- manifest[manifest$layer == "Boundary", ]
sink(file.path(out, "boundary_check.txt"))
tryCatch({
  b <- st_read(bm$full_path, quiet = TRUE)
  cat("Boundary file:", bm$path, "\nFeatures:", nrow(b), "\n")
  cat("Source CRS:", st_crs(b)$Name, " EPSG:", st_crs(b)$epsg, "\n")
  cat("Fields:", paste(setdiff(names(b), attr(b, "sf_column")), collapse = "; "), "\n")
  bp <- st_transform(st_make_valid(b), CRS_PROJ)
  a  <- as.numeric(st_area(bp)) / 1e6
  if (nrow(b) > 1) { cat("\nArea per feature (km2):\n"); print(round(a, 1)) }
  cat("\nTotal area (km2), EPSG", CRS_PROJ, "=", round(sum(a), 1), "\n")
  cat("Design-document value: ~21,335 km2\nDifference (km2):", round(sum(a) - 21335, 1), "\n")
}, error = function(e) cat("BOUNDARY CHECK FAILED:", conditionMessage(e), "\n"),
finally = sink())
cat("\n==== BOUNDARY CHECK ====\n")
cat(readLines(file.path(out, "boundary_check.txt")), sep = "\n")

# ---- 6. SITE TABLE PROFILE (CSV only) ----------------------------------------
sm <- manifest[manifest$layer == "SiteTable_csv", ]
sites <- tryCatch(
  read.csv(sm$full_path, stringsAsFactors = FALSE, check.names = FALSE,
           fileEncoding = "UTF-8-sig"),
  error = function(e) read.csv(sm$full_path, stringsAsFactors = FALSE, check.names = FALSE))

sink(file.path(out, "site_table_profile.txt"))
tryCatch({
  cat("Rows:", nrow(sites), " Cols:", ncol(sites), "\n\nColumn names:\n")
  print(names(sites))
  cat("\nStructure:\n"); str(sites, give.attr = FALSE)
  cat("\nNA / empty-string count per column:\n")
  print(sapply(sites, function(x) sum(is.na(x) | trimws(as.character(x)) == "")))
  cat("\nFirst 10 rows:\n"); print(head(sites, 10))
  cat("\nExact duplicate rows:", sum(duplicated(sites)), "\n")
  for (cn in names(sites)[sapply(sites, is.character)]) {
    if (length(unique(sites[[cn]])) <= 25) {
      cat("\nField '", cn, "' value counts:\n", sep = ""); print(table(sites[[cn]], useNA = "ifany"))
    }
  }
}, finally = sink())
cat("\n==== SITE TABLE PROFILE ====\n")
cat(readLines(file.path(out, "site_table_profile.txt")), sep = "\n")

# ---- 7. SESSION INFO ---------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_01.txt"))
cat("\nDONE. Outputs in:", out, "\n")