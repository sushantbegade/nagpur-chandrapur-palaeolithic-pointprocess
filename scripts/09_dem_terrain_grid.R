# =============================================================================
# 09_dem_terrain_grid.R   (v2: zero-fill outside window handled; OVERWRITES v1)
# Stage: 03_environment | Gate 4 (environment) support
# Purpose: diagnose zero-filled DEM cells outside the window, set them to NA,
#   compute terrain derivatives (slope, TRI) on the ORIGINAL EPSG:4326 grid,
#   reproject to EPSG:32644 at 50 m, build elevation + TPI(2 km) with a
#   valid-cell-normalised neighbourhood, aggregate to nested model grids
#   (250/500/1000 m), mask to window, summarise over window-CENTRE cells.
# Does NOT read site data. No site-by-covariate tables. No models.
# Outputs:
#   data_derived/env/<res>m/{elev,slope,tri,tpi2km,inwin,incenter}.tif  (git-ignored)
#   03_environment/outputs/: terrain_summary.csv, processing_log.csv,
#     grid_summary.txt, dem_zero_diagnostic.csv, sessionInfo_09.txt
#   docs/covariate_definitions.csv
# Not executed by assistant. Run top to bottom in RStudio. A few minutes.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
RES_SET     <- c(250, 500, 1000)       # model grids (m); 500 = primary
BASE_RES    <- 50                      # intermediate grid (m)
TPI_RADIUS  <- 2000                    # m
TPI_BASE    <- 250                     # m grid on which TPI is computed
MIN_VALID_W <- 0.05                    # TPI set NA if valid neighbourhood weight < this

while (sink.number() > 0) sink()
source(REPO_CONFIG)
suppressPackageStartupMessages({ library(terra); library(sf) })
sf_use_s2(FALSE)

out <- dir_of("03_environment"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
env_dir <- file.path(REPO_DIR, "data_derived", "env"); dir.create(env_dir, recursive = TRUE, showWarnings = FALSE)
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
getp <- function(l) file.path(RAW_DIR, man$path[man$layer == l])
stopifnot(all(RES_SET %% BASE_RES == 0), TPI_BASE %% BASE_RES == 0)

# ---- READ (before sink) ------------------------------------------------------
bnd    <- st_transform(st_make_valid(st_read(getp("Boundary"), quiet = TRUE)), CRS_PROJ)
win    <- st_union(st_geometry(bnd))
win_v  <- vect(st_sf(geometry = win))
win_ll <- vect(st_sf(geometry = st_transform(win, 4326)))
dem_raw <- rast(getp("DEM")); names(dem_raw) <- "elev"
bb <- st_bbox(win)
x0 <- floor(unname(bb["xmin"]) / 1000) * 1000; x1 <- ceiling(unname(bb["xmax"]) / 1000) * 1000
y0 <- floor(unname(bb["ymin"]) / 1000) * 1000; y1 <- ceiling(unname(bb["ymax"]) / 1000) * 1000
crs_txt <- paste0("EPSG:", CRS_PROJ)

logrows <- list()
lg <- function(step, detail) logrows[[length(logrows) + 1]] <<- data.frame(step = step, detail = detail, time = format(Sys.time()), stringsAsFactors = FALSE)
t0 <- Sys.time(); el <- function() round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
qs <- function(v) round(quantile(v, c(0, .01, .05, .25, .5, .75, .95, .99, 1), na.rm = TRUE), 3)

sink(file.path(out, "grid_summary.txt"), split = TRUE)
tryCatch({
  
  cat("==== 0. SETUP ====\nterra", as.character(packageVersion("terra")), "| CRS", crs_txt, "\n")
  cat("Master grid extent (m): x", x0, "-", x1, " y", y0, "-", y1, "\n")
  cat("DEM:", ncol(dem_raw), "x", nrow(dem_raw), "cells | res (deg):", signif(res(dem_raw)[1], 9), "\n")
  
  # ---- 0b. ZERO-FILL DIAGNOSTIC ------------------------------------------------
  cat("\n==== 0b. DEM ZERO DIAGNOSTIC (exact polygon on DEM grid) ====\n")
  inw_dem <- rasterize(win_ll, dem_raw, field = 1, touches = FALSE)
  ind  <- ifel(is.na(inw_dem), 0, 1)
  zero <- ifel(dem_raw == 0, 1, 0)
  n_all <- ncell(dem_raw)
  n_in  <- global(ind, "sum", na.rm = TRUE)[1, 1]
  z_all <- global(zero, "sum", na.rm = TRUE)[1, 1]
  z_in  <- global(zero * ind, "sum", na.rm = TRUE)[1, 1]
  min_in <- global(mask(dem_raw, inw_dem), "min", na.rm = TRUE)[1, 1]
  max_in <- global(mask(dem_raw, inw_dem), "max", na.rm = TRUE)[1, 1]
  n_out <- n_all - n_in; z_out <- z_all - z_in
  diag_tab <- data.frame(
    quantity = c("cells total", "cells inside window (centre)", "cells outside window",
                 "exact-zero cells total", "exact-zero cells inside window", "exact-zero cells outside window",
                 "nonzero cells outside window", "min elevation inside window (m)", "max elevation inside window (m)"),
    value = c(n_all, n_in, n_out, z_all, z_in, z_out, n_out - z_out, round(min_in, 2), round(max_in, 2)))
  print(diag_tab, row.names = FALSE)
  write.csv(diag_tab, file.path(out, "dem_zero_diagnostic.csv"), row.names = FALSE)
  cat(sprintf("Share of outside-window cells that are exactly zero: %.2f%%\n", 100 * z_out / n_out))
  if (z_in > 0) stop("Zero-valued DEM cells INSIDE the window (", z_in, "). Possible voids. Paste this output; nothing was modified.")
  dem <- ifel(is.na(inw_dem) & dem_raw == 0, NA, dem_raw); names(dem) <- "elev"
  lg("dem_clean", sprintf("%d exact-zero cells outside window set to NA; nonzero outside cells kept", z_out))
  cat("Zero cells outside the window set to NA:", z_out, "\n")
  
  # ---- 1. TERRAIN ON ORIGINAL GRID --------------------------------------------
  cat("\n==== 1. TERRAIN ON ORIGINAL (lon/lat) GRID ====\n")
  terr <- terrain(dem, v = c("slope", "TRI"), neighbors = 8, unit = "degrees")
  lg("terrain", "slope (Horn, 8 neighbours, degrees) and TRI (m) on cleaned original EPSG:4326 DEM")
  cat("Done. elapsed (min):", el(), "\n")
  
  # ---- 2. PROJECT TO 50 m UTM --------------------------------------------------
  cat("\n==== 2. PROJECT TO", BASE_RES, "m,", crs_txt, "====\n")
  tmpl <- rast(xmin = x0, xmax = x1, ymin = y0, ymax = y1, resolution = BASE_RES, crs = crs_txt)
  elev50 <- project(dem, tmpl, method = "bilinear")
  terr50 <- project(terr, tmpl, method = "bilinear")
  lg("project", sprintf("bilinear to %d m grid", BASE_RES))
  cat("Intermediate grid:", ncol(tmpl), "x", nrow(tmpl), "| elapsed (min):", el(), "\n")
  
  # ---- 3. TPI AND NESTED MODEL GRIDS ------------------------------------------
  cat("\n==== 3. NESTED MODEL GRIDS ====\n")
  e_tpi <- aggregate(elev50, fact = TPI_BASE / BASE_RES, fun = "mean", na.rm = TRUE)
  fm    <- focalMat(e_tpi, TPI_RADIUS, "circle")
  v_ok  <- ifel(is.na(e_tpi), 0, 1)
  e0    <- ifel(is.na(e_tpi), 0, e_tpi)
  num   <- focal(e0,   w = fm, fun = "sum", na.rm = TRUE)
  den   <- focal(v_ok, w = fm, fun = "sum", na.rm = TRUE)
  nb    <- ifel(den >= MIN_VALID_W, num / den, NA)
  tpi_b <- ifel(is.na(e_tpi), NA, e_tpi - nb); names(tpi_b) <- "tpi2km"
  cat(sprintf("TPI: %d m mean elevation minus VALID-cell mean in circle radius %d m (%dx%d weights); set NA where valid weight < %.2f\n",
              TPI_BASE, TPI_RADIUS, nrow(fm), ncol(fm), MIN_VALID_W))
  cat("Note: boundary cells use a truncated neighbourhood (inside window / valid cells only).\n")
  lg("tpi", sprintf("%d m mean elev minus valid-cell mean, circle %d m; min valid weight %.2f", TPI_BASE, TPI_RADIUS, MIN_VALID_W))
  
  summ <- list()
  for (res in RES_SET) {
    fct <- res / BASE_RES
    d <- file.path(env_dir, paste0(res, "m")); dir.create(d, showWarnings = FALSE)
    L <- list(
      elev   = aggregate(elev50, fact = fct, fun = "mean", na.rm = TRUE),
      slope  = aggregate(terr50[["slope"]], fact = fct, fun = "mean", na.rm = TRUE),
      tri    = aggregate(terr50[["TRI"]],   fact = fct, fun = "mean", na.rm = TRUE),
      tpi2km = if (res == TPI_BASE) tpi_b else aggregate(tpi_b, fact = res / TPI_BASE, fun = "mean", na.rm = TRUE))
    grid <- L$elev
    inwin    <- rasterize(win_v, grid, field = 1, touches = TRUE);  names(inwin) <- "inwin"
    incenter <- rasterize(win_v, grid, field = 1, touches = FALSE); names(incenter) <- "incenter"
    ic <- !is.na(values(incenter, mat = FALSE)); it <- !is.na(values(inwin, mat = FALSE))
    cat(sprintf("\n-- %d m grid: %d x %d | window-centre cells: %d (%.0f km2) | touching cells: %d (%.0f km2)\n",
                res, ncol(grid), nrow(grid), sum(ic), sum(ic) * res^2 / 1e6, sum(it), sum(it) * res^2 / 1e6))
    for (nm in names(L)) {
      r <- mask(L[[nm]], win_v, touches = TRUE); names(r) <- nm
      writeRaster(r, file.path(d, paste0(nm, ".tif")), overwrite = TRUE, datatype = "FLT4S")
      vv <- values(r, mat = FALSE)
      v <- vv[ic]; q <- qs(v)
      na_c <- round(100 * mean(is.na(v)), 3); na_t <- round(100 * mean(is.na(vv[it])), 3)
      cat(sprintf("%-7s NA centre %5.3f%% touch %5.3f%% | min %9.3f | p5 %9.3f | median %9.3f | p95 %9.3f | max %9.3f\n",
                  nm, na_c, na_t, q["0%"], q["5%"], q["50%"], q["95%"], q["100%"]))
      summ[[length(summ) + 1]] <- data.frame(res_m = res, layer = nm, n_centre_cells = sum(ic),
                                             pct_NA_centre = na_c, pct_NA_touch = na_t, min = q["0%"], p1 = q["1%"], p5 = q["5%"], p25 = q["25%"],
                                             median = q["50%"], p75 = q["75%"], p95 = q["95%"], p99 = q["99%"], max = q["100%"], row.names = NULL)
    }
    writeRaster(inwin,    file.path(d, "inwin.tif"),    overwrite = TRUE, datatype = "INT1U")
    writeRaster(incenter, file.path(d, "incenter.tif"), overwrite = TRUE, datatype = "INT1U")
    lg("aggregate", sprintf("%d m grid: mean aggregation by factor %d from %d m", res, fct, BASE_RES))
  }
  ts <- do.call(rbind, summ)
  write.csv(ts, file.path(out, "terrain_summary.csv"), row.names = FALSE)
  cat("\nElapsed (min):", el(), "\n")
  
  # ---- 4. COVARIATE DEFINITIONS (not frozen) ----------------------------------
  cat("\n==== 4. COVARIATE DEFINITIONS ====\n")
  defs <- data.frame(
    covariate = c("elev", "slope", "tpi2km", "tri", "ldist_riv3", "ldist_riv5", "lith5", "geomorph3"),
    domain = c("terrain (regional)", "terrain (local)", "terrain (local position)", "terrain (local; swap-in only)",
               "hydrology (local)", "hydrology (regional)", "geological / raw-material proxy", "geomorphology"),
    derivation = c("DEM (zero fill outside window removed), mean to grid", "DEM 30 m Horn slope on lon/lat grid, mean to grid",
                   sprintf("%d m mean elevation minus valid-cell mean in %d m circle", TPI_BASE, TPI_RADIUS),
                   "DEM TRI (mean abs difference to 8 neighbours), mean to grid",
                   "NGDR Rivers, Strahler order >= 3", "NGDR Rivers, Strahler order >= 5",
                   "NGDR Lithology, 5-class collapse", "NGDR Geomorphology, 3 classes; non-landform polygons recoded to nearest landform"),
    transform = c("z-score; linear + quadratic", "log(1+slope_deg), z-score [PROPOSED, pending clean-data check]", "z-score", "z-score",
                  "log(1+d_km), z-score", "log(1+d_km), z-score", "factor", "factor"),
    df = c(2, 1, 1, 0, 1, 1, 4, 2),
    E1 = c(1, 1, 1, 0, 1, 1, 0, 0), E2 = c(1, 1, 1, 0, 1, 1, 1, 0),
    E3 = c(1, 1, 1, 0, 1, 1, 1, 1), E4 = c(1, 1, 1, 0, 0, 0, 1, 0),
    caveat = c("modern DEM; Pleistocene surface not reconstructed; DEM product unknown", "30 m source; smoothed by aggregation",
               "scale chosen a priori; truncated neighbourhood at window edge", "enters only as slope replacement in a sensitivity variant",
               "modern mapped network; also encodes exposure/preservation", "modern mapped network; also encodes exposure/preservation",
               "lithology is not demonstrated raw-material availability", "alluvial class carries preservation role"),
    status = "APPROVED (author-delegated 2026-11-06); not FROZEN", stringsAsFactors = FALSE)
  write.csv(defs, file.path(REPO_DIR, "docs", "covariate_definitions.csv"), row.names = FALSE)
  cat("df totals: E1", sum(defs$df * defs$E1), "| E2", sum(defs$df * defs$E2), "| E3", sum(defs$df * defs$E3),
      "| E4", sum(defs$df * defs$E4), "(cap ~12)\n")
  write.csv(do.call(rbind, logrows), file.path(out, "processing_log.csv"), row.names = FALSE)
  cat("\nSaved rasters to data_derived/env/<res>m/ ; summary in terrain_summary.csv\n")
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_09.txt"))
cat("\nDONE. Outputs in:", out, "\n")