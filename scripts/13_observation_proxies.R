# =============================================================================
# 13_observation_proxies.R
# Stage: 04_observation | Gate 5 (observation) support
# Purpose: observation proxy layers from EXTERNAL layers only (OSM roads, ESA
#   WorldCover, window edge), O1-O3 log-offsets, composites for the free-coefficient
#   variant, overlap diagnostics with environmental covariates (window cells only).
# Does NOT read site data. No site-by-covariate tables. No models.
# Requires outputs of scripts 09-12 (env grids, z_* layers, WorldCover tiles, roads).
# Outputs:
#   data_derived/env/<res>m/{f_built,f_water,f_tree,f_crop,dist_road_km,dist_edge_km,
#     S_access,S_open,S_core,S2,S3,logE_O1,logE_O2,logE_O3,z_S2,z_S3}.tif
#   04_observation/outputs/: observation_report.txt, proxy_summary.csv, landcover_window.csv,
#     scenario_range.csv, proxy_env_correlations.csv, proxy_vif.csv,
#     processing_log_13.csv, sessionInfo_13.txt
# Not executed by assistant. Run top to bottom in RStudio. Long run (heavy rasters).
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
RES_SET <- c(250, 500, 1000); PRIMARY <- 500; BASE_RES <- 50; DIST_RES <- 250
ACCESS_DECAY_KM <- 2; EDGE_SAT_KM <- 20; TREE_WEIGHT <- 0.5; R_O2 <- 2; R_O3 <- 4
ROAD_MARGIN_M <- 20000; LC_AGG <- 6; LC_SNAP <- 0.0005; LC_MARGIN_DEG <- 0.02
VIF_MAX <- 5

while (sink.number() > 0) sink()
source(REPO_CONFIG)
suppressPackageStartupMessages({ library(terra); library(sf) })
sf_use_s2(FALSE)

out <- dir_of("04_observation"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
env_dir <- file.path(REPO_DIR, "data_derived", "env")
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
getp <- function(l) file.path(RAW_DIR, man$path[man$layer == l])
crs_txt <- paste0("EPSG:", CRS_PROJ)

# ---- READ (before sink) ------------------------------------------------------
bnd   <- st_transform(st_make_valid(st_read(getp("Boundary"), quiet = TRUE)), CRS_PROJ)
win   <- st_union(st_geometry(bnd)); win_v <- vect(st_sf(geometry = win))
bb    <- st_bbox(st_transform(win, 4326))
roads <- st_transform(st_read(getp("Roads_OSM"), quiet = TRUE), CRS_PROJ)
wc_files <- file.path(RAW_DIR, man$path[grepl("^WorldCover_", man$layer)])
stopifnot(length(wc_files) >= 1, all(file.exists(wc_files)))
log_f <- file.path(REPO_DIR, "docs", "external_download_log.csv")
dlog  <- if (file.exists(log_f)) read.csv(log_f, stringsAsFactors = FALSE) else NULL
for (r in RES_SET) {
  d <- file.path(env_dir, paste0(r, "m"))
  stopifnot(all(file.exists(file.path(d, paste0(c("inwin", "incenter", "z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km",
                                                  "z_ldist_riv3", "z_ldist_riv5", "lith5", "geomorph3"), ".tif")))))
}
logrows <- list()
lg <- function(step, detail) logrows[[length(logrows) + 1]] <<- data.frame(step = step, detail = detail, time = format(Sys.time()), stringsAsFactors = FALSE)
t0 <- Sys.time(); el <- function() round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
qs <- function(v) round(quantile(v, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE), 3)
vif_of <- function(M) { R <- tryCatch(cor(M), error = function(e) NULL)
if (is.null(R)) return(setNames(rep(NA_real_, ncol(M)), colnames(M)))
setNames(tryCatch(diag(solve(R)), error = function(e) rep(Inf, ncol(M))), colnames(M)) }
eta2 <- function(x, g) { m <- mean(x); st <- sum((x - m)^2); mg <- tapply(x, g, mean); ng <- tapply(x, g, length); sum(ng * (mg - m)^2) / st }

sink(file.path(out, "observation_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== 0. PARAMETERS AND PROVENANCE ====\n")
  cat(sprintf("access decay %g km | edge saturation %g km | tree weight %g | O2 factor %g | O3 factor %g\n",
              ACCESS_DECAY_KM, EDGE_SAT_KM, TREE_WEIGHT, R_O2, R_O3))
  if (!is.null(dlog)) { cat("External download log (checksums, notes):\n")
    print(data.frame(layer = dlog$layer, MB = round(dlog$bytes / 1e6, 1), md5 = dlog$md5, licence = dlog$licence,
                     note = substr(dlog$note, 1, 140)), row.names = FALSE) } else cat("External download log not found.\n")
  cat("Roads loaded:", nrow(roads), "ways | WorldCover files:", length(wc_files), "\n")
  
  # ---- 1. LAND-COVER FRACTIONS (lon/lat ~50 m) -----------------------------------
  cat("\n==== 1. WORLDCOVER FRACTIONS ====\n")
  mos <- tryCatch(vrt(wc_files), error = function(e) { cat("vrt() failed:", conditionMessage(e), "-> merge fallback\n"); NULL })
  if (is.null(mos)) mos <- do.call(merge, lapply(wc_files, rast))
  x0 <- floor((unname(bb["xmin"]) - LC_MARGIN_DEG) / LC_SNAP) * LC_SNAP; x1 <- ceiling((unname(bb["xmax"]) + LC_MARGIN_DEG) / LC_SNAP) * LC_SNAP
  y0 <- floor((unname(bb["ymin"]) - LC_MARGIN_DEG) / LC_SNAP) * LC_SNAP; y1 <- ceiling((unname(bb["ymax"]) + LC_MARGIN_DEG) / LC_SNAP) * LC_SNAP
  rc <- crop(mos, ext(x0, x1, y0, y1))
  cat("Cropped mosaic:", ncol(rc), "x", nrow(rc), "cells; extent", round(as.vector(ext(rc)), 4), "\n")
  CLS <- c(f_built = 50, f_water = 80, f_tree = 10, f_crop = 40)
  fr_ll <- rast(lapply(names(CLS), function(nm) { cat("  aggregating", nm, "... elapsed (min):", el(), "\n"); flush.console()
    a <- aggregate(rc == CLS[[nm]], fact = LC_AGG, fun = "mean", na.rm = TRUE); names(a) <- nm; a }))
  cat("Fraction raster (lon/lat):", ncol(fr_ll), "x", nrow(fr_ll), "| res (deg):", signif(res(fr_ll)[1], 5), "\n")
  lg("landcover", "WorldCover v200 2021 class fractions (built-up 50, water 80, tree 10, crop 40), 10 m -> factor 6 mean")
  
  r250 <- rast(file.path(env_dir, "250m", "inwin.tif")); e <- ext(r250)
  tmpl50 <- rast(xmin = xmin(e), xmax = xmax(e), ymin = ymin(e), ymax = ymax(e), resolution = BASE_RES, crs = crs_txt)
  fr50 <- project(fr_ll, tmpl50, method = "bilinear")
  cat("Projected to", BASE_RES, "m UTM. elapsed (min):", el(), "\n")
  
  # ---- 2. ROAD AND EDGE DISTANCE (250 m) --------------------------------------
  cat("\n==== 2. ROAD AND EDGE DISTANCE ====\n")
  tmpl_big <- rast(xmin = xmin(e) - ROAD_MARGIN_M, xmax = xmax(e) + ROAD_MARGIN_M,
                   ymin = ymin(e) - ROAD_MARGIN_M, ymax = ymax(e) + ROAD_MARGIN_M, resolution = DIST_RES, crs = crs_txt)
  rr <- rasterize(vect(st_sf(geometry = st_geometry(roads))), tmpl_big, field = 1, touches = TRUE)
  droad250 <- crop(distance(rr) / 1000, e)
  in250 <- rasterize(win_v, r250, field = 1, touches = FALSE)
  dedge250 <- distance(ifel(is.na(in250), 1, NA)) / 1000
  if (!compareGeom(droad250, r250, stopOnError = FALSE) || !compareGeom(dedge250, r250, stopOnError = FALSE))
    stop("Grid mismatch for 250 m distance layers. Paste this output.")
  cat("Road distance (250 m grid) and edge distance built. elapsed (min):", el(), "\n")
  lg("distance", "road distance with 20 km margin; edge distance from outside-cell mask")
  
  # ---- 3. PROXIES, OFFSETS, DIAGNOSTICS PER RESOLUTION ------------------------
  sum_rows <- list(); lc_rows <- list(); rng_rows <- list(); cor_rows <- list(); vif_rows <- list()
  for (res in RES_SET) {
    cat(sprintf("\n=================== %d m GRID ===================\n", res))
    d <- file.path(env_dir, paste0(res, "m")); rdl <- function(n) rast(file.path(d, paste0(n, ".tif")))
    inw <- rdl("inwin"); ic <- !is.na(values(rdl("incenter"), mat = FALSE))
    agg <- function(x, base) if (res > base) aggregate(x, fact = res / base, fun = "mean", na.rm = TRUE) else x
    fr <- agg(fr50, BASE_RES); droad <- agg(droad250, DIST_RES); dedge <- agg(dedge250, DIST_RES)
    for (x in list(fr, droad, dedge)) if (!compareGeom(x, inw, stopOnError = FALSE)) stop("Grid mismatch at ", res, " m. Paste this output.")
    names(droad) <- "dist_road_km"; names(dedge) <- "dist_edge_km"
    S_access <- exp(-droad / ACCESS_DECAY_KM); names(S_access) <- "S_access"
    S_open <- clamp(1 - (fr[["f_built"]] + fr[["f_water"]] + TREE_WEIGHT * fr[["f_tree"]]), 0, 1); names(S_open) <- "S_open"
    S_core <- clamp(dedge / EDGE_SAT_KM, 0, 1); names(S_core) <- "S_core"
    S2 <- (S_access + S_open) / 2; names(S2) <- "S2"
    S3 <- (S_access + S_open + S_core) / 3; names(S3) <- "S3"
    logE_O1 <- S2 * 0; names(logE_O1) <- "logE_O1"
    logE_O2 <- log(R_O2) * S2; names(logE_O2) <- "logE_O2"
    logE_O3 <- log(R_O3) * S3; names(logE_O3) <- "logE_O3"
    L <- c(as.list(fr), list(droad, dedge, S_access, S_open, S_core, S2, S3, logE_O1, logE_O2, logE_O3))
    names(L) <- c(names(fr), "dist_road_km", "dist_edge_km", "S_access", "S_open", "S_core", "S2", "S3", "logE_O1", "logE_O2", "logE_O3")
    vals <- list()
    for (nm in names(L)) {
      r <- mask(L[[nm]], inw); names(r) <- nm
      writeRaster(r, file.path(d, paste0(nm, ".tif")), overwrite = TRUE, datatype = "FLT4S")
      vals[[nm]] <- values(r, mat = FALSE)[ic]
      q <- qs(vals[[nm]])
      sum_rows[[length(sum_rows) + 1]] <- data.frame(res_m = res, layer = nm, NA_centre = sum(is.na(vals[[nm]])),
                                                     min = q["0%"], p5 = q["5%"], p25 = q["25%"], median = q["50%"], p75 = q["75%"], p95 = q["95%"], max = q["100%"], row.names = NULL)
    }
    for (nm in c("S2", "S3")) {
      v <- vals[[nm]]; z <- (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
      zr <- mask((L[[nm]] - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE), inw); names(zr) <- paste0("z_", nm)
      writeRaster(zr, file.path(d, paste0("z_", nm, ".tif")), overwrite = TRUE, datatype = "FLT4S")
      vals[[paste0("z_", nm)]] <- z
    }
    cat("NA in centre cells (any proxy layer):", sum(sapply(vals, function(x) sum(is.na(x)))), "\n")
    tb <- do.call(rbind, lapply(sum_rows[sapply(sum_rows, function(x) x$res_m[1] == res)], identity))
    print(tb[tb$layer %in% c("f_built", "f_water", "f_tree", "f_crop", "dist_road_km", "dist_edge_km", "S_access", "S_open", "S_core", "S2", "S3"),
             c("layer", "min", "p5", "median", "p95", "max")], row.names = FALSE)
    
    other <- 1 - (vals$f_built + vals$f_water + vals$f_tree + vals$f_crop)
    lc <- data.frame(res_m = res, class = c("tree cover", "cropland", "built-up", "permanent water", "other (grass/shrub/bare/wetland)"),
                     mean_fraction = round(c(mean(vals$f_tree), mean(vals$f_crop), mean(vals$f_built), mean(vals$f_water), mean(other)), 4))
    cat("\nWindow land-cover composition (mean fraction over centre cells):\n"); print(lc, row.names = FALSE)
    lc_rows[[length(lc_rows) + 1]] <- lc
    
    Eq <- function(lg) { e <- exp(lg); c(p5 = unname(quantile(e, .05, na.rm = TRUE)), p95 = unname(quantile(e, .95, na.rm = TRUE))) }
    rg <- do.call(rbind, lapply(c("logE_O1", "logE_O2", "logE_O3"), function(nm) { q <- Eq(vals[[nm]])
    data.frame(res_m = res, scenario = sub("logE_", "", nm), E_p5 = round(q["p5"], 3), E_p95 = round(q["p95"], 3),
               ratio_p95_p5 = round(q["p95"] / q["p5"], 3), E_max_over_min = round(exp(max(vals[[nm]], na.rm = TRUE) - min(vals[[nm]], na.rm = TRUE)), 3), row.names = NULL) }))
    cat("\nScenario offset range over centre cells (E = exp(logE)):\n"); print(rg, row.names = FALSE)
    rng_rows[[length(rng_rows) + 1]] <- rg
    
    # ---- overlap with environment ----
    zn <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5")
    D <- as.data.frame(lapply(zn, function(n) values(rdl(n), mat = FALSE)[ic])); names(D) <- zn
    lv <- values(rdl("lith5"), mat = FALSE)[ic]; gv <- values(rdl("geomorph3"), mat = FALSE)[ic]
    D$lith_basalt <- as.numeric(lv == 1); D$lith_sandstone <- as.numeric(lv == 3); D$lith_carbonate <- as.numeric(lv == 4)
    D$lith_other <- as.numeric(lv == 5); D$geom_dissected <- as.numeric(gv == 2); D$geom_alluvial <- as.numeric(gv == 3)
    E3 <- c(zn, "lith_basalt", "lith_sandstone", "lith_carbonate", "lith_other", "geom_dissected", "geom_alluvial")
    px <- c("S_access", "S_open", "S_core", "S2", "S3")
    Cm <- sapply(px, function(p) sapply(zn, function(z) cor(D[[z]], vals[[p]])))
    cat("\nCorrelation of proxies with environmental layers (rows) :\n"); print(round(Cm, 2))
    cor_rows[[length(cor_rows) + 1]] <- data.frame(res_m = res, layer = rownames(Cm), round(Cm, 3), row.names = NULL)
    et <- data.frame(proxy = px, eta2_lith5 = sapply(px, function(p) round(eta2(vals[[p]], lv), 3)),
                     eta2_geomorph3 = sapply(px, function(p) round(eta2(vals[[p]], gv), 3)), row.names = NULL)
    cat("\nEta-squared of proxies by categorical layers:\n"); print(et, row.names = FALSE)
    for (zs in c("z_S2", "z_S3")) {
      v <- vif_of(as.matrix(cbind(D[, E3], setNames(data.frame(vals[[zs]]), zs))))
      cat(sprintf("E3 + %s: max VIF %.2f (%s); VIF of %s = %.2f\n", zs, max(v), names(which.max(v)), zs, v[zs]))
      vif_rows[[length(vif_rows) + 1]] <- data.frame(res_m = res, model = paste0("E3+", zs), variable = names(v), vif = round(as.numeric(v), 3), row.names = NULL)
    }
    lg("proxies", sprintf("%d m: proxies, offsets, z composites saved; diagnostics computed", res))
  }
  write.csv(do.call(rbind, sum_rows), file.path(out, "proxy_summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, lc_rows), file.path(out, "landcover_window.csv"), row.names = FALSE)
  write.csv(do.call(rbind, rng_rows), file.path(out, "scenario_range.csv"), row.names = FALSE)
  write.csv(do.call(rbind, cor_rows), file.path(out, "proxy_env_correlations.csv"), row.names = FALSE)
  write.csv(do.call(rbind, vif_rows), file.path(out, "proxy_vif.csv"), row.names = FALSE)
  write.csv(do.call(rbind, logrows), file.path(out, "processing_log_13.csv"), row.names = FALSE)
  cat("\nElapsed (min):", el(), "\nSaved layers to data_derived/env/<res>m/ ; tables in 04_observation/outputs\n")
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_13.txt"))
cat("\nDONE. Outputs in:", out, "\n")