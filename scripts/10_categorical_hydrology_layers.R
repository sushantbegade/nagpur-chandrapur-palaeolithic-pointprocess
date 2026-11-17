# =============================================================================
# 10_categorical_hydrology_layers.R
# Stage: 03_environment | Gate 4 (environment) support
# Purpose: lithology 5-class layer, geomorphology 3-class layer (non-landform
#   polygons recoded by a pre-specified rule), distance to Strahler order >=3 and
#   >=5 streams, extreme-slope check against external layers.
# Does NOT read site data. No site-by-covariate tables. No models.
# Requires: 09_dem_terrain_grid.R outputs (data_derived/env/<res>m/*.tif) and
#           01_data_audit/outputs/polygon_class_areas.csv
# Outputs:
#   data_derived/env/<res>m/{lith5,geomorph3,dist_riv3_km,dist_riv5_km,ldist_riv3,ldist_riv5}.tif
#   03_environment/outputs/: categorical_hydrology_report.txt, class_summary.csv,
#     geomorph_recode_table.csv, distance_summary.csv, slope_extreme_check.csv,
#     processing_log_10.csv, sessionInfo_10.txt
#   docs/covariate_definitions.csv  (slope transform row updated)
# Not executed by assistant. Run top to bottom in RStudio. Several minutes.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
RES_SET   <- c(250, 500, 1000)
BASE_RES  <- 50
DIST_RES  <- 250
RIVER_ORDER_LOCAL    <- 3
RIVER_ORDER_REGIONAL <- 5
SLOPE_CHECK_DEG <- c(15, 20, 30)
ANTHRO_DIST_M   <- 1000
TIE_DIST_M      <- 1

while (sink.number() > 0) sink()
source(REPO_CONFIG)
suppressPackageStartupMessages({ library(terra); library(sf) })
sf_use_s2(FALSE)

out <- dir_of("03_environment"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
env_dir <- file.path(REPO_DIR, "data_derived", "env")
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
getp <- function(l) file.path(RAW_DIR, man$path[man$layer == l])
crs_txt <- paste0("EPSG:", CRS_PROJ)

# ---- READ (before sink) ------------------------------------------------------
bnd <- st_transform(st_make_valid(st_read(getp("Boundary"), quiet = TRUE)), CRS_PROJ)
win <- st_union(st_geometry(bnd)); win_v <- vect(st_sf(geometry = win))
prep <- function(l) st_make_valid(st_transform(st_read(getp(l), quiet = TRUE), CRS_PROJ))
lith <- prep("Lithology"); geo <- prep("Geomorphology"); riv <- prep("Rivers")
pca  <- read.csv(file.path(dir_of("01_data_audit"), "polygon_class_areas.csv"), stringsAsFactors = FALSE)
for (r in RES_SET) stopifnot(file.exists(file.path(env_dir, paste0(r, "m"), "inwin.tif")),
                             file.exists(file.path(env_dir, paste0(r, "m"), "incenter.tif")),
                             file.exists(file.path(env_dir, paste0(r, "m"), "slope.tif")))
defs_f <- file.path(REPO_DIR, "docs", "covariate_definitions.csv")
defs <- if (file.exists(defs_f)) read.csv(defs_f, stringsAsFactors = FALSE) else NULL

bb <- st_bbox(win)
x0 <- floor(unname(bb["xmin"]) / 1000) * 1000; x1 <- ceiling(unname(bb["xmax"]) / 1000) * 1000
y0 <- floor(unname(bb["ymin"]) / 1000) * 1000; y1 <- ceiling(unname(bb["ymax"]) / 1000) * 1000
tmpl50 <- rast(xmin = x0, xmax = x1, ymin = y0, ymax = y1, resolution = BASE_RES, crs = crs_txt)
tmplD  <- rast(xmin = x0, xmax = x1, ymin = y0, ymax = y1, resolution = DIST_RES, crs = crs_txt)

logrows <- list()
lg <- function(step, detail) logrows[[length(logrows) + 1]] <<- data.frame(step = step, detail = detail, time = format(Sys.time()), stringsAsFactors = FALSE)
norm <- function(x) toupper(trimws(gsub("\\s+", " ", x)))
qs <- function(v) round(quantile(v, c(0, .05, .25, .5, .75, .95, 1), na.rm = TRUE), 3)
vect_of <- function(sfx, col) vect(st_sf(setNames(data.frame(as.integer(sfx[[col]])), col), geometry = st_geometry(sfx)))

sink(file.path(out, "categorical_hydrology_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== 0. SELF-TESTS (terra behaviour) ====\n")
  chk <- rast(nrows = 5, ncols = 5, xmin = 0, xmax = 5000, ymin = 0, ymax = 5000, crs = crs_txt)
  values(chk) <- NA; chk[13] <- 1
  dv <- values(distance(chk), mat = FALSE)
  cat("distance(): neighbour-cell distance =", round(dv[8], 2), "(expect 1000); corner =", round(dv[1], 1), "(expect 2828.4)\n")
  if (abs(dv[8] - 1000) > 1 || abs(dv[1] - 2828.43) > 5) stop("distance() self-test FAILED. Paste this output.")
  m <- rast(nrows = 4, ncols = 4, vals = c(1,1,2,2, 1,1,2,2, 3,3,3,3, 3,3,3,4))
  am <- values(aggregate(m, fact = 2, fun = "modal", na.rm = TRUE), mat = FALSE)
  cat("aggregate modal =", am, "(expect 1 2 3 3)\n")
  if (!identical(as.numeric(am), c(1, 2, 3, 3))) stop("aggregate(modal) self-test FAILED. Paste this output.")
  f3 <- rast(nrows = 3, ncols = 3, vals = c(1,1,1, 1,NA,2, 1,2,2))
  fm3 <- values(focal(f3, w = 3, fun = "modal", na.policy = "only", na.rm = TRUE), mat = FALSE)
  cat("focal modal fill, centre =", fm3[5], "(expect 1)\n")
  if (is.na(fm3[5]) || fm3[5] != 1) stop("focal(modal) self-test FAILED. Paste this output.")
  cat("All self-tests passed. terra", as.character(packageVersion("terra")), "\n")
  
  # ---- 1. LITHOLOGY, 5 classes ---------------------------------------------------
  cat("\n==== 1. LITHOLOGY: 5-CLASS COLLAPSE ====\n")
  LITH_MAP <- list(
    basalt             = "BASALT",
    crystalline_felsic = c("GRANITE GNEISS/MIGMATITE", "GNEISS/MIGMATITE", "GRANITE", "FOLIATED GRANITE",
                           "GREY HORNBLENDE BIOTITE GNEISS", "GRANULITE"),
    sandstone_clastics = c("SANDSTONE", "FINE GRAINED SANDSTONE", "FERRUGINOUS SANDSTONE", "SANDSTONE WITH COAL"),
    carbonate          = c("LIMESTONE", "DOLOMITIC LIMESTONE", "CHERTY LIMESTONE", "MARBLE",
                           "MANGANIFEROUS MARBLE", "CALC GNEISS", "CALC SILICATE ROCK"))
  LITH_NAMES <- c(names(LITH_MAP), "other")
  code_lith <- function(x) { n <- norm(x); o <- rep(5L, length(n))
  for (i in seq_along(LITH_MAP)) o[n %in% norm(LITH_MAP[[i]])] <- i; o }
  lith$lith5 <- code_lith(lith$lithologic)
  cl <- pca[pca$layer == "Lithology", ]
  cl$lith5 <- code_lith(cl$class)
  ref_area <- tapply(cl$a_km2, factor(LITH_NAMES[cl$lith5], levels = LITH_NAMES), sum)
  cat("Classes in window per polygon_class_areas.csv:", nrow(cl), "| classes in shapefile attribute:",
      length(unique(norm(lith$lithologic))), "\n")
  cat("Shapefile classes absent from polygon_class_areas.csv (outside window or zero area):",
      paste(setdiff(unique(norm(lith$lithologic)), norm(cl$class)), collapse = "; "), "\n")
  cat("\nMapping to 'other' (every class listed):\n")
  print(data.frame(class = cl$class[cl$lith5 == 5], a_km2 = cl$a_km2[cl$lith5 == 5]), row.names = FALSE)
  cat("\nReference areas from polygon_class_areas.csv (km2, % of window 21262.3):\n")
  print(data.frame(class = LITH_NAMES, a_km2 = round(as.numeric(ref_area), 1),
                   pct = round(100 * as.numeric(ref_area) / 21262.3, 2)), row.names = FALSE)
  lg("lith5", "collapsed lithologic field into 5 classes; unlisted classes -> other")
  r50_l <- rasterize(vect_of(lith, "lith5"), tmpl50, field = "lith5", touches = FALSE)
  
  # ---- 2. GEOMORPHOLOGY, 3 classes + recode ---------------------------------------
  cat("\n==== 2. GEOMORPHOLOGY: 3 CLASSES + NON-LANDFORM RECODE ====\n")
  G_MAP <- c("PEDIMENT PEDIPLAIN COMPLEX" = 1L,
             "MODERATELY DISSECTED PLATEAU" = 2L, "LOW DISSECTED PLATEAU" = 2L,
             "MODERATELY DISSECTED HILLS AND VALLEYS" = 2L, "LOW DISSECTED HILLS AND VALLEYS" = 2L,
             "HIGHLY DISSECTED HILLS AND VALLEYS" = 2L,
             "ALLUVIAL PLAIN" = 3L, "FLOOD PLAIN" = 3L)
  G_NON <- c("WATERBODY - RIVER", "DAM AND RESERVOIR", "WATERBODIES-OTHER", "QUARRY AND MINE DUMP", "ANTHROPOGENIC TERRAIN")
  G_NAMES <- c("pediplain", "dissected_upland", "alluvial")
  leg <- norm(geo$legend_sho)
  unm <- setdiff(unique(leg), c(names(G_MAP), G_NON))
  if (length(unm)) stop("Unmapped geomorphology legend values: ", paste(unm, collapse = "; "))
  geo$g3  <- unname(G_MAP[leg])
  geo$non <- leg %in% G_NON
  geo$area_km2 <- as.numeric(st_area(geo)) / 1e6
  cat("Polygons:", nrow(geo), "| landform:", sum(!geo$non), "| non-landform:", sum(geo$non),
      "| non-landform area (km2):", round(sum(geo$area_km2[geo$non]), 1), "\n")
  
  land <- geo[!geo$non, ]; non <- geo[geo$non, ]
  land_g <- st_geometry(land); non_g <- st_geometry(non)
  cand <- st_is_within_distance(non, land, dist = TIE_DIST_M)
  rule <- character(nrow(non)); newg <- integer(nrow(non))
  for (i in seq_len(nrow(non))) {
    if (i %% 500 == 0) cat("  recoding", i, "of", nrow(non), "\n")
    j <- cand[[i]]
    if (!length(j)) {
      k <- st_nearest_feature(non_g[i], land_g); newg[i] <- land$g3[k]; rule[i] <- "nearest_beyond_1m"
    } else if (length(unique(land$g3[j])) == 1) {
      newg[i] <- land$g3[j[1]]; rule[i] <- "touching_single_class"
    } else {
      res_i <- tryCatch({
        ring <- st_difference(st_buffer(non_g[i], TIE_DIST_M), non_g[i])
        a <- vapply(j, function(k) sum(as.numeric(st_area(suppressWarnings(st_intersection(land_g[k], ring))))), 0)
        tb <- tapply(a, land$g3[j], sum); c(as.integer(names(which.max(tb))), 1L)
      }, error = function(e) c(NA_integer_, 0L))
      if (is.na(res_i[1])) { k <- st_nearest_feature(non_g[i], land_g); newg[i] <- land$g3[k]; rule[i] <- "tie_fallback_nearest" }
      else { newg[i] <- res_i[1]; rule[i] <- "touching_multi_boundary_length" }
    }
  }
  non$g3 <- newg; non$rule <- rule
  rt <- aggregate(area_km2 ~ rule, data = data.frame(rule = non$rule, area_km2 = non$area_km2), FUN = sum)
  rt$n_polygons <- as.integer(table(non$rule)[rt$rule]); rt$area_km2 <- round(rt$area_km2, 1)
  cat("\nRecode rule used (non-landform polygons):\n"); print(rt, row.names = FALSE)
  rec <- aggregate(area_km2 ~ legend_sho + g3, data = data.frame(legend_sho = non$legend_sho, g3 = G_NAMES[non$g3], area_km2 = non$area_km2), FUN = sum)
  rec$area_km2 <- round(rec$area_km2, 1)
  cat("\nNon-landform legend class -> assigned landform class (km2):\n"); print(rec, row.names = FALSE)
  write.csv(rec, file.path(out, "geomorph_recode_table.csv"), row.names = FALSE)
  cat("\nGeomorphology area before recode (km2): landform classes only\n")
  bef <- tapply(land$area_km2, factor(G_NAMES[land$g3], levels = G_NAMES), sum)
  aft <- tapply(c(land$area_km2, non$area_km2), factor(G_NAMES[c(land$g3, non$g3)], levels = G_NAMES), sum)
  print(data.frame(class = G_NAMES, before_km2 = round(as.numeric(bef), 1), after_km2 = round(as.numeric(aft), 1)), row.names = FALSE)
  
  geo_all <- rbind(data.frame(g3 = land$g3), data.frame(g3 = non$g3))
  geo_sf  <- st_sf(g3 = geo_all$g3, geometry = c(st_geometry(land), st_geometry(non)))
  r50_g <- rasterize(vect_of(geo_sf, "g3"), tmpl50, field = "g3", touches = FALSE)
  lg("geomorph3", "3 landform classes; non-landform polygons recoded by touching/boundary-length/nearest rule")
  
  # ---- 3. RIVERS: distance to order >= 3 and >= 5 -------------------------------
  cat("\n==== 3. RIVER DISTANCE LAYERS ====\n")
  cat("River features:", nrow(riv), "| Strahler order range:", range(riv$ordsh, na.rm = TRUE), "| order NA:", sum(is.na(riv$ordsh)), "\n")
  mk_dist <- function(k) {
    sel <- riv[!is.na(riv$ordsh) & riv$ordsh >= k, ]
    rr  <- rasterize(vect(st_sf(geometry = st_geometry(sel))), tmplD, field = 1, touches = TRUE)
    cat(sprintf("order >= %d: %d segments, length %.0f km, line cells at %d m: %d\n", k, nrow(sel),
                sum(as.numeric(st_length(sel))) / 1000, DIST_RES, global(!is.na(rr), "sum")[1, 1]))
    distance(rr) / 1000
  }
  dkm <- list(riv3 = mk_dist(RIVER_ORDER_LOCAL), riv5 = mk_dist(RIVER_ORDER_REGIONAL))
  cat("Note: the Rivers layer holds no features outside the window, so distances near the border\n",
      "are overestimated where the nearest river lies beyond the boundary.\n")
  lg("rivers", sprintf("distance at %d m from rasterised lines, order>=%d and >=%d, km", DIST_RES, RIVER_ORDER_LOCAL, RIVER_ORDER_REGIONAL))
  
  # ---- 4. WRITE LAYERS PER RESOLUTION ---------------------------------------------
  fill_cat <- function(r, inw, max_iter = 4) {
    n0 <- global(is.na(r) & !is.na(inw), "sum", na.rm = TRUE)[1, 1]
    for (k in seq_len(max_iter)) {
      if (global(is.na(r) & !is.na(inw), "sum", na.rm = TRUE)[1, 1] == 0) break
      r <- focal(r, w = 3, fun = "modal", na.policy = "only", na.rm = TRUE)
    }
    list(r = r, n_before = n0, n_after = global(is.na(r) & !is.na(inw), "sum", na.rm = TRUE)[1, 1])
  }
  cls_tab <- list(); dst_tab <- list()
  cat("\n==== 4. LAYERS PER RESOLUTION (window-centre cells) ====\n")
  for (res in RES_SET) {
    fct <- res / BASE_RES; d <- file.path(env_dir, paste0(res, "m"))
    inw <- rast(file.path(d, "inwin.tif")); ic <- !is.na(values(rast(file.path(d, "incenter.tif")), mat = FALSE))
    cat(sprintf("\n-- %d m grid: centre cells %d\n", res, sum(ic)))
    for (nm in c("lith5", "geomorph3")) {
      src <- if (nm == "lith5") r50_l else r50_g; labs <- if (nm == "lith5") LITH_NAMES else G_NAMES
      a  <- aggregate(src, fact = fct, fun = "modal", na.rm = TRUE)
      if (!compareGeom(a, inw, stopOnError = FALSE)) stop("Grid mismatch for ", nm, " at ", res, " m")
      fl <- fill_cat(a, inw)
      r  <- mask(fl$r, inw); names(r) <- nm
      writeRaster(r, file.path(d, paste0(nm, ".tif")), overwrite = TRUE, datatype = "INT1U")
      v  <- values(r, mat = FALSE)[ic]
      tb <- table(factor(v, levels = seq_along(labs)))
      cat(sprintf("%-9s touching-NA cells filled: %d -> %d remaining | NA in centre cells: %d\n",
                  nm, fl$n_before, fl$n_after, sum(is.na(v))))
      print(data.frame(class = labs, n_cells = as.integer(tb), area_km2 = round(as.integer(tb) * res^2 / 1e6, 1),
                       pct = round(100 * as.integer(tb) / sum(tb), 2)), row.names = FALSE)
      cls_tab[[length(cls_tab) + 1]] <- data.frame(res_m = res, layer = nm, class = labs, n_cells = as.integer(tb),
                                                   area_km2 = round(as.integer(tb) * res^2 / 1e6, 1), pct = round(100 * as.integer(tb) / sum(tb), 2),
                                                   NA_filled = fl$n_before, NA_remaining = fl$n_after)
    }
    for (k in names(dkm)) {
      dd <- if (res > DIST_RES) aggregate(dkm[[k]], fact = res / DIST_RES, fun = "mean", na.rm = TRUE) else dkm[[k]]
      if (!compareGeom(dd, inw, stopOnError = FALSE)) stop("Grid mismatch for distance at ", res, " m")
      dd <- mask(dd, inw); names(dd) <- paste0("dist_", k, "_km")
      ld <- log(1 + dd); names(ld) <- paste0("ldist_", k)
      writeRaster(dd, file.path(d, paste0("dist_", k, "_km.tif")), overwrite = TRUE, datatype = "FLT4S")
      writeRaster(ld, file.path(d, paste0("ldist_", k, ".tif")), overwrite = TRUE, datatype = "FLT4S")
      v <- values(dd, mat = FALSE)[ic]; q <- qs(v)
      cat(sprintf("dist_%s_km: NA %d | min %.3f | p5 %.3f | median %.3f | p95 %.3f | max %.3f\n",
                  k, sum(is.na(v)), q["0%"], q["5%"], q["50%"], q["95%"], q["100%"]))
      dst_tab[[length(dst_tab) + 1]] <- data.frame(res_m = res, layer = names(dd), NA_centre = sum(is.na(v)),
                                                   min = q["0%"], p5 = q["5%"], p25 = q["25%"], median = q["50%"], p75 = q["75%"], p95 = q["95%"], max = q["100%"], row.names = NULL)
    }
    lg("write", sprintf("%d m: lith5, geomorph3, dist/ldist riv3, riv5 written; touching-NA categorical cells filled by modal focal", res))
  }
  write.csv(do.call(rbind, cls_tab), file.path(out, "class_summary.csv"), row.names = FALSE)
  write.csv(do.call(rbind, dst_tab), file.path(out, "distance_summary.csv"), row.names = FALSE)
  
  # ---- 5. EXTREME-SLOPE CHECK (external layers only) ----------------------------
  cat("\n==== 5. EXTREME-SLOPE CHECK, 500 m GRID ====\n")
  d5 <- file.path(env_dir, "500m")
  sl <- rast(file.path(d5, "slope.tif")); el <- rast(file.path(d5, "elev.tif"))
  ic5 <- !is.na(values(rast(file.path(d5, "incenter.tif")), mat = FALSE))
  sv <- values(sl, mat = FALSE)
  anth_leg <- c("QUARRY AND MINE DUMP", "DAM AND RESERVOIR", "ANTHROPOGENIC TERRAIN")
  anth <- geo[leg %in% anth_leg, ]
  ar <- rasterize(vect(st_sf(geometry = st_geometry(anth))), sl, field = 1, touches = TRUE)
  da <- values(distance(ar), mat = FALSE)
  base_near <- mean(da[ic5] <= ANTHRO_DIST_M, na.rm = TRUE)
  cat(sprintf("Anthropogenic polygons (quarry/mine dump, dam/reservoir, anthropogenic terrain): %d. Share of ALL window cells within %d m of one: %.1f%%\n",
              nrow(anth), ANTHRO_DIST_M, 100 * base_near))
  sx <- do.call(rbind, lapply(SLOPE_CHECK_DEG, function(th) {
    h <- ic5 & !is.na(sv) & sv > th
    data.frame(slope_gt_deg = th, n_cells = sum(h), pct_of_cells = round(100 * mean(h[ic5]), 3),
               pct_within_1km_of_anthropogenic = round(100 * mean(da[h] <= ANTHRO_DIST_M, na.rm = TRUE), 1))
  }))
  print(sx, row.names = FALSE)
  write.csv(sx, file.path(out, "slope_extreme_check.csv"), row.names = FALSE)
  top <- order(sv, decreasing = TRUE, na.last = NA); top <- top[ic5[top]][1:10]
  xy <- xyFromCell(sl, top)
  ll <- st_coordinates(st_transform(st_as_sf(data.frame(x = xy[, 1], y = xy[, 2]), coords = c("x", "y"), crs = CRS_PROJ), 4326))
  cat("\nTen steepest 500 m window cells (inspect in GIS):\n")
  print(data.frame(x_m = round(xy[, 1]), y_m = round(xy[, 2]), lon = round(ll[, 1], 4), lat = round(ll[, 2], 4),
                   slope_deg = round(sv[top], 1), elev_m = round(values(el, mat = FALSE)[top]),
                   dist_to_anthropogenic_m = round(da[top])), row.names = FALSE)
  lg("slope_check", "extreme slope cells compared with quarry/mine/dam/anthropogenic polygons (external layers only)")
  
  # ---- 6. DEFINITIONS UPDATE ---------------------------------------------------------
  if (!is.null(defs)) {
    defs$transform[defs$covariate == "slope"] <- "log(1+slope_deg), z-score"
    defs$status <- "APPROVED (author-delegated 2026-11-06); not FROZEN"
    write.csv(defs, defs_f, row.names = FALSE)
    cat("\ncovariate_definitions.csv: slope transform set to log(1+slope_deg), z-score\n")
  }
  write.csv(do.call(rbind, logrows), file.path(out, "processing_log_10.csv"), row.names = FALSE)
  cat("\nSaved layers to data_derived/env/<res>m/ ; tables in 03_environment/outputs\n")
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_10.txt"))
cat("\nDONE. Outputs in:", out, "\n")