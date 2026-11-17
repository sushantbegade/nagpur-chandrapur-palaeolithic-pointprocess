# =============================================================================
# 12_download_observation_layers.R
# Stage: 04_observation | Gate 5 (observation) support
# Purpose: download ESA WorldCover 10 m 2021 v200 tiles covering the window and
#   OpenStreetMap roads (motorway..tertiary) clipped to window + margin; verify
#   files; log provenance (source, date, bytes, MD5, licence); extend the input manifest.
# Does NOT read site data. No models. Writes raw downloads to RAW_DIR (outside repo).
# Outputs:
#   <RAW_DIR>/Observation Layers (External)/ESA_WorldCover_10m_2021_v200_<tile>_Map.tif
#   <RAW_DIR>/Observation Layers (External)/roads_osm_motorway_to_tertiary.gpkg
#   docs/external_download_log.csv  (created / updated by layer name)
#   docs/input_manifest.csv         (rows appended if missing)
#   04_observation/outputs/download_report.txt, sessionInfo_12.txt
# Not executed by assistant. Run top to bottom in RStudio. Needs internet.
# Large downloads; check free disk space first. Do not re-run the Overpass part
# repeatedly (public servers are rate-limited): existing files are skipped.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
MARGIN_DEG       <- 0.2
WC_BASE          <- "https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
HIGHWAY_CLASSES  <- c("motorway", "trunk", "primary", "secondary", "tertiary")
OSM_CHUNKS       <- c(2, 2)          # columns x rows
OVERPASS_TIMEOUT <- 300
DL_SUBDIR        <- "Observation Layers (External)"

while (sink.number() > 0) sink()
source(REPO_CONFIG)
suppressPackageStartupMessages({ library(sf); library(terra) })
sf_use_s2(FALSE)

out <- dir_of("04_observation"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
dl_dir <- file.path(RAW_DIR, DL_SUBDIR); dir.create(dl_dir, showWarnings = FALSE, recursive = TRUE)
man_f <- file.path(REPO_DIR, "docs", "input_manifest.csv")
man   <- read.csv(man_f, stringsAsFactors = FALSE)
getp  <- function(l) file.path(RAW_DIR, man$path[man$layer == l])
log_f <- file.path(REPO_DIR, "docs", "external_download_log.csv")

# ---- READ (before sink) ------------------------------------------------------
bnd <- st_transform(st_make_valid(st_read(getp("Boundary"), quiet = TRUE)), CRS_PROJ)
win <- st_union(st_geometry(bnd))
bb  <- st_bbox(st_transform(win, 4326))
xmin <- unname(bb["xmin"]) - MARGIN_DEG; xmax <- unname(bb["xmax"]) + MARGIN_DEG
ymin <- unname(bb["ymin"]) - MARGIN_DEG; ymax <- unname(bb["ymax"]) + MARGIN_DEG
lon0 <- floor(xmin / 3) * 3; lon1 <- floor(xmax / 3) * 3
lat0 <- floor(ymin / 3) * 3; lat1 <- floor(ymax / 3) * 3
tiles <- expand.grid(lat = seq(lat0, lat1, by = 3), lon = seq(lon0, lon1, by = 3))
stopifnot(all(tiles$lat >= 0), all(tiles$lon >= 0))
tiles$id   <- sprintf("N%02dE%03d", tiles$lat, tiles$lon)
tiles$file <- sprintf("ESA_WorldCover_10m_2021_v200_%s_Map.tif", tiles$id)
tiles$url  <- paste0(WC_BASE, tiles$file)
tiles$dest <- file.path(dl_dir, tiles$file)

dl_file <- function(url, dest, tries = 3) {
  if (file.exists(dest) && file.size(dest) > 1e6) return("exists")
  options(timeout = 7200)
  for (k in seq_len(tries)) {
    ok <- tryCatch({ download.file(url, dest, mode = "wb", quiet = FALSE, method = "libcurl"); TRUE },
                   error   = function(e) { cat("download error:", conditionMessage(e), "\n"); FALSE },
                   warning = function(w) { cat("download warning:", conditionMessage(w), "\n"); FALSE })
    if (ok && file.exists(dest) && file.size(dest) > 1e6) return("downloaded")
    if (file.exists(dest)) file.remove(dest)
    Sys.sleep(10)
  }
  "failed"
}
logrows <- list()
add_log <- function(layer, url, path, bytes, md5, licence, attribution, note)
  logrows[[length(logrows) + 1]] <<- data.frame(layer = layer, source_url = url, local_path = path,
                                                downloaded = format(Sys.time()), bytes = bytes, md5 = md5, licence = licence,
                                                attribution = attribution, note = note, stringsAsFactors = FALSE)

sink(file.path(out, "download_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== 0. PLAN ====\n")
  cat("Window + margin (lon/lat):", round(xmin, 3), round(ymin, 3), round(xmax, 3), round(ymax, 3), "\n")
  cat("WorldCover tiles needed:", paste(tiles$id, collapse = ", "), "\n")
  cat("Download folder:", dl_dir, "\n")
  
  # ---- 1. WORLDCOVER -----------------------------------------------------------
  cat("\n==== 1. ESA WORLDCOVER 10 m 2021 v200 ====\n")
  WC_CLASSES <- c("10" = "tree cover", "20" = "shrubland", "30" = "grassland", "40" = "cropland", "50" = "built-up",
                  "60" = "bare/sparse", "70" = "snow/ice", "80" = "permanent water", "90" = "herbaceous wetland",
                  "95" = "mangroves", "100" = "moss/lichen")
  wc_ok <- character(0)
  for (i in seq_len(nrow(tiles))) {
    cat("\n--", tiles$id[i], "--\n")
    st <- dl_file(tiles$url[i], tiles$dest[i]); cat("status:", st, "\n")
    if (st == "failed") { cat("FAILED:", tiles$url[i], "\n"); next }
    r <- tryCatch(rast(tiles$dest[i]), error = function(e) { cat("rast() failed:", conditionMessage(e), "\n"); NULL })
    if (is.null(r)) next
    cd <- tryCatch(crs(r, describe = TRUE), error = function(e) data.frame(code = NA, name = NA))
    cat("size:", ncol(r), "x", nrow(r), "| res (deg):", signif(res(r)[1], 6), "| CRS EPSG:", cd$code, "\n")
    cat("extent (lon/lat):", round(as.vector(ext(r)), 3), "\n")
    sm <- spatSample(r[[1]], 400000, method = "regular", warn = FALSE)[, 1]
    tb <- table(factor(sm, levels = as.integer(names(WC_CLASSES))))
    cat("class sample (whole tile, incl. outside window):\n")
    print(data.frame(code = names(WC_CLASSES), class = unname(WC_CLASSES), pct = round(100 * as.numeric(tb) / sum(tb), 2)), row.names = FALSE)
    md <- unname(tools::md5sum(tiles$dest[i]))
    add_log(paste0("WorldCover_", tiles$id[i]), tiles$url[i], tiles$dest[i], file.size(tiles$dest[i]), md,
            "CC BY 4.0", "© ESA WorldCover project / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium",
            "ESA WorldCover 10 m 2021 v200; DOI 10.5281/zenodo.7254221; copy full author list from Zenodo record")
    wc_ok <- c(wc_ok, i)
  }
  
  # ---- 2. OSM ROADS ------------------------------------------------------------
  cat("\n==== 2. OPENSTREETMAP ROADS ====\n")
  roads_f <- file.path(dl_dir, "roads_osm_motorway_to_tertiary.gpkg")
  if (file.exists(roads_f)) {
    cat("Roads file exists, skipping download:", roads_f, "\n")
  } else if (!requireNamespace("osmdata", quietly = TRUE)) {
    cat("Package osmdata missing. Run install.packages('osmdata') and re-run. Roads skipped.\n")
  } else {
    suppressPackageStartupMessages(library(osmdata))
    xs <- seq(xmin, xmax, length.out = OSM_CHUNKS[1] + 1); ys <- seq(ymin, ymax, length.out = OSM_CHUNKS[2] + 1)
    parts <- list(); stamps <- character(0); failed <- 0
    for (i in seq_len(OSM_CHUNKS[1])) for (j in seq_len(OSM_CHUNKS[2])) {
      bbm <- matrix(c(xs[i], xs[i + 1], ys[j], ys[j + 1]), nrow = 2, byrow = TRUE,
                    dimnames = list(c("x", "y"), c("min", "max")))
      q <- add_osm_feature(opq(bbox = bbm, timeout = OVERPASS_TIMEOUT), key = "highway", value = HIGHWAY_CLASSES)
      qs <- opq_string(q)
      nums <- as.numeric(strsplit(gsub("[()]", "", regmatches(qs, regexpr("\\(-?[0-9.]+,-?[0-9.]+,-?[0-9.]+,-?[0-9.]+\\)", qs))), ",")[[1]])
      cat(sprintf("chunk %d,%d  query box (south,west,north,east) = %s\n", i, j, paste(round(nums, 3), collapse = ", ")))
      if (length(nums) != 4 || abs(nums[1] - ys[j]) > 0.01 || abs(nums[2] - xs[i]) > 0.01)
        stop("Overpass bbox axis order is not (south,west,north,east). Nothing was downloaded. Paste this output.")
      res <- NULL
      for (att in 1:3) {
        res <- tryCatch(osmdata_sf(q, quiet = FALSE), error = function(e) { cat("Overpass error, attempt", att, ":", conditionMessage(e), "\n"); NULL })
        if (!is.null(res)) break; Sys.sleep(30)
      }
      if (is.null(res)) { failed <- failed + 1; next }
      ln <- res$osm_lines
      if (!is.null(ln) && nrow(ln)) { keep <- intersect(c("osm_id", "highway", "name", "ref"), names(ln)); parts[[length(parts) + 1]] <- ln[, keep] }
      stamps <- c(stamps, tryCatch(as.character(res$meta$timestamp), error = function(e) NA_character_))
      cat("  ways in chunk:", if (is.null(ln)) 0 else nrow(ln), "\n")
    }
    if (failed > 0 || !length(parts)) {
      cat("Overpass chunks failed:", failed, "of", prod(OSM_CHUNKS), ". Roads NOT saved (incomplete data would bias the proxy). Re-run later or paste this output.\n")
    } else {
      cols <- Reduce(intersect, lapply(parts, names))
      roads <- do.call(rbind, lapply(parts, function(p) p[, cols]))
      roads <- roads[!duplicated(roads$osm_id), ]
      st_write(roads, roads_f, delete_dsn = TRUE, quiet = TRUE)
      cat("Saved:", roads_f, "| unique ways:", nrow(roads), "\n")
      add_log("Roads_OSM", "Overpass API via osmdata", roads_f, file.size(roads_f), unname(tools::md5sum(roads_f)),
              "ODbL 1.0", "© OpenStreetMap contributors",
              paste0("highway in {", paste(HIGHWAY_CLASSES, collapse = ","), "}; Overpass timestamp(s): ",
                     paste(unique(stamps[!is.na(stamps)]), collapse = " | "), "; query date ", format(Sys.Date())))
    }
  }
  if (file.exists(roads_f)) {
    rd <- tryCatch(st_read(roads_f, quiet = TRUE), error = function(e) NULL)
    if (!is.null(rd)) {
      ru <- st_transform(rd, CRS_PROJ)
      ins <- lengths(st_intersects(ru, win)) > 0
      clip <- tryCatch(suppressWarnings(st_intersection(ru[ins, ], win)), error = function(e) NULL)
      cat("\nRoads in file:", nrow(rd), "| ways touching window:", sum(ins), "\n")
      if (!is.null(clip)) {
        clip$len_km <- as.numeric(st_length(clip)) / 1000
        cat("Length inside window by class (km):\n"); print(round(tapply(clip$len_km, clip$highway, sum)))
        cat("Total length inside window (km):", round(sum(clip$len_km)), "| density (km/km2):",
            round(sum(clip$len_km) / (sum(as.numeric(st_area(win))) / 1e6), 4), "\n")
      }
    }
  }
  
  # ---- 3. LOG + MANIFEST --------------------------------------------------------
  if (length(logrows)) {
    new <- do.call(rbind, logrows)
    old <- if (file.exists(log_f)) read.csv(log_f, stringsAsFactors = FALSE) else NULL
    if (!is.null(old)) old <- old[!old$layer %in% new$layer, ]
    write.csv(rbind(old, new), log_f, row.names = FALSE)
    cat("\nLogged", nrow(new), "layer(s) to docs/external_download_log.csv\n")
  }
  add <- data.frame(stringsAsFactors = FALSE, type = character(), layer = character(), domain = character(), path = character())
  for (i in seq_len(nrow(tiles))) if (file.exists(tiles$dest[i]))
    add <- rbind(add, data.frame(type = "raster", layer = paste0("WorldCover_", tiles$id[i]), domain = "observation_external",
                                 path = file.path(DL_SUBDIR, tiles$file[i])))
  if (file.exists(roads_f)) add <- rbind(add, data.frame(type = "vector", layer = "Roads_OSM", domain = "observation_external",
                                                         path = file.path(DL_SUBDIR, basename(roads_f))))
  add <- add[!add$layer %in% man$layer, ]
  if (nrow(add)) { write.csv(rbind(man[, c("type", "layer", "domain", "path")], add), man_f, row.names = FALSE)
    cat("Manifest extended with:", paste(add$layer, collapse = ", "), "\n") } else cat("Manifest unchanged.\n")
  
  cat("\nAttribution to carry into the paper:\n")
  cat(" WorldCover: \u00a9 ESA WorldCover project / Contains modified Copernicus Sentinel data (2021) processed by ESA WorldCover consortium\n")
  cat(" Roads: \u00a9 OpenStreetMap contributors (ODbL). Share-alike applies to derived databases. Check before archiving.\n")
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_12.txt"))
cat("\nDONE. Outputs in:", out, "\n")