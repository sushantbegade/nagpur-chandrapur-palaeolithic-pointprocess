# =============================================================================
# 08_window_variants_benchmark.R
# Stage: 02_spatial_window | Gate 3 support: computational window choice
# Purpose: compare exact vs simplified (polygon) vs pixel-mask windows:
#   geometry change, sites outside, K(r) stability, and run time.
# Uses ONLY the point pattern and window. No covariates, no models.
# Outputs (02_spatial_window/outputs/): window_variants.csv, window_variants_report.txt,
#   sessionInfo_08.txt
# Not executed by assistant. Run top to bottom in RStudio.
# Expect a few minutes. Progress prints per variant. If the exact-window
# translation benchmark stalls >10 min, press Stop and set TRANSLATE_N <- 0.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
TOLERANCES  <- c(50, 100, 250, 500)      # polygon simplification tolerances (m)
MASK_EPS    <- c(250, 500, 1000)          # pixel sizes for mask windows (m)
TRANSLATE_N <- 15                         # points used for translation-correction timing (0 = skip)

while (sink.number() > 0) sink()
source(REPO_CONFIG)
need <- c("spatstat.geom", "spatstat.explore", "spatstat.random", "sf")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: install.packages(c(",
                       paste(sprintf('"%s"', miss), collapse = ", "), "))")
suppressPackageStartupMessages({ library(sf); library(spatstat.geom)
  library(spatstat.explore); library(spatstat.random) })
sf_use_s2(FALSE)

out <- dir_of("02_spatial_window"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
man <- read.csv(file.path(REPO_DIR, "docs", "input_manifest.csv"), stringsAsFactors = FALSE)
getp <- function(l) file.path(RAW_DIR, man$path[man$layer == l])

# ---- READ (before sink) ------------------------------------------------------
bnd <- st_transform(st_make_valid(st_read(getp("Boundary"), quiet = TRUE)), CRS_PROJ)
win_sf <- st_union(st_geometry(bnd))
obj <- readRDS(file.path(REPO_DIR, "data_derived", "ppp_NA_window.rds"))
X0 <- obj$X_m; W0 <- obj$W_m
xy <- cbind(X0$x, X0$y); n <- nrow(xy)
A0 <- sum(as.numeric(st_area(win_sf)))
lam_m2 <- n / A0
rgrid <- seq(0, 15000, by = 250)
set.seed(MASTER_SEED)
sub_idx <- if (TRANSLATE_N > 0) sample(n, TRANSLATE_N) else integer(0)

nvert_owin <- function(W) if (W$type == "polygonal") sum(vapply(W$bdry, function(b) length(b$x), 1L)) else NA
tm <- function(expr) tryCatch(unname(system.time(expr)["elapsed"]), error = function(e) NA_real_)

sink(file.path(out, "window_variants_report.txt"), split = TRUE)
tryCatch({
  cat("==== WINDOW VARIANTS: geometry, K stability, timing ====\n")
  cat("n =", n, "| exact area (km2):", round(A0 / 1e6, 1), "| subset for translation timing:", TRANSLATE_N, "\n")
  cat("spatstat.geom", as.character(packageVersion("spatstat.geom")),
      "| spatstat.explore", as.character(packageVersion("spatstat.explore")), "\n\n")
  
  variants <- list(list(name = "exact", type = "polygon", par = 0, sfc = win_sf, W = W0))
  for (tl in TOLERANCES) {
    g <- tryCatch(st_union(st_make_valid(st_simplify(win_sf, preserveTopology = TRUE, dTolerance = tl))),
                  error = function(e) NULL)
    if (is.null(g)) { cat("simplify failed at", tl, "\n"); next }
    Wv <- tryCatch(as.owin(g), error = function(e) NULL)
    if (!is.null(Wv)) { unitname(Wv) <- c("metre", "metres")
    variants[[length(variants) + 1]] <- list(name = paste0("simplify_", tl, "m"), type = "polygon",
                                             par = tl, sfc = g, W = Wv) }
  }
  for (e in MASK_EPS) {
    Wv <- tryCatch(as.mask(W0, eps = e), error = function(er) NULL)
    if (!is.null(Wv)) variants[[length(variants) + 1]] <- list(name = paste0("mask_", e, "m"),
                                                               type = "mask", par = e, sfc = NULL, W = Wv)
  }
  
  K_ref <- NULL; rows <- list()
  for (v in variants) {
    cat("---", v$name, "---\n"); flush.console()
    Wv <- v$W
    ins <- inside.owin(xy[, 1], xy[, 2], Wv)
    Av  <- area.owin(Wv)
    symd <- if (!is.null(v$sfc) && v$par > 0)
      tryCatch(sum(as.numeric(st_area(st_sym_difference(win_sf, v$sfc)))) / 1e6, error = function(e) NA) else
        if (v$par == 0) 0 else NA
    cat("vertices/pixels:", if (v$type == "polygon") nvert_owin(Wv) else prod(Wv$dim),
        "| area (km2):", round(Av / 1e6, 1), "| sites outside:", sum(!ins), "\n")
    if (any(!ins)) cat("Outside IDs:", obj$marks$site_id[!ins], "\n")
    
    Xv <- ppp(xy[ins, 1], xy[ins, 2], window = Wv); unitname(Xv) <- c("metre", "metres")
    t_iso <- tm(Kiso <- Kest(Xv, r = rgrid, correction = "isotropic"))
    t_bor <- tm(Kbor <- Kest(Xv, r = rgrid, correction = "border"))
    t_sim <- tm(rpoispp(lam_m2, win = Wv, nsim = 20))
    t_trn <- NA_real_; t_trn_full <- NA_real_
    if (TRANSLATE_N > 0) {
      Xs <- ppp(xy[sub_idx[ins[sub_idx]], 1], xy[sub_idx[ins[sub_idx]], 2], window = Wv)
      t_trn <- tm(Kest(Xs, r = rgrid, correction = "translate"))
      t_trn_full <- t_trn * (n / npoints(Xs))^2 / 60
    }
    dK <- NA_real_
    if (v$name == "exact") K_ref <- Kiso
    if (!is.null(K_ref) && !inherits(Kiso, "try-error") && exists("Kiso")) {
      k  <- Kiso$iso; kr <- K_ref$iso; sel <- Kiso$r >= 1000 & is.finite(k) & is.finite(kr) & kr > 0
      dK <- if (any(sel)) max(abs(k[sel] - kr[sel]) / kr[sel]) else NA_real_
    }
    cat(sprintf("time (s): K iso %.1f | K border %.1f | 20 rpoispp %.1f | translate on %d pts %.1f -> est. full-n %.1f min\n",
                t_iso, t_bor, t_sim, length(sub_idx), t_trn, t_trn_full))
    cat("max relative |dK_iso| vs exact, r>=1 km:", signif(dK, 3), "\n\n")
    rows[[length(rows) + 1]] <- data.frame(
      variant = v$name, type = v$type, param_m = v$par,
      vertices_or_pixels = if (v$type == "polygon") nvert_owin(Wv) else prod(Wv$dim),
      area_km2 = round(Av / 1e6, 1), area_diff_pct = round(100 * (Av - A0) / A0, 3),
      symdiff_km2 = round(symd, 2), sites_outside = sum(!ins),
      t_Kiso_s = round(t_iso, 2), t_Kborder_s = round(t_bor, 2), t_rpoispp20_s = round(t_sim, 2),
      t_translate_subset_s = round(t_trn, 2), est_translate_full_min = round(t_trn_full, 1),
      max_rel_dK_iso_ge1km = signif(dK, 3), stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, rows)
  write.csv(res, file.path(out, "window_variants.csv"), row.names = FALSE)
  cat("==== SUMMARY TABLE ====\n"); options(width = 230); print(res, row.names = FALSE)
  cat("\nNo window is frozen by this script. The choice is made after reading this table.\n")
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_08.txt"))
cat("\nDONE. Outputs in:", out, "\n")