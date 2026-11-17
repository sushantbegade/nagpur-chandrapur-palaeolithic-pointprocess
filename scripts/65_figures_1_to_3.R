# =====================================================================
# 65_figures_1_to_3.R
# Stage   : 14_figures
# Folder  : 14_figures/scripts/   (NEW; creates 14_figures/outputs)
# Purpose : Figure 1 (study region, window, sites by coordinate source), Figure 2 (database
#           reconciliation flow, values read from T1a_audit.csv), Figure 3 (environmental covariates
#           and the fixed observation offset O2; no sites drawn: circularity). PNG (300 dpi) + PDF.
# Outputs : outputs/Fig01_study_region.*, Fig02_database_reconciliation.*, Fig03_covariates_observation.*,
#           65_report.txt, 65_figure_register.csv, 65_sessionInfo.txt, protocol_deviation_log_before_65.csv;
#           docs/protocol_deviation_log.csv (1 row)
# Notes   : No analysis is run; figures only. Colours: Okabe-Ito (colour-blind safe). ASCII labels only.
#           Each figure is built inside tryCatch; failures are reported, not hidden.
#           Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)
CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))

OUT <- dir_of("14_figures", "outputs"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
D15 <- dir_of("15_tables", "outputs")

# ---- open inputs BEFORE sink ----
f_t1 <- file.path(D15, "T1a_audit.csv"); stopifnot(file.exists(f_t1)); T1a <- read.csv(f_t1, stringsAsFactors = FALSE)
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
SITES <- read.csv(f_sites, stringsAsFactors = FALSE)
stopifnot(nrow(SITES) == 197, max(abs(X_obs_km$x - SITES$x_m / 1000), abs(X_obs_km$y - SITES$y_m / 1000)) < 1e-3)
f_obs <- file.path(REPO_DIR, "data_derived", "recovery_obs_cache.rds"); stopifnot(file.exists(f_obs)); OBS <- readRDS(f_obs)
f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE); LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")

ROW <- data.frame(date = "2026-11-15",
                  item = "Clarification: figure list - Figure 7 (Geyer diagnostics) replaced by the Stage 1 classification figure; supplementary figures added",
                  reason = "The hierarchy ended at M3 (design Class B) and Stage 1 showed Gibbs parameters not identifiable (Gate 6 FAIL), so M4 was not fitted and Geyer generative diagnostics do not exist. Figure 7 shows the Stage 1 replicate classification per recovery cell instead. Supplementary figures added: ERL p versus nsim and self-calibration (scripts 51-52), source structure of close pairs (scripts 61-62), block validation (scripts 56-57).",
                  effect_on_inference = "Presentation only; no analysis changes. Figure content follows the tables (15_tables), which are generated from the saved outputs.",
                  approved_by = "author-delegated (assistant decision; clarification, no frozen item changed)", stringsAsFactors = FALSE)
if (identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT, "protocol_deviation_log_before_65.csv"), overwrite = TRUE)
  if (!(ROW$item %in% LOG0$item)) write.csv(rbind(LOG0, ROW), f_log, row.names = FALSE); LOG_MSG <- "Log row appended (backup in outputs)."
} else { write.csv(ROW, file.path(OUT, "65_log_row_for_manual_append.csv"), row.names = FALSE); LOG_MSG <- "Log columns differ; row saved for manual append." }

# ---- helpers ----
OI <- c(black = "#000000", orange = "#E69F00", sky = "#56B4E9", green = "#009E73", yellow = "#F0E442", blue = "#0072B2", verm = "#D55E00", pink = "#CC79A7", grey = "#BBBBBB")
FR <- Frame(W_km); XR <- FR$xrange; YR <- FR$yrange
save_fig <- function(name, w_mm, h_mm, fun) {
  for (dev in c("png", "pdf")) {
    f <- file.path(OUT, paste0(name, ".", dev))
    if (dev == "png") png(f, width = round(w_mm / 25.4 * 300), height = round(h_mm / 25.4 * 300), res = 300) else pdf(f, width = w_mm / 25.4, height = h_mm / 25.4)
    tryCatch(fun(), finally = dev.off())
  }
  invisible(NULL)
}
outline <- function(lwd = 0.4, col = "grey20", fill = NA) for (b in W_km$bdry) polygon(b$x, b$y, col = fill, border = col, lwd = lwd)
tn <- function(item) { v <- T1a$n[T1a$item == item]; if (length(v) != 1) stop("T1a item not found: ", item); v }

REG <- data.frame()
REPORT <- file.path(OUT, "65_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("65_figures_1_to_3 | ", format(Sys.time()), "\n", LOG_MSG, "\n", sep = "")
  
  ## ---------------- Figure 1 ----------------
  cat("\n== Figure 1: study region ==\n")
  tryCatch({
    cls <- SITES$coordinate_source; lv <- c("FIELD", "AUTHOR_REPORTED", "GIS_MAP")
    cols <- c(FIELD = OI[["green"]], AUTHOR_REPORTED = OI[["orange"]], GIS_MAP = OI[["blue"]]); pchs <- c(FIELD = 16, AUTHOR_REPORTED = 17, GIS_MAP = 15)
    cat("sites by class:", paste(sprintf("%s %d", lv, vapply(lv, function(k) sum(cls == k), 1L)), collapse = ", "), "| window area", sprintf("%.1f km2", area.owin(W_km)), "\n")
    save_fig("Fig01_study_region", 115, 150, function() {
      par(mar = c(3.4, 3.8, 0.6, 0.6), mgp = c(2.2, 0.6, 0), cex = 0.8, las = 1)
      plot(NA, xlim = XR, ylim = YR, asp = 1, xlab = "Easting (km, UTM zone 44N)", ylab = "Northing (km, UTM zone 44N)", bty = "l")
      outline(0.5, "grey30", "grey96")
      for (k in lv) points(SITES$x_m[cls == k] / 1000, SITES$y_m[cls == k] / 1000, pch = pchs[[k]], col = cols[[k]], cex = 0.7)
      legend("bottomleft", legend = sprintf("%s (n = %d)", c("Field", "Author reported", "GIS map"), vapply(lv, function(k) sum(cls == k), 1L)), pch = pchs[lv], col = cols[lv], bty = "n", cex = 0.8, title = "Coordinate source", title.adj = 0)
      x0 <- XR[2] - 70; y0 <- YR[1] + 8; segments(x0, y0, x0 + 50, y0, lwd = 2); text(x0 + 25, y0 + 5, "50 km", cex = 0.7)
      xa <- XR[2] - 12; arrows(xa, YR[2] - 28, xa, YR[2] - 8, length = 0.08, lwd = 1.5); text(xa, YR[2] - 4, "N", cex = 0.8, font = 2)
    })
    REG <- rbind(REG, data.frame(figure = "Fig01_study_region", content = "window and sites by coordinate source", sites_plotted = 197))
  }, error = function(e) cat("Figure 1 FAILED:", conditionMessage(e), "\n"))
  
  ## ---------------- Figure 2 ----------------
  cat("\n== Figure 2: database reconciliation ==\n")
  tryCatch({
    box_txt <- list(
      c("Master inventory", sprintf("N = %d records", tn("records in sites_master_decided.csv"))),
      c("Identity / duplicate audit", sprintf("%d candidate pairs reviewed: all confirmed distinct (no merge)", tn("candidate pairs reviewed and confirmed distinct (no merge)")),
        sprintf("%d records excluded as duplicates", tn("records excluded as duplicates"))),
      c("Provenance audit", sprintf("%d author-stated decisions; %d default inclusions (not individually confirmed)", tn("decisions stated by the author (AUTHOR_STATED)"), tn("decisions by default inclusion (DEFAULT_UNQUEUED, not individually confirmed)")),
        sprintf("%d author-corrected fields; %d redeposited sites retained (sensitivity set)", tn("records with author-corrected fields, included"), tn("records flagged redeposited (included; sensitivity set)"))),
      c("Coordinate audit", sprintf("Field %d | Author reported %d | GIS map %d (datum unknown)", tn("FIELD coordinates"), tn("AUTHOR_REPORTED coordinates"), tn("GIS_MAP coordinates (datum unknown)"))),
      c("Inclusion criteria and window screening", sprintf("%d records excluded; %d inside the analytical window", tn("records excluded"), tn("records inside the analytical window"))),
      c("Analytical sample", sprintf("N_A = %d sites (dataset A)", tn("sites in primary dataset A"))))
    cat(paste(vapply(box_txt, function(b) paste(b, collapse = " / "), ""), collapse = "\n"), "\n")
    save_fig("Fig02_database_reconciliation", 150, 150, function() {
      par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10.4))
      nb <- length(box_txt); h <- 1.3; gap <- 0.35; fills <- c("#E8E8E8", rep("#DCEBF7", nb - 2), "#FCE9C8")
      for (i in seq_len(nb)) {
        top <- 10.2 - (i - 1) * (h + gap); bot <- top - h
        rect(0.6, bot, 9.4, top, col = fills[i], border = "grey30", lwd = 0.8)
        text(5, top - 0.28, box_txt[[i]][1], font = 2, cex = 0.85)
        if (length(box_txt[[i]]) > 1) for (j in 2:length(box_txt[[i]])) text(5, top - 0.28 - 0.3 * (j - 1) - 0.12, box_txt[[i]][j], cex = 0.68)
        if (i < nb) arrows(5, bot, 5, bot - gap + 0.02, length = 0.07, lwd = 1.2)
      }
    })
    REG <- rbind(REG, data.frame(figure = "Fig02_database_reconciliation", content = "197 to N_A flow from T1a", sites_plotted = NA))
  }, error = function(e) cat("Figure 2 FAILED:", conditionMessage(e), "\n"))
  
  ## ---------------- Figure 3 ----------------
  cat("\n== Figure 3: covariates and observation offset ==\n")
  tryCatch({
    cont_cols <- hcl.colors(100, "viridis")
    cb <- function(zlim, cols) {
      n <- length(cols); xs <- seq(XR[1] + 0.1 * diff(XR), XR[1] + 0.9 * diff(XR), length.out = n + 1); y0 <- YR[1] - 0.075 * diff(YR); y1 <- YR[1] - 0.035 * diff(YR)
      par(xpd = NA); rect(xs[-(n + 1)], y0, xs[-1], y1, col = cols, border = NA)
      text(xs[1], y0, sprintf("%.1f", zlim[1]), adj = c(0, 1.3), cex = 0.6); text(xs[n + 1], y0, sprintf("%.1f", zlim[2]), adj = c(1, 1.3), cex = 0.6); par(xpd = FALSE)
    }
    pn <- function(im, title, letter, cat_cols = NULL, cat_labs = NULL, z = NULL) {
      plot(NA, xlim = XR, ylim = YR, asp = 1, axes = FALSE, xlab = "", ylab = ""); v <- if (is.null(z)) im$v else z
      if (is.null(cat_cols)) { zl <- quantile(v, c(0.005, 0.995), na.rm = TRUE); v <- pmin(pmax(v, zl[1]), zl[2]); image(x = im$xcol, y = im$yrow, z = t(v), col = cont_cols, zlim = zl, add = TRUE); outline(); cb(zl, cont_cols)
      cat(sprintf("  (%s) %-44s range %.2f to %.2f | valid pixels %d\n", letter, title, min(v, na.rm = TRUE), max(v, na.rm = TRUE), sum(!is.na(v))))
      } else { k <- length(cat_cols); image(x = im$xcol, y = im$yrow, z = t(v), col = cat_cols, breaks = seq(0.5, k + 0.5, by = 1), add = TRUE); outline()
      legend("bottomleft", legend = cat_labs, fill = cat_cols, bty = "n", cex = 0.58, border = NA, inset = c(0, -0.02))
      cat(sprintf("  (%s) %-44s classes: %s\n", letter, title, paste(sprintf("%s %d", cat_labs, as.integer(table(factor(v[!is.na(v)], levels = seq_len(k))))), collapse = ", "))) }
      mtext(sprintf("(%s) %s", letter, title), side = 3, adj = 0, line = 0.1, cex = 0.6, font = 2)
    }
    lith <- 1 + covs$lith_basalt$v + 2 * covs$lith_sandstone$v + 3 * covs$lith_carbonate$v + 4 * covs$lith_other$v
    geom <- 1 + covs$geom_dissected$v + 2 * covs$geom_alluvial$v
    save_fig("Fig03_covariates_observation", 180, 205, function() {
      layout(matrix(1:9, 3, byrow = TRUE)); par(mar = c(2.4, 0.3, 1.3, 0.3), cex = 0.8)
      pn(covs$z_elev_p1, "Elevation, linear term (z)", "a"); pn(covs$z_elev_p2, "Elevation, quadratic term (z)", "b"); pn(covs$z_lslope, "log(1 + slope), z", "c")
      pn(covs$z_tpi2km, "TPI, 2 km (z)", "d"); pn(covs$z_ldist_riv3, "log distance to order >= 3 rivers (z)", "e"); pn(covs$z_ldist_riv5, "log distance to order >= 5 rivers (z)", "f")
      pn(covs$z_elev_p1, "Lithology (5 classes)", "g", cat_cols = c(OI[["grey"]], OI[["orange"]], OI[["sky"]], OI[["green"]], OI[["pink"]]), cat_labs = c("felsic (ref.)", "basalt", "sandstone", "carbonate", "other"), z = lith)
      pn(covs$z_elev_p1, "Geomorphology (3 classes)", "h", cat_cols = c(OI[["grey"]], OI[["verm"]], OI[["blue"]]), cat_labs = c("pediplain (ref.)", "dissected", "alluvial"), z = geom)
      pn(OBS$logE_O2, "Fixed observation offset O2 (log E, assumed)", "i")
    })
    REG <- rbind(REG, data.frame(figure = "Fig03_covariates_observation", content = "8 covariate panels + O2 offset; no sites", sites_plotted = 0))
  }, error = function(e) cat("Figure 3 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Files ==\n")
  fl <- list.files(OUT, pattern = "^Fig0[1-3].*\\.(png|pdf)$", full.names = TRUE)
  FI <- data.frame(file = basename(fl), kb = round(file.size(fl) / 1024, 1)); print(FI, row.names = FALSE)
  write.csv(REG, file.path(OUT, "65_figure_register.csv"), row.names = FALSE)
  cat("expected 6 files (3 figures x png/pdf); found", nrow(FI), "\n\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT, "65_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")