# =====================================================================
# 67_figures_7_9_10_11.R
# Stage   : 14_figures
# Folder  : 14_figures/scripts/   (NEW)
# Purpose : Figure 7 (Stage 1 classification per recovery cell; replaces the Geyer-diagnostics figure),
#           Figure 9 (Stage 1 recovery / power metrics with Wilson intervals vs frozen thresholds),
#           Figure 10 (robustness profile by assumption category: p on a log axis, pair ratios, jitter
#           retention), Figure 11 (conceptual diagram: hierarchy with results; interpretation levels).
#           All values are read from 15_tables CSVs (and 62_author_share_bands.csv for the author share).
# Outputs : outputs/Fig07_stage1_classification.*, Fig09_recovery_power.*, Fig10_robustness.*,
#           Fig11_conceptual_diagram.*, 67_report.txt, 67_figure_register.csv, 67_sessionInfo.txt
# Notes   : No analysis is run. Okabe-Ito palette, ASCII labels. Figure 11 states no behavioural
#           mechanism. A file listing of ALL figure files (Fig01-Fig11) is printed at the end.
#           Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)
source(file.path(REPO_DIR, "R", "freeze_v1_params.R"))

D15 <- dir_of("15_tables", "outputs"); D13 <- dir_of("13_sensitivity", "outputs")
OUT <- dir_of("14_figures", "outputs"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- open inputs BEFORE sink ----
need <- c(file.path(D15, c("T6_stage1_recovery.csv", "T7a_robustness_M1_M3.csv", "T7b_jitter.csv", "T5_generative_gof.csv", "T4b_M3_parameters.csv")), file.path(D13, "62_author_share_bands.csv"))
stopifnot(all(file.exists(need)))
T6 <- read.csv(need[1], stringsAsFactors = FALSE); T7a <- read.csv(need[2], stringsAsFactors = FALSE); T7b <- read.csv(need[3], stringsAsFactors = FALSE)
T5 <- read.csv(need[4], stringsAsFactors = FALSE); T4b <- read.csv(need[5], stringsAsFactors = FALSE); AU <- read.csv(need[6], stringsAsFactors = FALSE)

OI <- c(black = "#000000", orange = "#E69F00", sky = "#56B4E9", green = "#009E73", yellow = "#F0E442", blue = "#0072B2", verm = "#D55E00", pink = "#CC79A7", grey = "#BBBBBB")
save_fig <- function(name, w_mm, h_mm, fun) {
  for (dev in c("png", "pdf")) {
    f <- file.path(OUT, paste0(name, ".", dev))
    if (dev == "png") png(f, width = round(w_mm / 25.4 * 300), height = round(h_mm / 25.4 * 300), res = 300) else pdf(f, width = w_mm / 25.4, height = h_mm / 25.4)
    tryCatch(fun(), finally = dev.off())
  }
  invisible(NULL)
}
wil <- function(k, n) as.numeric(suppressWarnings(prop.test(round(k), n, correct = FALSE)$conf.int))
tag <- function(txt) mtext(txt, side = 3, adj = 0, line = 0.2, cex = 0.8, font = 2)
REG <- data.frame()

REPORT <- file.path(OUT, "67_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("67_figures_7_9_10_11 | ", format(Sys.time()), "\n", sep = "")
  ord <- c("P1", "P2", "N1", "N2", "N3", "N4", "L1", "L2", "L3", "L4", "N5", "C2", "C4", "G1", "G2", "G3", "G4", "C1", "C3")
  stopifnot(setequal(ord, T6$cell)); d6 <- T6[match(ord, T6$cell), ]
  
  cat("\n== Figure 7: Stage 1 classification ==\n")
  tryCatch({
    M <- rbind(Poisson = d6$class_Poisson, LGCP = d6$class_LGCP, Geyer = d6$class_Geyer, inadequate = d6$class_inadequate) / d6$replicates
    cat("class totals over 1900 replicates:", paste(sprintf("%s %d", c("Poisson", "LGCP", "Geyer", "inadequate"), c(sum(d6$class_Poisson), sum(d6$class_LGCP), sum(d6$class_Geyer), sum(d6$class_inadequate))), collapse = ", "), "\n")
    cat("column sums of the plotted proportions (all must be 1):", paste(round(colSums(M), 3), collapse = " "), "\n")
    lab <- sprintf("%s  %s%s", d6$cell, d6$truth, ifelse(d6$n_target < 197, sprintf(" (n %d)", d6$n_target), ""))
    save_fig("Fig07_stage1_classification", 170, 140, function() {
      par(mar = c(3.6, 8.2, 1.6, 9.6), mgp = c(2.3, 0.6, 0), cex = 0.8, las = 1)
      bp <- barplot(M[, rev(seq_len(ncol(M)))], horiz = TRUE, col = c(OI[["grey"]], OI[["sky"]], OI[["verm"]], OI[["yellow"]]), border = "white", names.arg = rev(lab), xlab = "Share of replicates (100 per cell, nsim 99)", cex.names = 0.7, xlim = c(0, 1))
      par(xpd = NA); text(1.03, bp, sprintf("M1 %.2f | M3 %d | M4 %d", rev(d6$M1_rejection), rev(d6$reached_M3), rev(d6$reached_M4)), adj = 0, cex = 0.55)
      text(1.03, max(bp) + 0.95, "M1 rejection | reached M3 | reached M4", adj = 0, cex = 0.55, font = 2); par(xpd = FALSE)
      legend("bottom", inset = c(0, -0.2), horiz = TRUE, legend = c("Poisson", "LGCP", "Geyer", "inadequate"), fill = c(OI[["grey"]], OI[["sky"]], OI[["verm"]], OI[["yellow"]]), bty = "n", cex = 0.62, xpd = NA); tag("Classification of simulated patterns by the three-test hierarchy")
    })
    REG <- rbind(REG, data.frame(figure = "Fig07_stage1_classification", content = "class shares per recovery cell"))
  }, error = function(e) cat("Figure 7 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Figure 9: recovery and power metrics ==\n")
  tryCatch({
    panel_ci <- function(d, est, lo, hi, thr, thr_lab, ylab_txt, ttl, ylim = c(0, 1)) {
      n <- nrow(d); x <- seq_len(n); plot(NA, xlim = c(0.4, n + 0.6), ylim = ylim, xaxt = "n", xlab = "", ylab = ylab_txt, bty = "l")
      axis(1, at = x, labels = d$cell, cex.axis = 0.7); if (!is.na(thr)) { abline(h = thr, lty = 2, col = OI[["verm"]]); text(n + 0.55, thr, thr_lab, adj = c(1, -0.4), cex = 0.55, col = OI[["verm"]]) }
      segments(x, d[[lo]], x, d[[hi]], lwd = 1.5, col = OI[["blue"]]); points(x, d[[est]], pch = 19, col = OI[["blue"]], cex = 0.9); tag(ttl)
    }
    ga <- d6[d6$truth == "Geyer", ]; lg <- d6[d6$truth == "LGCP", ]; nb <- d6[d6$truth != "Geyer", ]
    cat(sprintf("Geyer power: %s | LGCP recovery: %s | Geyer-call rate in non-Gibbs cells (13 cells): %d of %d replicates\n", paste(sprintf("%s %.2f", ga$cell, ga$P_Geyer), collapse = ", "),
                paste(sprintf("%s %.2f", lg$cell, lg$P_correct_class), collapse = ", "), sum(nb$class_Geyer), sum(nb$replicates)))
    m1 <- do.call(rbind, lapply(seq_len(nrow(d6)), function(i) { ci <- wil(d6$M1_rejection[i] * d6$replicates[i], d6$replicates[i]); data.frame(cell = d6$cell[i], est = d6$M1_rejection[i], lo = ci[1], hi = ci[2]) }))
    save_fig("Fig09_recovery_power", 180, 130, function() {
      par(mfrow = c(2, 2), mar = c(3.2, 3.9, 1.6, 0.6), mgp = c(2.3, 0.6, 0), cex = 0.8, las = 1)
      panel_ci(ga, "P_Geyer", "Geyer_lo", "Geyer_hi", DECISION$power_min, "frozen threshold 0.80", "Power: P(Geyer | Gibbs truth)", "(a) Gibbs-truth cells")
      panel_ci(lg, "P_correct_class", "correct_lo", "correct_hi", DECISION$mechanism_recovery_min, "frozen threshold 0.70", "Recovery: P(LGCP | LGCP truth)", "(b) LGCP-truth cells")
      panel_ci(nb, "P_Geyer", "Geyer_lo", "Geyer_hi", DECISION$fpr_max, "frozen limit 0.10", "False-positive rate: P(Geyer | no Gibbs)", "(c) non-Gibbs cells", ylim = c(0, 0.2))
      panel_ci(m1, "est", "lo", "hi", NA, "", "M1 rejection rate (alpha 0.05)", "(d) all cells, first test of the hierarchy")
    })
    REG <- rbind(REG, data.frame(figure = "Fig09_recovery_power", content = "power, recovery, false positives, M1 rejection with Wilson intervals"))
  }, error = function(e) cat("Figure 9 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Figure 10: robustness profile ==\n")
  tryCatch({
    catmap <- c(R00_baseline = "control (new seed)", R01_redep = "site definition", R02_datasetC = "site definition", R03_LP = "chronology", R04_MP = "chronology", R05_UP = "chronology",
                R06_eroded5 = "window", R18_exact_window = "window", R07_E1 = "environmental specification", R08_E2 = "environmental specification", R09_E4 = "environmental specification",
                R10_grid250 = "resolution", R11_grid1000 = "resolution", R12_quad05 = "resolution", R13_edge_translate = "edge correction", R14_edge_border_L = "edge correction",
                R15_zS2 = "observation proxy / offset", R16_zS3 = "observation proxy / offset", M2_O2 = "observation proxy / offset", M2_O3 = "observation proxy / offset")
    r <- T7a[T7a$test == "M1", c("setting", "n", "nsim", "p", "mc_se", "ratio_2", "ratio_3", "ratio_5")]
    o2 <- T5[T5$model == "M2 (O2)" & T5$setting == "precision rerun", ]; o3 <- T5[T5$model == "M2 (O3)", ]
    r <- rbind(r, data.frame(setting = c("M2_O2", "M2_O3"), n = 197, nsim = c(o2$nsim, o3$nsim), p = c(o2$p, o3$p), mc_se = c(o2$mc_se, o3$mc_se), ratio_2 = NA, ratio_3 = NA, ratio_5 = NA))
    stopifnot(all(r$setting %in% names(catmap))); r$cat <- catmap[r$setting]
    corder <- c("control (new seed)", "site definition", "chronology", "window", "environmental specification", "resolution", "edge correction", "observation proxy / offset")
    r <- r[order(match(r$cat, corder)), ]; ccol <- setNames(c(OI[["black"]], OI[["orange"]], OI[["green"]], OI[["blue"]], OI[["pink"]], OI[["sky"]], OI[["verm"]], "grey40"), corder)
    r$lab <- sprintf("%s (n %d)", gsub("_", " ", sub("^R[0-9]+_", "", r$setting)), r$n); ny <- nrow(r); y <- rev(seq_len(ny))
    cat(sprintf("settings plotted: %d | p < 0.05: %d | p >= 0.05: %s\n", ny, sum(r$p < 0.05), paste(r$setting[r$p >= 0.05], collapse = ", ")))
    save_fig("Fig10_robustness", 180, 150, function() {
      layout(matrix(c(1, 2, 3, 3), 2, byrow = TRUE), widths = c(1.55, 1), heights = c(2.2, 1)); par(cex = 0.8, las = 1, mgp = c(2.3, 0.6, 0))
      par(mar = c(3.4, 11.5, 1.6, 0.6)); plot(NA, xlim = c(5e-4, 1), ylim = c(0.5, ny + 0.5), log = "x", yaxt = "n", xlab = "ERL p-value (log scale)", ylab = "", bty = "l")
      axis(2, at = y, labels = r$lab, cex.axis = 0.58, tick = FALSE); abline(v = 0.05, lty = 2, col = "grey30")
      segments(pmax(r$p - r$mc_se, 5e-4), y, pmin(r$p + r$mc_se, 1), y, col = ccol[r$cat], lwd = 1.5); points(r$p, y, pch = ifelse(r$nsim >= 10000, 18, 19), col = ccol[r$cat], cex = 1.1)
      text(0.05, ny + 0.45, "0.05", cex = 0.55, adj = c(-0.2, 0)); tag("(a) M1 test under each assumption")
      legend("bottomleft", legend = corder, col = ccol, pch = 19, bty = "n", cex = 0.5, ncol = 1); legend("bottomright", c("nsim 999", "nsim 10000"), pch = c(19, 18), bty = "n", cex = 0.55, col = "grey30")
      par(mar = c(3.4, 0.6, 1.6, 0.6)); rr <- r[, c("ratio_2", "ratio_3", "ratio_5")]; plot(NA, xlim = c(0.8, max(rr, na.rm = TRUE) * 1.05), ylim = c(0.5, ny + 0.5), yaxt = "n", xlab = "Observed / simulated pairs", ylab = "", bty = "l")
      abline(v = 1, lty = 3, col = "grey50"); points(rr$ratio_2, y, pch = 1, cex = 0.9); points(rr$ratio_3, y, pch = 2, cex = 0.9); points(rr$ratio_5, y, pch = 0, cex = 0.9); tag("(b) pair-count ratios")
      legend("topright", c("2 km", "3 km", "5 km"), pch = c(1, 2, 0), bty = "n", cex = 0.6)
      par(mar = c(5.2, 3.9, 1.6, 0.6)); nj <- nrow(T7b); xj <- seq_len(nj); plot(NA, xlim = c(0.4, nj + 0.6), ylim = c(0.6, 1), xaxt = "n", xlab = "", ylab = "Share of replicates with p < 0.05", bty = "l")
      jl <- sub(" \\(FIELD 10 m\\)", "", sub("primary jitter \\(FIELD 10, AUTHOR 25, GIS 250 m\\)", "primary jitter", sub(" \\(nominal coordinates\\)", "", T7b$scenario)))
      axis(1, at = xj, labels = sprintf("%s\n(%d / %d m)", jl, T7b$author_m, T7b$gis_m), cex.axis = 0.52, padj = 0.5)
      abline(h = 0.8, lty = 2, col = OI[["verm"]]); segments(xj, T7b$ci_lo, xj, T7b$ci_hi, lwd = 1.5, col = OI[["blue"]]); points(xj, T7b$retention, pch = 19, col = OI[["blue"]]); tag("(c) coordinate jitter (author / GIS sigma)")
      text(nj + 0.55, 0.8, "frozen criterion 0.80", adj = c(1, 1.5), cex = 0.55, col = OI[["verm"]])
    })
    REG <- rbind(REG, data.frame(figure = "Fig10_robustness", content = "M1 p by assumption category; pair ratios; jitter retention"))
  }, error = function(e) cat("Figure 10 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Figure 11: conceptual diagram ==\n")
  tryCatch({
    p1 <- T5$p[T5$model == "M1 (O1)" & T5$setting == "simplified window"]; p1x <- T5$p[T5$model == "M1 (O1)" & T5$setting == "exact window"]
    p2a <- T5$p[T5$model == "M2 (O2)" & T5$setting == "precision rerun"]; n2a <- T5$nsim[T5$model == "M2 (O2)" & T5$setting == "precision rerun"]; p2b <- T5$p[T5$model == "M2 (O3)"]
    p3 <- T5$p[grepl("^M3", T5$model)]; s2 <- T4b$estimate[1]; al <- T4b$estimate[2]
    base <- T7a[T7a$setting == "R00_baseline" & T7a$test == "M1", ]; rmin <- min(base$ratio_2, base$ratio_3, base$ratio_5); rmax <- max(base$ratio_2, base$ratio_3, base$ratio_5)
    ga <- T6[T6$truth == "Geyer", ]; sa <- sum(AU$same_author[AU$band_km %in% c("0-2", "2-5")]); np <- sum(AU$pairs[AU$band_km %in% c("0-2", "2-5")])
    lgc <- c(sum(T6$class_LGCP[T6$truth == "LGCP"]), sum(T6$class_LGCP[T6$truth == "Geyer"]), sum(T6$class_LGCP[T6$truth == "Poisson"]))
    cat(sprintf("values used: M1 p %.3f (exact %.3f) | M2 O2 %.4f (nsim %d), O3 %.3f | M3 p %.3f, sigma2 %.2f, alpha %.1f | ratios %.1f-%.1f | Geyer power %.2f-%.2f | same-author share of pairs within 5 km %.2f (%d of %d) | LGCP-classed replicates %d / %d / %d\n",
                p1, p1x, p2a, n2a, p2b, p3, s2, al, rmin, rmax, min(ga$P_Geyer), max(ga$P_Geyer), sa / np, sa, np, lgc[1], lgc[2], lgc[3]))
    box <- function(x0, x1, ytop, title, body, fill, cexb = 0.6) {
      wrap <- function(txt) { words <- strsplit(txt, " ")[[1]]; out <- character(0); cur <- ""
      for (w in words) { tr <- if (nchar(cur)) paste(cur, w) else w; if (strwidth(tr, cex = cexb, units = "user") <= (x1 - x0 - 0.4) || !nchar(cur)) cur <- tr else { out <- c(out, cur); cur <- w } }; c(out, cur) }
      L <- unlist(lapply(body, wrap)); h <- 0.55 + 0.27 * length(L); rect(x0, ytop - h, x1, ytop, col = fill, border = "grey30", lwd = 0.8); text((x0 + x1) / 2, ytop - 0.27, title, font = 2, cex = 0.7)
      for (j in seq_along(L)) text(x0 + 0.2, ytop - 0.27 - 0.27 * j, L[j], adj = 0, cex = cexb); invisible(h)
    }
    save_fig("Fig11_conceptual_diagram", 180, 150, function() {
      par(mar = c(0.2, 0.2, 0.2, 0.2)); plot.new(); plot.window(xlim = c(0, 20), ylim = c(0, 17)); text(5.2, 16.7, "Model hierarchy and results", font = 2, cex = 0.85); text(15.2, 16.7, "Interpretation levels", font = 2, cex = 0.85)
      y <- 16.2; fills <- c("#DCEBF7", "#DCEBF7", "#DCEBF7", "#E8E8E8")
      L <- list(list("M1: measured environment (E3 trend)", sprintf("Not sufficient. Observed site pairs are %.1f-%.1f times the simulated mean at 2-5 km. Global ERL p = %.3f (simplified window), %.3f (exact window).", rmin, rmax, p1, p1x)),
                list("M2: + observation offsets (assumed strengths)", sprintf("Not sufficient. O2 p = %.3f (nsim %d); O3 p = %.3f.", p2a, n2a, p2b)),
                list("M3: + latent spatial field (LGCP)", sprintf("Not rejected (p = %.3f). Estimator outputs: sigma2 %.2f, alpha %.1f km (not a distance; wide bootstrap interval).", p3, s2, al)),
                list("M4: + Gibbs interaction (Geyer)", sprintf("Not fitted. Stage 1 recovery power for Gibbs interaction at N_A = 197 was %.2f-%.2f (Gate 6 fails).", min(ga$P_Geyer), max(ga$P_Geyer))))
      for (i in seq_along(L)) { h <- box(0.3, 10.1, y, L[[i]][[1]], L[[i]][2], fills[i]); if (i < length(L)) arrows(5.2, y - h, 5.2, y - h - 0.55, length = 0.08, lwd = 1.2); if (i < length(L)) text(5.5, y - h - 0.28, "if not sufficient", adj = 0, cex = 0.55); y <- y - h - 0.6 }
      yy <- 16.2
      h <- box(10.7, 19.8, yy, "Level 1: observation", c(sprintf("The recorded pattern has an excess of close pairs at 2-5 km (no excess at 1 km). %.0f%% of site pairs within 5 km share a lead author.", 100 * sa / np)), "#FCE9C8"); yy <- yy - h - 0.5
      h <- box(10.7, 19.8, yy, "Level 2: statistical inference", c(sprintf("Environment and assumed observation offsets do not reproduce the pair excess; a latent-field model does. An LGCP verdict arose from Gibbs-truth patterns almost as often (%d) as from LGCP-truth patterns (%d), so latent clustering and interaction cannot be separated here.", lgc[2], lgc[1])), "#FCE9C8"); yy <- yy - h - 0.5
      box(10.7, 19.8, yy, "Level 3: archaeological interpretation (not tested)", c("Competing explanations: survey and discovery structure, repeated records of one locality, shared landscapes or recurrent use, chronological superposition, preservation, site spacing.", "Statistical dependence is a property of the recorded pattern, not proof of a behavioural mechanism."), "#EFEFEF")
    })
    REG <- rbind(REG, data.frame(figure = "Fig11_conceptual_diagram", content = "hierarchy with results and interpretation levels"))
  }, error = function(e) cat("Figure 11 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== All figure files ==\n")
  fl <- list.files(OUT, pattern = "^Fig(0[1-9]|1[01]).*\\.(png|pdf)$", full.names = TRUE); FI <- data.frame(file = basename(fl), kb = round(file.size(fl) / 1024, 1)); print(FI, row.names = FALSE)
  write.csv(REG, file.path(OUT, "67_figure_register.csv"), row.names = FALSE)
  cat("figures present (1-11):", paste(sort(unique(sub("^Fig(\\d+)_.*", "\\1", FI$file))), collapse = " "), "| files:", nrow(FI), "(expected 20 for 10 figures; Fig 7 here, Fig 1-6, 8-11)\n\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT, "67_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")