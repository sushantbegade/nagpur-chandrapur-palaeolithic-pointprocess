# =====================================================================
# 68_supplementary_figures_and_captions.R   (v2: par() style applied through do.call; caption numbers listed)
# Stage   : 14_figures
# Folder  : 14_figures/scripts/   (OVERWRITES v1)
# Purpose : supplementary figures S1-S5 and a caption draft (Figures 1-11, S1-S5) from saved outputs:
#           S1 ERL p versus nsim and self-calibration (scripts 51, 52, 50), S2 source structure of close
#           pairs (scripts 60, 62, 63), S3 block validation and coefficient stability (scripts 56-57),
#           S4 observed curves under four robustness settings (script 49), S5 anisotropy sectors and
#           GIS-class composition of close pairs (scripts 58, 60).
# Outputs : outputs/FigS1_p_vs_nsim.*, FigS2_source_structure.*, FigS3_block_validation.*,
#           FigS4_robustness_curves.*, FigS5_anisotropy_composition.*, figure_captions_draft.txt,
#           68_report.txt, 68_figure_register.csv, 68_sessionInfo.txt
# Notes   : No analysis is run. v1 passed the style list PAR to par() unnamed, so R ignored it;
#           v2 uses setpar() (do.call). Captions are generated from table values and carry the caveats;
#           they are a DRAFT for the author to edit. All numeric tokens in the captions are printed so
#           they can be checked against the tables. Okabe-Ito palette, ASCII only.
#           Restart R before running. Main session only.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)
source(file.path(REPO_DIR, "R", "freeze_v1_params.R"))

D11 <- dir_of("11_generative_gof", "outputs"); D12 <- dir_of("12_validation", "outputs"); D13 <- dir_of("13_sensitivity", "outputs"); D15 <- dir_of("15_tables", "outputs")
OUT <- dir_of("14_figures", "outputs"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- open inputs BEFORE sink ----
need <- c(file.path(D11, c("51_nsim_consistency.csv", "52_selfcal.csv", "50_jitter_replicates.csv", "49_curves_all.csv")), file.path(D12, c("57_blockcv_exploratory.csv", "57_coef_stability.csv")),
          file.path(D13, c("62_author_share_bands.csv", "62_block_perm.csv", "62_loao_result.csv", "63_matched_control_summary.csv", "63_matched_control_replicates.csv", "60_pair_composition.csv", "58_anisotropy_sectors.csv")),
          file.path(D15, c("T5_generative_gof.csv", "T4b_M3_parameters.csv", "T6_stage1_recovery.csv", "T7b_jitter.csv", "T7a_robustness_M1_M3.csv")))
stopifnot(all(file.exists(need))); rdc <- function(p) read.csv(p, stringsAsFactors = FALSE)
S51 <- rdc(need[1]); S52 <- rdc(need[2]); J50 <- rdc(need[3]); C49 <- rdc(need[4]); V57 <- rdc(need[5]); CS57 <- rdc(need[6]); AU <- rdc(need[7]); BP <- rdc(need[8]); LO <- rdc(need[9])
CM <- rdc(need[10]); CR <- rdc(need[11]); PC <- rdc(need[12]); SEC <- rdc(need[13]); T5 <- rdc(need[14]); T4b <- rdc(need[15]); T6 <- rdc(need[16]); T7b <- rdc(need[17]); T7a <- rdc(need[18])

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
tag <- function(txt) mtext(txt, side = 3, adj = 0, line = 0.2, cex = 0.75, font = 2)
ribbon <- function(d, col) polygon(c(d$r_km, rev(d$r_km)), c(d$lo, rev(d$hi)), col = adjustcolor(col, 0.35), border = NA)
PAR <- list(mar = c(3.4, 3.9, 1.5, 0.6), mgp = c(2.3, 0.6, 0), cex = 0.78, las = 1)
setpar <- function(mfrow, mar = PAR$mar) { a <- PAR; a$mar <- mar; a$mfrow <- mfrow; do.call(par, a) }
REG <- data.frame()

REPORT <- file.path(OUT, "68_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("68_supplementary_figures_and_captions v2 | ", format(Sys.time()), "\n", sep = "")
  
  cat("\n== S1: ERL p versus nsim ==\n")
  tryCatch({
    pse <- S52[S52$level == "pseudo", ]; re <- S52[S52$level == "real", ]; bs <- J50[J50$type == "base" & !is.na(J50$p), ]
    ideal <- function(N) (ceiling(0.05 * (N + 1)) - 1) / (N + 1)
    sc <- do.call(rbind, lapply(sort(unique(pse$N)), function(N) { d <- pse[pse$N == N & !is.na(pse$p), ]; k <- sum(d$p < 0.05); ci <- wil(k, nrow(d)); data.frame(N = N, rate = k / nrow(d), lo = ci[1], hi = ci[2], ideal = ideal(N)) }))
    mp <- tapply(re$p, re$N, mean); p999 <- T5$p[T5$model == "M1 (O1)" & T5$setting == "simplified window"]
    print(sc, digits = 3, row.names = FALSE); cat("real observed, mean p by N:", paste(sprintf("N %s %.5f", names(mp), mp), collapse = " | "), "| nsim 999 primary", p999, "\n")
    cat(sprintf("script-50 baseline (nsim 199): %d replicates, median p %.3f, share < 0.05 %.2f\n", nrow(bs), median(bs$p), mean(bs$p < 0.05)))
    save_fig("FigS1_p_vs_nsim", 180, 130, function() {
      setpar(c(2, 2))
      boxplot(comb ~ N, data = S51, xlab = "Simulations per subset (N)", ylab = "Combined ERL p (random subsets of the 999)", col = OI[["sky"]], border = "grey30", outline = FALSE); abline(h = c(0.031, 0.05), lty = c(1, 2), col = c(OI[["verm"]], "grey40")); tag("(a)")
      text(par("usr")[2], 0.031, "full set 0.031", adj = c(1.05, -0.5), cex = 0.6, col = OI[["verm"]])
      xs <- seq_len(nrow(sc)); plot(NA, xlim = c(0.5, nrow(sc) + 0.5), ylim = c(0, 0.12), xaxt = "n", xlab = "Simulations per test (N)", ylab = "P(p < 0.05) for exchangeable curves", bty = "l"); axis(1, at = xs, labels = sc$N)
      segments(xs, sc$lo, xs, sc$hi, lwd = 1.5, col = OI[["blue"]]); points(xs, sc$rate, pch = 19, col = OI[["blue"]]); points(xs, sc$ideal, pch = 4, col = OI[["verm"]], cex = 1.2); tag("(b)")
      legend("topleft", c("observed rate (Wilson 95%)", "ideal exact-MC level"), pch = c(19, 4), col = c(OI[["blue"]], OI[["verm"]]), bty = "n", cex = 0.6)
      hist(bs$p, breaks = seq(0, max(0.2, max(bs$p) + 0.01), by = 0.005), col = OI[["sky"]], border = "white", main = "", xlab = "p of fresh baseline replicates (nsim 199)", ylab = "Replicates"); abline(v = c(0.031, 0.05), lty = c(1, 2), col = c(OI[["verm"]], "grey40")); tag("(c)")
      xn <- c(as.numeric(names(mp)), 999); yn <- c(as.numeric(mp), p999); plot(xn, yn, log = "x", pch = 19, col = OI[["blue"]], ylim = c(0, 0.05), xlab = "Simulations (N, log scale)", ylab = "Mean p of the observed pattern", bty = "l", type = "b"); abline(h = 0.05, lty = 2, col = "grey40"); tag("(d)")
    })
    REG <- rbind(REG, data.frame(figure = "FigS1_p_vs_nsim", content = "p versus nsim; self-calibration; baseline p"))
  }, error = function(e) cat("S1 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== S2: source structure ==\n")
  tryCatch({
    cat("share same author by band:", paste(sprintf("%s %.3f", AU$band_km, AU$share_same_author), collapse = " | "), "\n")
    cat("LOAO p:", paste(sprintf("%s %.3f", LO$author_removed, LO$p), collapse = " | "), "| matched control retention G", CM$retention[CM$group == "G"], "R", CM$retention[CM$group == "R"], "\n")
    save_fig("FigS2_source_structure", 180, 130, function() {
      setpar(c(2, 2)); x <- seq_len(nrow(AU))
      plot(NA, xlim = c(0.7, nrow(AU) + 0.3), ylim = c(0, 0.75), xaxt = "n", xlab = "Distance band (km)", ylab = "Share of site pairs", bty = "l"); axis(1, at = x, labels = AU$band_km, cex.axis = 0.6)
      lines(x, AU$share_same_author, type = "b", pch = 19, col = OI[["blue"]]); lines(x, AU$share_same_reference, type = "b", pch = 17, col = OI[["orange"]]); tag("(a)")
      legend("topright", c("same lead author", "same source reference"), pch = c(19, 17), col = c(OI[["blue"]], OI[["orange"]]), lty = 1, bty = "n", cex = 0.6)
      bl <- sort(unique(BP$block_km)); bands <- unique(BP$band_km); plot(NA, xlim = c(0.5, length(bands) + 0.5), ylim = range(c(BP$perm_lo, BP$perm_hi, BP$observed_same_author)), xaxt = "n", xlab = "Distance band (km)", ylab = "Same-author pairs", bty = "l"); axis(1, at = seq_along(bands), labels = bands)
      for (k in seq_along(bl)) { d <- BP[BP$block_km == bl[k], ]; xp <- seq_along(bands) + (k - 1.5) * 0.18; segments(xp, d$perm_lo, xp, d$perm_hi, lwd = 3, col = c(OI[["sky"]], OI[["green"]])[k]); points(xp, d$observed_same_author, pch = 19) }
      legend("topleft", c("observed", sprintf("shuffled within %g km blocks (95%%)", bl)), pch = c(19, NA, NA), lty = c(NA, 1, 1), lwd = c(NA, 3, 3), col = c("black", OI[["sky"]], OI[["green"]]), bty = "n", cex = 0.58); tag("(b)")
      yy <- rev(seq_len(nrow(LO))); plot(NA, xlim = c(1e-3, 1), ylim = c(0.5, nrow(LO) + 0.5), log = "x", yaxt = "n", xlab = "M1 ERL p (log scale), nsim 999", ylab = "", bty = "l"); axis(2, at = yy, labels = paste("minus", LO$author_removed), cex.axis = 0.65, tick = FALSE)
      segments(pmax(LO$p - LO$mc_se, 1e-3), yy, pmin(LO$p + LO$mc_se, 1), yy, col = OI[["blue"]], lwd = 1.5); points(LO$p, yy, pch = 19, col = OI[["blue"]]); abline(v = c(0.031, 0.05), lty = c(1, 2), col = c(OI[["verm"]], "grey40")); tag("(c) leave one author out")
      eg <- ecdf(CR$p[CR$grp == "G"]); er <- ecdf(CR$p[CR$grp == "R"]); pp <- seq(0, 0.3, by = 0.002); plot(pp, eg(pp), type = "s", col = OI[["verm"]], lwd = 2, ylim = c(0, 1), xlab = "M1 ERL p (nsim 199)", ylab = "Cumulative share of replicates", bty = "l")
      lines(pp, er(pp), type = "s", col = OI[["blue"]], lwd = 2); abline(v = 0.05, lty = 2, col = "grey40"); tag("(d) matched control")
      legend("bottomright", c(sprintf("GIS_MAP sites removed (retention %.2f)", CM$retention[CM$group == "G"]), sprintf("60 random sites removed (retention %.2f)", CM$retention[CM$group == "R"])), lty = 1, lwd = 2, col = c(OI[["verm"]], OI[["blue"]]), bty = "n", cex = 0.58)
    })
    REG <- rbind(REG, data.frame(figure = "FigS2_source_structure", content = "same-author shares; block permutation; leave-one-author-out; matched control"))
  }, error = function(e) cat("S2 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== S3: block validation ==\n")
  tryCatch({
    pool <- function(d, v) sum(d[[v]] * d$n_test) / sum(d$n_test)
    cmb <- expand.grid(model = c("E1", "E3"), block_km = sort(unique(V57$block_km)), stringsAsFactors = FALSE)
    gl <- lapply(seq_len(nrow(cmb)), function(i) { d <- V57[V57$model == cmb$model[i] & V57$block_km == cmb$block_km[i], ]; vapply(split(d, d$rep), pool, 1, v = "gain") })
    pl <- lapply(seq_len(nrow(cmb)), function(i) { d <- V57[V57$model == cmb$model[i] & V57$block_km == cmb$block_km[i], ]; vapply(split(d, d$rep), pool, 1, v = "pct") })
    nm <- sprintf("%s\n%g km", cmb$model, cmb$block_km); cc <- ifelse(cmb$model == "E3", OI[["blue"]], OI[["orange"]])
    cat("median gain:", paste(sprintf("%s %.3f", gsub("\n", " ", nm), vapply(gl, median, 1)), collapse = " | "), "\nmedian percentile:", paste(sprintf("%s %.3f", gsub("\n", " ", nm), vapply(pl, median, 1)), collapse = " | "), "\n")
    save_fig("FigS3_block_validation", 180, 80, function() {
      setpar(c(1, 3), mar = c(4.4, 3.9, 1.5, 0.6))
      boxplot(gl, names = nm, col = cc, border = "grey30", ylab = "Pooled held-out log-likelihood gain per site", cex.axis = 0.55, outline = FALSE); abline(h = 0, lty = 2, col = "grey40"); tag("(a)")
      boxplot(pl, names = nm, col = cc, border = "grey30", ylab = "Mean percentile of held-out sites (0.5 = chance)", cex.axis = 0.55, outline = FALSE); abline(h = 0.5, lty = 2, col = "grey40"); tag("(b)")
      o <- order(CS57$same_sign_share); par(mar = c(4.4, 7.2, 1.5, 0.6)); barplot(CS57$same_sign_share[o], horiz = TRUE, names.arg = CS57$term[o], col = ifelse(CS57$same_sign_share[o] >= 0.8, OI[["green"]], OI[["verm"]]), border = NA, xlab = "Share of training fits with the full-data sign", cex.names = 0.55, xlim = c(0, 1)); abline(v = 0.8, lty = 2, col = "grey40"); tag("(c)")
    })
    REG <- rbind(REG, data.frame(figure = "FigS3_block_validation", content = "gain and percentile by block size; E3 coefficient sign stability"))
  }, error = function(e) cat("S3 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== S4: observed curves under four robustness settings ==\n")
  tryCatch({
    sets <- c("R04_MP", "R16_zS3", "R06_eroded5", "R02_datasetC"); lab <- c("MP subset (n 114)", "free z_S3 (adds core term)", "eroded window 5 km (n 169)", "Dataset C (n 97)")
    stopifnot(all(sets %in% C49$setting)); pv <- setNames(T7a$p[match(sets, T7a$setting)], sets)
    cat("p:", paste(sprintf("%s %.4f", sets, pv), collapse = " | "), "\n")
    save_fig("FigS4_robustness_curves", 180, 100, function() {
      setpar(c(2, 4), mar = c(3.2, 3.6, 1.6, 0.4))
      for (cv in c("L_inhom_minus_r", "g_inhom")) for (k in seq_along(sets)) { d <- C49[C49$setting == sets[k] & C49$curve == cv, ]
      plot(NA, xlim = range(d$r_km), ylim = range(c(d$lo, d$hi, d$obs)), xlab = "r (km)", ylab = if (cv == "g_inhom") "g_inhom(r)" else "L_inhom(r) - r", bty = "l"); ribbon(d, OI[["blue"]]); lines(d$r_km, d$sim_mean, lty = 2, col = OI[["blue"]]); lines(d$r_km, d$obs, lwd = 2)
      abline(h = if (cv == "g_inhom") 1 else 0, lty = 3, col = "grey50"); if (cv == "L_inhom_minus_r") mtext(sprintf("%s  p = %.3f", lab[k], pv[k]), side = 3, adj = 0, line = 0.2, cex = 0.55, font = 2) }
    })
    REG <- rbind(REG, data.frame(figure = "FigS4_robustness_curves", content = "L and g under MP, z_S3, eroded window, Dataset C"))
  }, error = function(e) cat("S4 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== S5: anisotropy sectors and GIS-class composition ==\n")
  tryCatch({
    cat("sector shares:", paste(sprintf("%s %.3f", SEC$sector_deg, SEC$obs_share), collapse = " | "), "\n")
    d2 <- PC[PC$band != "0-5 km", ]; d2$lab <- sprintf("%s\n%s", sub(" km", "", d2$band), ifelse(grepl("both", d2$statistic), "both GIS", "one GIS")); cat("composition rows:", nrow(d2), "\n")
    save_fig("FigS5_anisotropy_composition", 180, 75, function() {
      setpar(c(1, 2), mar = c(4.2, 3.9, 1.5, 0.6)); x <- seq_len(nrow(SEC))
      plot(NA, xlim = c(0.5, nrow(SEC) + 0.5), ylim = range(c(SEC$M3_lo, SEC$M3_hi, SEC$M1_lo, SEC$M1_hi, SEC$obs_share)), xaxt = "n", xlab = "Pair orientation (degrees from east, modulo 180)", ylab = "Share of 2-10 km pairs", bty = "l"); axis(1, at = x, labels = SEC$sector_deg, cex.axis = 0.55)
      segments(x - 0.12, SEC$M1_lo, x - 0.12, SEC$M1_hi, lwd = 3, col = OI[["sky"]]); segments(x + 0.12, SEC$M3_lo, x + 0.12, SEC$M3_hi, lwd = 3, col = OI[["green"]]); points(x, SEC$obs_share, pch = 19); tag("(a)")
      legend("topright", c("observed", "M1 95% range", "M3 95% range"), pch = c(19, NA, NA), lty = c(NA, 1, 1), lwd = c(NA, 3, 3), col = c("black", OI[["sky"]], OI[["green"]]), bty = "n", cex = 0.58)
      xx <- seq_len(nrow(d2)); plot(NA, xlim = c(0.5, nrow(d2) + 0.5), ylim = range(c(d2$perm_lo, d2$perm_hi, d2$observed)), xaxt = "n", xlab = "Distance band and statistic", ylab = "Number of pairs", bty = "l"); axis(1, at = xx, labels = d2$lab, cex.axis = 0.55, padj = 0.4)
      segments(xx, d2$perm_lo, xx, d2$perm_hi, lwd = 3, col = OI[["grey"]]); points(xx, d2$observed, pch = 19); tag("(b)")
      legend("topright", c("observed", "label permutations (95%)"), pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 3), col = c("black", OI[["grey"]]), bty = "n", cex = 0.58)
    })
    REG <- rbind(REG, data.frame(figure = "FigS5_anisotropy_composition", content = "sector shares vs nulls; GIS-class pair composition"))
  }, error = function(e) cat("S5 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Caption draft ==\n")
  tryCatch({
    pm <- function(mod, set) T5$p[T5$model == mod & T5$setting == set]
    p1 <- pm("M1 (O1)", "simplified window"); p1x <- pm("M1 (O1)", "exact window"); p2 <- pm("M2 (O2)", "precision rerun"); p2b <- T5$p[T5$model == "M2 (O3)"]; p3 <- T5$p[grepl("^M3", T5$model)]
    s2 <- T4b$estimate[1]; al <- T4b$estimate[2]; ga <- T6[T6$truth == "Geyer", ]
    cap <- c("FIGURE CAPTION DRAFT (generated from saved outputs; edit before use)", "",
             "Figure 1. Study area (Nagpur and Chandrapur districts, Vidarbha), analytical window (21,262 km2) and the 197 Palaeolithic sites by coordinate source (field 64, author reported 73, GIS map 60; GIS-map datum unknown). UTM zone 44N.", "",
             "Figure 2. Reconciliation of the master inventory (197 records) to the analytical sample N_A = 197. No record was excluded or merged; 30 inclusion decisions are stated by the author and 167 are default inclusions that were not individually confirmed.", "",
             "Figure 3. Environmental covariates (present-day proxies; terrain, hydrology and geomorphology also carry preservation and visibility roles) and the fixed observation offset O2 (an assumed-strength proxy from roads and land cover, not an estimated survey-effort surface). Colour scales are clipped at the 0.5% and 99.5% quantiles. Elevation panels show orthogonal polynomial terms, not metres. No sites are drawn.", "",
             sprintf("Figure 4. Observed inhomogeneous L(r) - r, g(r) and K(r) (K derived as pi L^2; dashed line = pi (mean L)^2, approximate) against 999 simulations of the environmental Poisson model M1. The band is the pointwise 95%% range, uncorrected and for display only; inference uses the global envelope test (combined ERL p = %.3f).", p1), "",
             "Figure 5. M1 diagnostics. (a) Standardised deviations of the observed curves from the simulation mean; the +/-2 lines are a Gaussian reference only, because the null is heavy tailed. (b) Distribution of the per-simulation maximum |z| with the observed maximum marked.", "",
             sprintf("Figure 6. M3 (environment plus latent log-Gaussian field; conditional simulation) diagnostics. (a, b) Observed curves against the M3 pointwise band (combined ERL p = %.3f). (c, d) Parametric bootstrap of the estimator outputs sigma2 (estimate %.2f) and alpha (estimate %.1f km). These are estimator outputs, not physical parameters or distances; alpha is poorly determined and is not translated into territory or mobility.", p3, s2, al), "",
             "Figure 7. Classification of simulated patterns by the three-test hierarchy (100 replicates per recovery cell, nsim 99; M4 was not fitted to the data). Right-hand text: M1 rejection rate, replicates that reached M3 and M4.", "",
             sprintf("Figure 8. Observed against simulated structure under M1, M2 (fixed offset O2) and M3. The two observed curves differ because the inhomogeneous intensity is refitted under each model. (c) Raw maximum |z| of the observed curves with the global p and nsim (M1 %.3f, M2-O2 %.4f at nsim 10000, M2-O3 %.3f, M3 %.3f).", p1, p2, p2b, p3), "",
             sprintf("Figure 9. Stage 1 recovery: (a) power for Gibbs interaction (%.2f-%.2f; frozen threshold 0.80), (b) recovery of the LGCP mechanism, (c) false-positive rate in non-Gibbs cells (0 of 1300 replicates; frozen limit 0.10), (d) first-test (M1) rejection rate. Wilson 95%% intervals. The tested interaction strength stops at gamma = 2.5.", min(ga$P_Geyer), max(ga$P_Geyer)), "",
             "Figure 10. Robustness of the M1 result: (a) global p by assumption category (bars: Monte Carlo standard error; nsim 999 unless marked), (b) observed / simulated pair-count ratios, (c) retention of the result under coordinate jitter (author / GIS-map sigma in metres) against the frozen 0.80 criterion. Descriptive profile, not a count of tests passed.", "",
             "Figure 11. Model hierarchy with results and interpretation levels. Level 3 explanations were not tested; statistical dependence is a property of the recorded pattern and is not evidence of a behavioural mechanism.", "",
             "Figure S1. Behaviour of the combined ERL p-value across the number of simulations. (a) Random subsets of the 999 simulations. (b) Self-calibration on exchangeable simulated curves against the ideal exact Monte Carlo level. (c) Fresh baseline replicates at nsim 199. (d) Mean p of the observed pattern against N.", "",
             "Figure S2. Source structure of close pairs (exploratory, post hoc). (a) Share of site pairs sharing a lead author or source reference by distance. (b) Same-author pairs against lead-author labels shuffled within 30 and 60 km blocks. (c) M1 p with each large author removed in turn. (d) Matched control: GIS-map sites removed against 60 random sites removed.", "",
             "Figure S3. Spatial block validation (first-order structure only). (a) Pooled held-out log-likelihood gain per site over a homogeneous intensity. (b) Mean percentile of held-out sites by predicted intensity (0.5 = chance). (c) Sign stability of the E3 coefficients across training fits at 60 km. A good rank score does not support M1 or M3 as complete explanations.", "",
             "Figure S4. Observed curves under four robustness settings (pointwise bands as in Figure 4). Panel titles give the global p at nsim 999 (Dataset C: nsim 10000 value in Table 7a).", "",
             "Figure S5. (a) Pair orientations (2-10 km) against the M1 and M3 null ranges: no directional structure was detected, and no anisotropic model was fitted. (b) Pairs involving GIS-map sites against label permutations (exploratory). Orientation is a description of pair vectors, not a corridor or behavioural axis.")
    writeLines(cap, file.path(OUT, "figure_captions_draft.txt")); cat("lines written:", length(cap), "\n")
    tok <- unique(unlist(regmatches(cap, gregexpr("[0-9]+(?:[.,][0-9]+)*", cap, perl = TRUE))))
    cat("numeric tokens in the captions (check each against the tables):\n", paste(tok, collapse = " "), "\n")
  }, error = function(e) cat("Caption draft FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Supplementary files ==\n")
  fl <- list.files(OUT, pattern = "^FigS[1-5].*\\.(png|pdf)$", full.names = TRUE); FI <- data.frame(file = basename(fl), kb = round(file.size(fl) / 1024, 1)); print(FI, row.names = FALSE)
  write.csv(REG, file.path(OUT, "68_figure_register.csv"), row.names = FALSE); cat("expected 10 files; found", nrow(FI), "| caption draft present:", file.exists(file.path(OUT, "figure_captions_draft.txt")), "\n\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT, "68_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")