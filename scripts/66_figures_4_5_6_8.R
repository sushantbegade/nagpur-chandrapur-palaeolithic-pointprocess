# =====================================================================
# 66_figures_4_5_6_8.R
# Stage   : 14_figures
# Folder  : 14_figures/scripts/   (NEW)
# Purpose : Figure 4 (observed L_inhom - r, g_inhom, K_inhom vs the M1 pointwise band), Figure 5 (M1
#           generative diagnostics: standardised deviations; null distribution of the per-simulation
#           max |z| with the observed value), Figure 6 (M3 diagnostics: L and g vs M3 band; bootstrap
#           distributions of sigma2 and alpha), Figure 8 (observed vs simulated structure across models;
#           raw max |z| by model). Also a text-fit check of Figure 2 and a wrapped re-render of it.
# Outputs : outputs/Fig04_empirical_statistics.*, Fig05_M1_diagnostics.*, Fig06_M3_diagnostics.*,
#           Fig08_model_comparison.*, Fig02_database_reconciliation.* (re-rendered), 66_report.txt,
#           66_figure_register.csv, 66_sessionInfo.txt
# Notes   : No analysis is run: figures from saved outputs (42/43/44/46 curves, 43 simulations, 46
#           bootstrap, T5, T1a). Pointwise bands are display only; inference uses the global ERL test.
#           K_inhom is derived as pi * L^2. Okabe-Ito palette, ASCII labels. Restart R before running.
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)
source(file.path(REPO_DIR, "R", "freeze_v1_params.R"))

D11 <- dir_of("11_generative_gof", "outputs"); D15 <- dir_of("15_tables", "outputs")
OUT <- dir_of("14_figures", "outputs"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# ---- open inputs BEFORE sink ----
need <- c(file.path(D11, c("42_m1_curves.csv", "42_m1_result.csv", "44_m2_curves.csv", "46_m3_curves.csv", "46_m3_boot.csv", "46_m3_result.csv", "43_m1_sims.rds")),
          file.path(D15, c("T5_generative_gof.csv", "T1a_audit.csv")))
stopifnot(all(file.exists(need)))
C42 <- read.csv(file.path(D11, "42_m1_curves.csv"), stringsAsFactors = FALSE); R42 <- read.csv(file.path(D11, "42_m1_result.csv"), stringsAsFactors = FALSE)
C44 <- read.csv(file.path(D11, "44_m2_curves.csv"), stringsAsFactors = FALSE); C46 <- read.csv(file.path(D11, "46_m3_curves.csv"), stringsAsFactors = FALSE)
B46 <- read.csv(file.path(D11, "46_m3_boot.csv"), stringsAsFactors = FALSE); R46 <- read.csv(file.path(D11, "46_m3_result.csv"), stringsAsFactors = FALSE)
SM <- readRDS(file.path(D11, "43_m1_sims.rds"))
T5 <- read.csv(file.path(D15, "T5_generative_gof.csv"), stringsAsFactors = FALSE); T1a <- read.csv(file.path(D15, "T1a_audit.csv"), stringsAsFactors = FALSE)

OI <- c(black = "#000000", orange = "#E69F00", sky = "#56B4E9", green = "#009E73", yellow = "#F0E442", blue = "#0072B2", verm = "#D55E00", pink = "#CC79A7", grey = "#BBBBBB")
save_fig <- function(name, w_mm, h_mm, fun) {
  for (dev in c("png", "pdf")) {
    f <- file.path(OUT, paste0(name, ".", dev))
    if (dev == "png") png(f, width = round(w_mm / 25.4 * 300), height = round(h_mm / 25.4 * 300), res = 300) else pdf(f, width = w_mm / 25.4, height = h_mm / 25.4)
    tryCatch(fun(), finally = dev.off())
  }
  invisible(NULL)
}
getc <- function(df, cv) df[df$curve == cv, ]
tn <- function(item) { v <- T1a$n[T1a$item == item]; if (length(v) != 1) stop("T1a item not found: ", item); v }
ribbon <- function(d, col) polygon(c(d$r_km, rev(d$r_km)), c(d$lo, rev(d$hi)), col = adjustcolor(col, 0.35), border = NA)
tag <- function(txt) mtext(txt, side = 3, adj = 0, line = 0.2, cex = 0.8, font = 2)
REG <- data.frame()

REPORT <- file.path(OUT, "66_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("66_figures_4_5_6_8 | ", format(Sys.time()), "\n", sep = "")
  stopifnot(setequal(unique(C42$curve), c("L_inhom_minus_r", "g_inhom")), all(c("scenario", "curve", "r_km", "obs", "sim_mean", "lo", "hi") %in% names(C44)))
  dL <- getc(C42, "L_inhom_minus_r"); dG <- getc(C42, "g_inhom")
  
  cat("\n== Figure 4: empirical statistics vs the M1 band ==\n")
  tryCatch({
    kf <- function(d) { K <- function(y) pi * (y + d$r_km)^2; d$obs <- K(d$obs); d$sim_mean <- K(d$sim_mean); d$lo <- K(d$lo); d$hi <- K(d$hi); d }
    dK <- kf(dL)
    chk <- do.call(rbind, lapply(c(2, 3, 5, 10), function(r) data.frame(r_km = r, L_obs = dL$obs[dL$r_km == r], L_hi = dL$hi[dL$r_km == r], g_obs = if (r %in% dG$r_km) dG$obs[dG$r_km == r] else NA,
                                                                        g_hi = if (r %in% dG$r_km) dG$hi[dG$r_km == r] else NA, K_obs = dK$obs[dK$r_km == r], K_sim_mean = dK$sim_mean[dK$r_km == r])))
    print(chk, digits = 4, row.names = FALSE)
    save_fig("Fig04_empirical_statistics", 180, 68, function() {
      par(mfrow = c(1, 3), mar = c(3.4, 3.9, 1.5, 0.6), mgp = c(2.3, 0.6, 0), cex = 0.8, las = 1)
      plot(NA, xlim = range(dL$r_km), ylim = range(c(dL$lo, dL$hi, dL$obs)), xlab = "Distance r (km)", ylab = "L_inhom(r) - r", bty = "l"); ribbon(dL, OI[["blue"]])
      lines(dL$r_km, dL$sim_mean, lty = 2, col = OI[["blue"]], lwd = 1.2); lines(dL$r_km, dL$obs, lwd = 2); abline(h = 0, col = "grey50", lty = 3); tag("(a)")
      legend("topleft", c("Observed", "M1 simulation mean", "M1 95% pointwise band"), lty = c(1, 2, NA), lwd = c(2, 1.2, NA), pch = c(NA, NA, 15), col = c("black", OI[["blue"]], adjustcolor(OI[["blue"]], 0.35)), bty = "n", cex = 0.62)
      plot(NA, xlim = range(dG$r_km), ylim = range(c(dG$lo, dG$hi, dG$obs)), xlab = "Distance r (km)", ylab = "g_inhom(r)", bty = "l"); ribbon(dG, OI[["blue"]])
      lines(dG$r_km, dG$sim_mean, lty = 2, col = OI[["blue"]], lwd = 1.2); lines(dG$r_km, dG$obs, lwd = 2); abline(h = 1, col = "grey50", lty = 3); tag("(b)")
      text(par("usr")[2], par("usr")[4], sprintf("global ERL p = %.3f (nsim %d)", R42$p_erl_nstep2, R42$nsim), adj = c(1.02, 1.5), cex = 0.62)
      plot(NA, xlim = range(dK$r_km), ylim = range(c(dK$lo, dK$hi, dK$obs)), xlab = "Distance r (km)", ylab = "K_inhom(r) (km^2)", bty = "l"); ribbon(dK, OI[["blue"]])
      lines(dK$r_km, dK$sim_mean, lty = 2, col = OI[["blue"]], lwd = 1.2); lines(dK$r_km, dK$obs, lwd = 2); tag("(c)")
    })
    REG <- rbind(REG, data.frame(figure = "Fig04_empirical_statistics", content = "L, g, K (derived) vs M1 band"))
  }, error = function(e) cat("Figure 4 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Figure 5: M1 generative diagnostics ==\n")
  tryCatch({
    zv <- function(o, M) (o - rowMeans(M)) / apply(M, 1, sd)
    zL <- zv(SM$ob$L, SM$sL); zG <- zv(SM$ob$g, SM$sG)
    zsim <- rbind(SM$sL - rowMeans(SM$sL), SM$sG - rowMeans(SM$sG)) / rbind(matrix(apply(SM$sL, 1, sd), nrow(SM$sL), ncol(SM$sL)), matrix(apply(SM$sG, 1, sd), nrow(SM$sG), ncol(SM$sG)))
    mz <- apply(abs(zsim), 2, max); obsmax <- max(abs(c(zL, zG)))
    cat(sprintf("observed max |z| %.2f | simulations: median %.2f, 95th %.2f, 99th %.2f, max %.2f | share of simulations with max |z| >= observed: %.3f (not the ERL p)\n",
                obsmax, median(mz), quantile(mz, .95), quantile(mz, .99), max(mz), mean(mz >= obsmax)))
    save_fig("Fig05_M1_diagnostics", 180, 70, function() {
      par(mfrow = c(1, 2), mar = c(3.4, 3.9, 1.5, 0.6), mgp = c(2.3, 0.6, 0), cex = 0.8, las = 1)
      plot(NA, xlim = range(SM$rL), ylim = range(c(zL, zG, -2, 2)), xlab = "Distance r (km)", ylab = "Standardised deviation (observed - mean) / sd", bty = "l")
      lines(SM$rL, zL, col = OI[["blue"]], lwd = 2); lines(SM$rG, zG, col = OI[["verm"]], lwd = 2); abline(h = c(-2, 2), lty = 3, col = "grey50"); abline(h = 0, col = "grey70"); tag("(a)")
      legend("topleft", c("L_inhom - r", "g_inhom", "+/- 2 (reference only)"), lty = c(1, 1, 3), lwd = c(2, 2, 1), col = c(OI[["blue"]], OI[["verm"]], "grey50"), bty = "n", cex = 0.62)
      hist(mz, breaks = 40, col = OI[["sky"]], border = "white", main = "", xlab = "Maximum |z| over r (per simulation)", ylab = "Simulations"); abline(v = obsmax, col = OI[["verm"]], lwd = 2); tag("(b)")
      text(par("usr")[2], par("usr")[4], sprintf("observed %.2f\nsimulation 95th percentile %.2f", obsmax, quantile(mz, .95)), adj = c(1.02, 1.3), cex = 0.62)
    })
    REG <- rbind(REG, data.frame(figure = "Fig05_M1_diagnostics", content = "standardised deviations; null max |z| distribution"))
  }, error = function(e) cat("Figure 5 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Figure 6: M3 diagnostics ==\n")
  tryCatch({
    eL <- getc(C46, "L_inhom_minus_r"); eG <- getc(C46, "g_inhom")
    cat(sprintf("M3 curves: L %d points, g %d points | observed L at 3 km %.3f (M3 band %.3f-%.3f) | g at 2 km %.3f (band %.3f-%.3f)\n", nrow(eL), nrow(eG), eL$obs[eL$r_km == 3], eL$lo[eL$r_km == 3], eL$hi[eL$r_km == 3],
                eG$obs[eG$r_km == 2], eG$lo[eG$r_km == 2], eG$hi[eG$r_km == 2]))
    nclip <- sum(B46$alpha >= 30, na.rm = TRUE)
    cat(sprintf("bootstrap refits %d | sigma2*: median %.3f | alpha*: median %.2f, share > %g km %.3f, clipped at 30 km: %d\n", nrow(B46), median(B46$sigma2, na.rm = TRUE), median(B46$alpha, na.rm = TRUE), R_MAX_DECISION_KM, mean(B46$alpha > R_MAX_DECISION_KM, na.rm = TRUE), nclip))
    save_fig("Fig06_M3_diagnostics", 180, 140, function() {
      par(mfrow = c(2, 2), mar = c(3.4, 3.9, 1.5, 0.6), mgp = c(2.3, 0.6, 0), cex = 0.8, las = 1)
      plot(NA, xlim = range(eL$r_km), ylim = range(c(eL$lo, eL$hi, eL$obs)), xlab = "Distance r (km)", ylab = "L_inhom(r) - r", bty = "l"); ribbon(eL, OI[["green"]])
      lines(eL$r_km, eL$sim_mean, lty = 2, col = OI[["green"]], lwd = 1.2); lines(eL$r_km, eL$obs, lwd = 2); abline(h = 0, col = "grey50", lty = 3); tag("(a)")
      legend("topleft", c("Observed", "M3 simulation mean", "M3 95% pointwise band"), lty = c(1, 2, NA), lwd = c(2, 1.2, NA), pch = c(NA, NA, 15), col = c("black", OI[["green"]], adjustcolor(OI[["green"]], 0.35)), bty = "n", cex = 0.62)
      plot(NA, xlim = range(eG$r_km), ylim = range(c(eG$lo, eG$hi, eG$obs)), xlab = "Distance r (km)", ylab = "g_inhom(r)", bty = "l"); ribbon(eG, OI[["green"]])
      lines(eG$r_km, eG$sim_mean, lty = 2, col = OI[["green"]], lwd = 1.2); lines(eG$r_km, eG$obs, lwd = 2); abline(h = 1, col = "grey50", lty = 3); tag("(b)")
      text(par("usr")[2], par("usr")[4], sprintf("global ERL p = %.3f (nsim %d, conditional)", R46$p_erl_nstep2, R46$nsim), adj = c(1.02, 1.5), cex = 0.62)
      hist(B46$sigma2, breaks = 25, col = OI[["sky"]], border = "white", main = "", xlab = "Bootstrap sigma2 (latent variance estimate)", ylab = "Refits"); abline(v = R46$sigma2_hat, col = OI[["verm"]], lwd = 2); tag("(c)")
      xa <- pmin(B46$alpha, 30); hist(xa, breaks = seq(0, 30, by = 1), col = OI[["sky"]], border = "white", main = "", xlab = "Bootstrap alpha (km; clipped at 30)", ylab = "Refits")
      abline(v = R46$alpha_hat, col = OI[["verm"]], lwd = 2); abline(v = R_MAX_DECISION_KM, lty = 2, col = "grey30"); tag("(d)")
      legend("topright", c("point estimate", "10 km decision range"), lty = c(1, 2), lwd = c(2, 1), col = c(OI[["verm"]], "grey30"), bty = "n", cex = 0.62)
    })
    REG <- rbind(REG, data.frame(figure = "Fig06_M3_diagnostics", content = "L, g vs M3 band; sigma2 and alpha bootstrap"))
  }, error = function(e) cat("Figure 6 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Figure 8: observed vs simulated structure across models ==\n")
  tryCatch({
    o2L <- C44[C44$scenario == "O2" & C44$curve == "L_inhom_minus_r", ]; o2G <- C44[C44$scenario == "O2" & C44$curve == "g_inhom", ]
    z <- T5[!is.na(T5$raw_max_z), ]; z$label <- ifelse(grepl("^M1", z$model), "M1", ifelse(grepl("^M3", z$model), "M3", sub("M2 \\((O.)\\)", "M2-\\1", z$model)))
    cat("raw max |z| by model:", paste(sprintf("%s %.2f", z$label, z$raw_max_z), collapse = ", "), "\n")
    save_fig("Fig08_model_comparison", 180, 70, function() {
      par(mfrow = c(1, 3), mar = c(3.4, 3.9, 1.5, 0.6), mgp = c(2.3, 0.6, 0), cex = 0.8, las = 1)
      plot(NA, xlim = range(dL$r_km), ylim = range(c(dL$lo, dL$hi, o2L$lo, o2L$hi, eL$lo, eL$hi, dL$obs)), xlab = "Distance r (km)", ylab = "L_inhom(r) - r", bty = "l")
      ribbon(dL, OI[["blue"]]); ribbon(o2L, OI[["orange"]]); ribbon(eL, OI[["green"]]); lines(dL$r_km, dL$obs, lwd = 2); lines(o2L$r_km, o2L$obs, lwd = 1.5, lty = 2, col = OI[["orange"]]); abline(h = 0, col = "grey50", lty = 3); tag("(a)")
      legend("topleft", c("Observed (E3 intensity)", "Observed (E3 + O2 offset)", "M1 band", "M2 (O2) band", "M3 band"), lty = c(1, 2, NA, NA, NA), lwd = c(2, 1.5, NA, NA, NA), pch = c(NA, NA, 15, 15, 15),
             col = c("black", OI[["orange"]], adjustcolor(OI[["blue"]], 0.5), adjustcolor(OI[["orange"]], 0.5), adjustcolor(OI[["green"]], 0.5)), bty = "n", cex = 0.55)
      plot(NA, xlim = range(dG$r_km), ylim = range(c(dG$lo, dG$hi, o2G$lo, o2G$hi, eG$lo, eG$hi, dG$obs)), xlab = "Distance r (km)", ylab = "g_inhom(r)", bty = "l")
      ribbon(dG, OI[["blue"]]); ribbon(o2G, OI[["orange"]]); ribbon(eG, OI[["green"]]); lines(dG$r_km, dG$obs, lwd = 2); lines(o2G$r_km, o2G$obs, lwd = 1.5, lty = 2, col = OI[["orange"]]); abline(h = 1, col = "grey50", lty = 3); tag("(b)")
      bp <- barplot(z$raw_max_z, names.arg = z$label, col = c(OI[["blue"]], OI[["orange"]], OI[["orange"]], OI[["green"]])[seq_len(nrow(z))], ylab = "Raw max |z| of the observed curves", las = 1, cex.names = 0.7, border = NA)
      text(bp, z$raw_max_z, sprintf("p %.3f\nnsim %d", z$p, z$nsim), pos = 3, cex = 0.55, xpd = NA); tag("(c)")
    })
    REG <- rbind(REG, data.frame(figure = "Fig08_model_comparison", content = "M1, M2(O2), M3 bands; raw max |z| by model"))
  }, error = function(e) cat("Figure 8 FAILED:", conditionMessage(e), "\n"))
  
  cat("\n== Figure 2: text-fit check and wrapped re-render ==\n")
  tryCatch({
    boxes <- list(
      c("Master inventory", sprintf("N = %d records", tn("records in sites_master_decided.csv"))),
      c("Identity / duplicate audit", sprintf("%d candidate pairs reviewed: all confirmed distinct (no merge)", tn("candidate pairs reviewed and confirmed distinct (no merge)")), sprintf("%d records excluded as duplicates", tn("records excluded as duplicates"))),
      c("Provenance audit", sprintf("%d author-stated decisions; %d default inclusions (not individually confirmed)", tn("decisions stated by the author (AUTHOR_STATED)"), tn("decisions by default inclusion (DEFAULT_UNQUEUED, not individually confirmed)")),
        sprintf("%d author-corrected fields; %d redeposited sites retained (sensitivity set)", tn("records with author-corrected fields, included"), tn("records flagged redeposited (included; sensitivity set)"))),
      c("Coordinate audit", sprintf("Field %d | Author reported %d | GIS map %d (datum unknown)", tn("FIELD coordinates"), tn("AUTHOR_REPORTED coordinates"), tn("GIS_MAP coordinates (datum unknown)"))),
      c("Inclusion criteria and window screening", sprintf("%d records excluded; %d inside the analytical window", tn("records excluded"), tn("records inside the analytical window"))),
      c("Analytical sample", sprintf("N_A = %d sites (dataset A)", tn("sites in primary dataset A"))))
    WMAX <- 8.6
    pdf(NULL, width = 150 / 25.4, height = 150 / 25.4); par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10.4))
    allb <- unlist(lapply(boxes, function(b) if (length(b) > 1) b[-1] else NULL)); w0 <- strwidth(allb, cex = 0.68, units = "user"); dev.off()
    cat(sprintf("original layout (cex 0.68, usable width %.1f user units): %d body lines, widest %.2f, lines exceeding: %d\n", WMAX, length(allb), max(w0), sum(w0 > WMAX)))
    if (any(w0 > WMAX)) cat("overflowing:", paste(substr(allb[w0 > WMAX], 1, 60), collapse = " | "), "\n")
    fig2 <- function() {
      par(mar = c(0.3, 0.3, 0.3, 0.3)); plot.new(); plot.window(xlim = c(0, 10), ylim = c(0, 10.4))
      wrap <- function(txt, cex) { words <- strsplit(txt, " ")[[1]]; out <- character(0); cur <- ""
      for (w in words) { trial <- if (nchar(cur)) paste(cur, w) else w; if (strwidth(trial, cex = cex, units = "user") <= WMAX || !nchar(cur)) cur <- trial else { out <- c(out, cur); cur <- w } }; c(out, cur) }
      cex_b <- 0.68
      repeat { LL <- lapply(boxes, function(b) list(title = b[1], body = if (length(b) > 1) unlist(lapply(b[-1], wrap, cex = cex_b)) else character(0)))
      hh <- vapply(LL, function(l) 0.55 + 0.3 * length(l$body), 1); tot <- sum(hh) + 0.3 * (length(LL) - 1); if (tot <= 10.1 || cex_b <= 0.5) break; cex_b <- cex_b - 0.04 }
      cat(sprintf("re-render: body text cex %.2f, total height %.2f of 10.1, widest wrapped line %.2f (limit %.1f)\n", cex_b, tot, max(vapply(LL, function(l) if (length(l$body)) max(strwidth(l$body, cex = cex_b, units = "user")) else 0, 1)), WMAX))
      ytop <- 10.2; fills <- c("#E8E8E8", rep("#DCEBF7", length(LL) - 2), "#FCE9C8")
      for (i in seq_along(LL)) { bot <- ytop - hh[i]; rect(0.6, bot, 9.4, ytop, col = fills[i], border = "grey30", lwd = 0.8); text(5, ytop - 0.28, LL[[i]]$title, font = 2, cex = 0.85)
      for (j in seq_along(LL[[i]]$body)) text(5, ytop - 0.28 - 0.3 * j, LL[[i]]$body[j], cex = cex_b)
      if (i < length(LL)) arrows(5, bot, 5, bot - 0.28, length = 0.07, lwd = 1.2); ytop <- bot - 0.3 }
    }
    save_fig("Fig02_database_reconciliation", 150, 150, fig2)
    REG <- rbind(REG, data.frame(figure = "Fig02_database_reconciliation", content = "re-rendered with wrapped text"))
  }, error = function(e) { cat("Figure 2 re-render FAILED:", conditionMessage(e), "\n"); while (dev.cur() > 1) dev.off() })
  
  cat("\n== Files ==\n")
  fl <- list.files(OUT, pattern = "^Fig0[2468].*\\.(png|pdf)$", full.names = TRUE); FI <- data.frame(file = basename(fl), kb = round(file.size(fl) / 1024, 1)); print(FI, row.names = FALSE)
  write.csv(REG, file.path(OUT, "66_figure_register.csv"), row.names = FALSE)
  cat("expected 10 files (Fig02, 04, 05, 06, 08 x png/pdf); found", nrow(FI), "\n\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT, "66_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")