# =============================================================================
# 11_covariate_stack_audit.R
# Stage: 03_environment | Gate 4 (environment) closure support
# Purpose: fill touching-edge NA, transform (log slope, log TRI), winsorise by rule,
#   orthogonal quadratic elevation, standardise over window-centre cells, save
#   z-scored layers, collinearity audit (correlations, eta2, Cramer's V, VIF) for
#   E1-E4 at 250/500/1000 m, simulate the pre-specified drop rule.
# Does NOT read site data. No site-by-covariate tables. No models.
# Requires outputs of scripts 09 and 10.
# Outputs:
#   data_derived/env/<res>m/z_*.tif, transform_params.rds
#   03_environment/outputs/: covariate_stack_report.txt, transform_params.csv,
#     corr_matrix_<res>m.csv, vif_by_set.csv, categorical_association.csv,
#     drop_rule_result.csv, processing_log_11.csv, sessionInfo_11.txt
#   docs/covariate_definitions.csv (transform rows updated)
# Not executed by assistant. Run top to bottom in RStudio.
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
RES_SET   <- c(250, 500, 1000)
PRIMARY   <- 500
Z_BEYOND  <- 5          # |z| threshold
WINS_FRAC <- 0.001      # winsorise if more than this fraction of centre cells beyond |z|>Z_BEYOND
WINS_Q    <- c(0.005, 0.995)
VIF_MAX   <- 5
R_FLAG    <- 0.7
DROP_ORDER <- c("z_ldist_riv5", "z_elev_p2", "z_tpi2km", "z_ldist_riv3")

while (sink.number() > 0) sink()
source(REPO_CONFIG)
suppressPackageStartupMessages(library(terra))

out <- dir_of("03_environment"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
env_dir <- file.path(REPO_DIR, "data_derived", "env")
for (r in RES_SET) {
  d <- file.path(env_dir, paste0(r, "m"))
  stopifnot(all(file.exists(file.path(d, paste0(c("elev", "slope", "tri", "tpi2km", "ldist_riv3", "ldist_riv5",
                                                  "lith5", "geomorph3", "inwin", "incenter"), ".tif")))))
}
defs_f <- file.path(REPO_DIR, "docs", "covariate_definitions.csv")
defs <- if (file.exists(defs_f)) read.csv(defs_f, stringsAsFactors = FALSE) else NULL

logrows <- list()
lg <- function(step, detail) logrows[[length(logrows) + 1]] <<- data.frame(step = step, detail = detail, time = format(Sys.time()), stringsAsFactors = FALSE)
options(width = 200)

fill_mean <- function(r, inw, max_iter = 4) {
  cnt <- function(x) global(is.na(x) & !is.na(inw), "sum", na.rm = TRUE)[1, 1]
  n0 <- cnt(r)
  for (k in seq_len(max_iter)) { if (cnt(r) == 0) break
    r <- focal(r, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE) }
  list(r = r, n0 = n0, n1 = cnt(r))
}
vif_of <- function(M) {
  R <- tryCatch(cor(M), error = function(e) NULL)
  if (is.null(R)) return(setNames(rep(NA_real_, ncol(M)), colnames(M)))
  v <- tryCatch(diag(solve(R)), error = function(e) rep(Inf, ncol(M)))
  setNames(v, colnames(M))
}
apply_rule <- function(M) {
  dropped <- character(0)
  repeat {
    v <- vif_of(M)
    if (all(is.finite(v)) && max(v) <= VIF_MAX) break
    cand <- intersect(DROP_ORDER, colnames(M)); if (!length(cand)) break
    M <- M[, setdiff(colnames(M), cand[1]), drop = FALSE]; dropped <- c(dropped, cand[1])
  }
  list(M = M, dropped = dropped, vif = vif_of(M))
}
eta2 <- function(x, g) { m <- mean(x); st <- sum((x - m)^2); mg <- tapply(x, g, mean); ng <- tapply(x, g, length); sum(ng * (mg - m)^2) / st }

sets <- list(
  E1 = c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5"),
  E2 = c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5",
         "lith_basalt", "lith_sandstone", "lith_carbonate", "lith_other"),
  E3 = c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5",
         "lith_basalt", "lith_sandstone", "lith_carbonate", "lith_other", "geom_dissected", "geom_alluvial"),
  E4 = c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km",
         "lith_basalt", "lith_sandstone", "lith_carbonate", "lith_other"))

sink(file.path(out, "covariate_stack_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== 0. SELF-TEST ====\n")
  f <- rast(nrows = 3, ncols = 3, vals = 1:9); f[5] <- NA
  fv <- values(focal(f, w = 3, fun = "mean", na.policy = "only", na.rm = TRUE), mat = FALSE)
  cat("focal mean fill, centre =", fv[5], "(expect 5)\n")
  if (is.na(fv[5]) || abs(fv[5] - 5) > 1e-9) stop("focal(mean, na.policy='only') self-test FAILED. Paste this output.")
  cat("Passed. terra", as.character(packageVersion("terra")), "\n")
  cat("Rules: winsorise at q", WINS_Q, "if >", 100 * WINS_FRAC, "% of centre cells have |z| >", Z_BEYOND,
      "| VIF limit", VIF_MAX, "| |r| flag", R_FLAG, "| drop order:", paste(sub("^z_", "", DROP_ORDER), collapse = " > "), "\n")
  
  par_rows <- list(); vif_rows <- list(); drop_rows <- list(); cat_rows <- list()
  
  for (res in RES_SET) {
    cat(sprintf("\n=================== %d m GRID ===================\n", res))
    d  <- file.path(env_dir, paste0(res, "m")); rd <- function(n) rast(file.path(d, paste0(n, ".tif")))
    inw <- rd("inwin"); ic <- !is.na(values(rd("incenter"), mat = FALSE))
    cat("Window-centre cells:", sum(ic), "\n")
    
    # ---- fill ----
    raw <- list(elev = rd("elev"), slope = rd("slope"), tri = rd("tri"), tpi2km = rd("tpi2km"),
                ldist_riv3 = rd("ldist_riv3"), ldist_riv5 = rd("ldist_riv5"))
    nfill <- list()
    for (nm in names(raw)) { fl <- fill_mean(raw[[nm]], inw); raw[[nm]] <- fl$r; nfill[[nm]] <- c(fl$n0, fl$n1)
    cat(sprintf("fill %-10s touching-NA cells %5d -> remaining %d\n", nm, fl$n0, fl$n1)) }
    
    # ---- transforms ----
    cont <- list(elev = raw$elev, lslope = log(1 + raw$slope), tpi2km = raw$tpi2km, ltri = log(1 + raw$tri),
                 ldist_riv3 = raw$ldist_riv3, ldist_riv5 = raw$ldist_riv5)
    
    # ---- winsorise rule ----
    cat("\nWinsorise check (centre cells):\n")
    for (nm in names(cont)) {
      v <- values(cont[[nm]], mat = FALSE); vc <- v[ic]
      z <- (vc - mean(vc, na.rm = TRUE)) / sd(vc, na.rm = TRUE)
      fr <- mean(abs(z) > Z_BEYOND, na.rm = TRUE); do_w <- fr > WINS_FRAC
      qb <- quantile(vc, WINS_Q, na.rm = TRUE)
      if (do_w) cont[[nm]] <- clamp(cont[[nm]], qb[1], qb[2], values = TRUE)
      cat(sprintf("%-10s share |z|>%d: %.4f%% | max |z| %.2f | winsorised: %s\n", nm, Z_BEYOND, 100 * fr,
                  max(abs(z), na.rm = TRUE), ifelse(do_w, "YES", "no")))
      par_rows[[length(par_rows) + 1]] <- data.frame(res_m = res, layer = nm, frac_beyond_z = fr,
                                                     max_abs_z = max(abs(z), na.rm = TRUE), winsorised = do_w,
                                                     lo = if (do_w) qb[1] else NA, hi = if (do_w) qb[2] else NA,
                                                     n_fill_before = nfill[[ifelse(nm == "lslope", "slope", ifelse(nm == "ltri", "tri", nm))]][1], row.names = NULL)
    }
    
    # ---- orthogonal quadratic elevation ----
    ve <- values(cont$elev, mat = FALSE); okv <- !is.na(ve)
    pp <- poly(ve[ic & okv], degree = 2); cf <- attr(pp, "coefs")
    P <- matrix(NA_real_, nrow = length(ve), ncol = 2)
    P[okv, ] <- poly(ve[okv], degree = 2, coefs = cf)
    stack_in <- list(elev_p1 = setValues(cont$elev, P[, 1]), elev_p2 = setValues(cont$elev, P[, 2]),
                     lslope = cont$lslope, tpi2km = cont$tpi2km, ltri = cont$ltri,
                     ldist_riv3 = cont$ldist_riv3, ldist_riv5 = cont$ldist_riv5)
    
    # ---- standardise + save ----
    zr <- list(); Zc <- list()
    for (nm in names(stack_in)) {
      v <- values(stack_in[[nm]], mat = FALSE); m <- mean(v[ic], na.rm = TRUE); s <- sd(v[ic], na.rm = TRUE)
      z <- (stack_in[[nm]] - m) / s; names(z) <- paste0("z_", nm)
      zr[[paste0("z_", nm)]] <- mask(z, inw)
      writeRaster(zr[[paste0("z_", nm)]], file.path(d, paste0("z_", nm, ".tif")), overwrite = TRUE, datatype = "FLT4S")
      Zc[[paste0("z_", nm)]] <- values(zr[[paste0("z_", nm)]], mat = FALSE)[ic]
      par_rows[[length(par_rows) + 1]] <- data.frame(res_m = res, layer = paste0("z_", nm), frac_beyond_z = NA, max_abs_z = NA,
                                                     winsorised = NA, lo = m, hi = s, n_fill_before = NA, row.names = NULL)   # lo = mean, hi = sd for z_ rows
    }
    saveRDS(list(res = res, poly_coefs = cf, note = "z rows in transform_params.csv: lo = mean, hi = sd"),
            file.path(d, "transform_params.rds"))
    na_centre <- sapply(Zc, function(x) sum(is.na(x)))
    cat("\nNA in centre cells, z-layers:", paste(names(na_centre), na_centre, collapse = "; "), "\n")
    lg("stack", sprintf("%d m: fill, transform, winsorise rule, poly(elev,2), standardise, saved", res))
    
    # ---- design matrices ----
    lv <- values(rd("lith5"), mat = FALSE)[ic]; gv <- values(rd("geomorph3"), mat = FALSE)[ic]
    D <- data.frame(Zc, lith_basalt = as.numeric(lv == 1), lith_sandstone = as.numeric(lv == 3),
                    lith_carbonate = as.numeric(lv == 4), lith_other = as.numeric(lv == 5),
                    geom_dissected = as.numeric(gv == 2), geom_alluvial = as.numeric(gv == 3))
    cat("Reference levels: lithology = crystalline_felsic, geomorphology = pediplain\n")
    
    # ---- correlations ----
    ccols <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ltri", "z_ldist_riv3", "z_ldist_riv5")
    Rm <- cor(D[, ccols]); write.csv(round(Rm, 3), file.path(out, sprintf("corr_matrix_%dm.csv", res)))
    if (res == PRIMARY) { cat("\nCorrelation matrix, continuous layers (window-centre cells):\n"); print(round(Rm, 2)) }
    up <- which(upper.tri(Rm) & abs(Rm) > R_FLAG, arr.ind = TRUE)
    cat("Pairs with |r| >", R_FLAG, ":",
        if (nrow(up)) paste(sprintf("%s~%s (%.2f)", rownames(Rm)[up[, 1]], colnames(Rm)[up[, 2]], Rm[up]), collapse = "; ") else "none", "\n")
    
    # ---- categorical association ----
    cat("\nEta-squared (share of variance in continuous layer explained by category):\n")
    et <- data.frame(layer = ccols,
                     lith5 = sapply(ccols, function(c) eta2(D[[c]], lv)),
                     geomorph3 = sapply(ccols, function(c) eta2(D[[c]], gv)), row.names = NULL)
    print(cbind(et[1], round(et[-1], 3)), row.names = FALSE)
    tb <- table(lith5 = lv, geomorph3 = gv)
    chi <- suppressWarnings(chisq.test(tb)$statistic); V <- sqrt(chi / (sum(tb) * (min(dim(tb)) - 1)))
    cat("Cramer's V, lithology x geomorphology:", round(as.numeric(V), 3), "(descriptive; cells are spatially dependent)\n")
    if (res == PRIMARY) { cat("Cross-tabulation, % of window-centre cells:\n"); print(round(100 * tb / sum(tb), 2)) }
    cat_rows[[length(cat_rows) + 1]] <- cbind(res_m = res, et, cramers_V = as.numeric(V))
    
    # ---- VIF by set + drop rule ----
    cat("\nVIF by model set:\n")
    for (s in names(sets)) {
      M <- as.matrix(D[, sets[[s]]]); v <- vif_of(M)
      ar <- apply_rule(M)
      cat(sprintf("%s: max VIF %.2f (%s)", s, max(v), names(which.max(v))))
      cat(sprintf(" | drop rule: %s\n", if (length(ar$dropped)) paste("drop", paste(sub("^z_", "", ar$dropped), collapse = ", "),
                                                                      sprintf("-> max VIF %.2f", max(ar$vif))) else "no drop needed"))
      vif_rows[[length(vif_rows) + 1]] <- data.frame(res_m = res, set = s, variable = names(v), vif = round(as.numeric(v), 3), row.names = NULL)
      drop_rows[[length(drop_rows) + 1]] <- data.frame(res_m = res, set = s, max_vif_before = round(max(v), 3),
                                                       dropped = paste(ar$dropped, collapse = ";"), max_vif_after = round(max(ar$vif), 3), row.names = NULL)
      if (res == PRIMARY) print(round(v, 2))
    }
  }
  
  write.csv(do.call(rbind, par_rows), file.path(out, "transform_params.csv"), row.names = FALSE)
  write.csv(do.call(rbind, vif_rows), file.path(out, "vif_by_set.csv"), row.names = FALSE)
  write.csv(do.call(rbind, drop_rows), file.path(out, "drop_rule_result.csv"), row.names = FALSE)
  write.csv(do.call(rbind, cat_rows), file.path(out, "categorical_association.csv"), row.names = FALSE)
  
  if (!is.null(defs)) {
    defs$transform[defs$covariate == "elev"]   <- "orthogonal degree-2 polynomial over window-centre cells, z-score"
    defs$transform[defs$covariate == "slope"]  <- "log(1+slope_deg), z-score"
    defs$transform[defs$covariate == "tri"]    <- "log(1+TRI_m), z-score"
    write.csv(defs, defs_f, row.names = FALSE)
    cat("\ncovariate_definitions.csv updated (elev, slope, tri transforms)\n")
  }
  write.csv(do.call(rbind, logrows), file.path(out, "processing_log_11.csv"), row.names = FALSE)
  cat("\nSaved z-layers to data_derived/env/<res>m/ ; tables in 03_environment/outputs\n")
  
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_11.txt"))
cat("\nDONE. Outputs in:", out, "\n")