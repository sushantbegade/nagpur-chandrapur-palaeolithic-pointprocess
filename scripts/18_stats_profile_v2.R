# =============================================================================
# 18_stats_profile_v2.R
# Stage: 10_simulation_recovery | cost reduction, corrected re-run of script 17 B-F
# Purpose: (A) pcfinhom bandwidth semantics; (B) leave-one-out vs plain intensity,
#   timing and curve differences; (C) quadrature spacing; (D) Geyer fit vs quadrature;
#   (E) rmh fixed overhead and convergence; (F) corrected pipeline timing + parallel.
# SYNTHETIC patterns only. Requires R/recovery_setup.R (written by script 17).
# Outputs: 10_simulation_recovery/outputs/stats_profile_v2_report.txt,
#   stats_profile_v2.csv, rmh_convergence.csv, sessionInfo_18.txt
# Not executed by assistant. Run top to bottom in RStudio (15-25 min).
# =============================================================================

REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
MULT <- 0.175
G_ARGS <- list(kernel = "epanechnikov", stoyan = 0.15, zerocor = "convolution", divisor = "r", renormalise = TRUE, normpower = 2)
L_ARGS <- list(renormalise = TRUE, normpower = 1)

while (sink.number() > 0) sink()
source(REPO_CONFIG)
source(file.path(REPO_DIR, "R", "recovery_setup.R"))
out <- dir_of("10_simulation_recovery"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(width = 200)
maxdiff <- function(a, b) suppressWarnings(max(abs(a - b), na.rm = TRUE))
tab <- list(); add <- function(step, variant, sec = NA, value = NA) tab[[length(tab) + 1]] <<- data.frame(step = step, variant = variant, sec = sec, value = value)

truth <- make_truth()
set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 1)
Xp <- to_poly(rpoispp(truth$lamW))

mk_stats <- function(loo, quad = QUAD_KM) function(Xq) {
  fit <- ppm(Xq, trend = formula_E3, covariates = covs, quad.args = list(eps = quad))
  lam <- as.numeric(fitted(fit, dataonly = TRUE, leaveoneout = loo))
  Lh <- do.call(Linhom, c(list(Xq, lambda = lam, r = RGRID, correction = "isotropic"), L_ARGS))
  gh <- do.call(pcfinhom, c(list(Xq, lambda = lam, r = RGRID, correction = "Ripley"), G_ARGS))
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(rG, 6), round(gh$r, 6))])
}

sink(file.path(out, "stats_profile_v2_report.txt"), split = TRUE)
tryCatch({
  
  cat("==== versions ====\n")
  for (p in c("spatstat.geom", "spatstat.explore", "spatstat.model", "spatstat.random", "GET")) cat(sprintf("%-18s %s\n", p, as.character(packageVersion(p))))
  cat("pinned: g_inhom", paste(names(G_ARGS), unlist(G_ARGS), sep = "=", collapse = ", "), "| L_inhom", paste(names(L_ARGS), unlist(L_ARGS), sep = "=", collapse = ", "), "\n")
  cat("synthetic Poisson pattern n =", npoints(Xp), "| intensity (per km2):", round(npoints(Xp) / area.owin(W_km), 5), "\n")
  
  cat("\n==== B0. FIT ====\n")
  t_fit <- system.time(fit <- ppm(Xp, trend = formula_E3, covariates = covs, quad.args = list(eps = QUAD_KM)))["elapsed"]
  cat(sprintf("trend-only ppm fit: %.2f s\n", t_fit)); add("ppm fit", "quad 1 km", t_fit)
  stopifnot(exists("fit"))
  
  sec("A bandwidth", {
    cat("\n==== A. pcfinhom BANDWIDTH SEMANTICS ====\n")
    lam_plain <- as.numeric(fitted(fit, dataonly = TRUE))
    g0 <- do.call(pcfinhom, c(list(Xp, lambda = lam_plain, r = RGRID, correction = "Ripley"), G_ARGS))
    cat("attributes of the result:", paste(names(attributes(g0)), collapse = ", "), "\n")
    info <- attr(g0, "info"); if (!is.null(info)) { cat("info elements:", paste(names(info), collapse = ", "), "\n")
      for (nm in intersect(c("bw", "h", "bandwidth", "stoyan"), names(info))) cat(" info$", nm, "=", paste(unlist(info[[nm]]), collapse = " "), "\n") }
    b1 <- tryCatch(as.numeric(bw.stoyan(Xp, co = 0.15)), error = function(e) { cat("bw.stoyan(co=) failed:", conditionMessage(e), "\n"); tryCatch(as.numeric(bw.stoyan(Xp)), error = function(e2) NA) })
    lam_bar <- npoints(Xp) / area.owin(W_km)
    cat(sprintf("bw.stoyan(Xp) = %.4f km | 0.15/sqrt(lambda) = %.4f km | that / sqrt(5) = %.4f km\n", b1, 0.15 / sqrt(lam_bar), 0.15 / sqrt(5 * lam_bar)))
    G2 <- G_ARGS; G2$stoyan <- NULL
    ga <- do.call(pcfinhom, c(list(Xp, lambda = lam_plain, r = RGRID, correction = "Ripley", bw = b1), G2))
    gh <- do.call(pcfinhom, c(list(Xp, lambda = lam_plain, r = RGRID, correction = "Ripley", h = b1), G2))
    gc <- function(f) col_of(f, c("iso", "trans"))
    d_bw <- maxdiff(gc(ga), gc(g0)); d_h <- maxdiff(gc(gh), gc(g0))
    cat(sprintf("default vs explicit bw = bw.stoyan: max |diff| %.5f | default vs explicit h = bw.stoyan: max |diff| %.5f\n", d_bw, d_h))
    cat("Reading: the setting with the SMALLER difference is how the default interprets the Stoyan value.\n")
    add("bandwidth test", "maxdiff default vs bw=b", NA, d_bw); add("bandwidth test", "maxdiff default vs h=b", NA, d_h)
  })
  
  sec("B intensity", {
    cat("\n==== B. FITTED INTENSITY AT DATA POINTS: leave-one-out vs plain ====\n")
    t_loo <- system.time(lam_loo <- as.numeric(fitted(fit, dataonly = TRUE, leaveoneout = TRUE)))["elapsed"]
    t_pl  <- system.time(lam_pl  <- as.numeric(fitted(fit, dataonly = TRUE, leaveoneout = FALSE)))["elapsed"]
    cat(sprintf("leave-one-out %.2f s | plain %.2f s | ratio loo/plain at data points: min %.3f median %.3f max %.3f\n", t_loo, t_pl,
                min(lam_loo / lam_pl), median(lam_loo / lam_pl), max(lam_loo / lam_pl)))
    add("fitted()", "leaveoneout", t_loo); add("fitted()", "plain", t_pl)
    cat("\nCurves from the two intensities (same pattern):\n")
    cv <- function(lam) { Lh <- do.call(Linhom, c(list(Xp, lambda = lam, r = RGRID, correction = "isotropic"), L_ARGS))
    gh <- do.call(pcfinhom, c(list(Xp, lambda = lam, r = RGRID, correction = "Ripley"), G_ARGS))
    list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(rL, 6), round(Lh$r, 6))], g = col_of(gh, c("iso", "trans"))[match(round(rG, 6), round(gh$r, 6))]) }
    tl <- system.time(c_loo <- cv(lam_loo))["elapsed"]; tp <- system.time(c_pl <- cv(lam_pl))["elapsed"]
    cat(sprintf("L_inhom + g_inhom time: loo-intensity %.2f s | plain-intensity %.2f s\n", tl, tp))
    cat(sprintf("max |dL| %.4f (L-r ranges %.2f..%.2f) | max |dg| %.4f (g ranges %.2f..%.2f)\n", maxdiff(c_loo$L, c_pl$L), min(c_loo$L), max(c_loo$L),
                maxdiff(c_loo$g, c_pl$g), min(c_loo$g), max(c_loo$g)))
    add("L/g curves", "loo intensity", tl); add("L/g curves", "plain intensity", tp); add("curve difference", "max |dL| loo vs plain", NA, maxdiff(c_loo$L, c_pl$L)); add("curve difference", "max |dg| loo vs plain", NA, maxdiff(c_loo$g, c_pl$g))
  })
  
  sec("C quadrature", {
    cat("\n==== C. QUADRATURE SPACING (Poisson trend-only fit) ====\n")
    rows <- list(); ref <- NULL; ref_ll <- NA
    for (eps in c(1, 2, 0.5)) {
      t <- system.time(f <- ppm(Xp, trend = formula_E3, covariates = covs, quad.args = list(eps = eps)))["elapsed"]
      nq <- tryCatch(npoints(f$Q$dummy) + npoints(f$Q$data), error = function(e) NA)
      if (eps == 1) { ref <- coef(f); ref_ll <- as.numeric(logLik(f)) }
      rows[[length(rows) + 1]] <- data.frame(eps_km = eps, quad_points = nq, sec = round(t, 2), max_coef_diff_vs_1km = round(maxdiff(coef(f), ref), 4), logLik = round(as.numeric(logLik(f)), 2))
    }
    print(do.call(rbind, rows), row.names = FALSE)
    for (r in rows) add("ppm quadrature", paste0("eps ", r$eps_km, " km"), r$sec, r$max_coef_diff_vs_1km)
  })
  
  sec("E rmh", {
    cat("\n==== E. rmh FIXED OVERHEAD AND CONVERGENCE (Geyer gamma 2.5, r 3, sat 3, trend multiplier", MULT, ") ====\n")
    t_mod <- system.time(mod <- rmhmodel(cif = "geyer", par = list(beta = 1, gamma = 2.5, r = 3, sat = 3), w = W_km, trend = truth$lamW * MULT))["elapsed"]
    cat(sprintf("rmhmodel construction: %.2f s\n", t_mod)); add("rmhmodel", "construct", t_mod)
    conv <- list(); Xg <- NULL
    for (nr in c(1e3, 1e5, 5e5, 2e6)) for (s in 1:5) {
      set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 400 + s)
      t <- system.time(P <- rmh(mod, start = list(n.start = 197), control = list(nrep = nr, expand = 1), verbose = FALSE))["elapsed"]
      P <- P[W_km]; if (nr == 5e5 && s == 1) Xg <- P
      conv[[length(conv) + 1]] <- data.frame(nrep = nr, rep = s, n = npoints(P), sec = round(t, 2))
    }
    conv <- do.call(rbind, conv); write.csv(conv, file.path(out, "rmh_convergence.csv"), row.names = FALSE)
    ag <- do.call(rbind, lapply(split(conv, conv$nrep), function(d) data.frame(nrep = d$nrep[1], mean_n = round(mean(d$n), 1), sd_n = round(sd(d$n), 1), mean_sec = round(mean(d$sec), 2))))
    print(ag, row.names = FALSE)
    for (i in seq_len(nrow(ag))) add("rmh", paste("nrep", ag$nrep[i]), ag$mean_sec[i], ag$mean_n[i])
    set.seed(SEEDS$master + SEEDS$offsets["recovery"] + 500)
    t_multi <- system.time(PM <- rmh(mod, start = list(n.start = 197), control = list(nrep = 5e5, expand = 1), nsim = 3, verbose = FALSE))["elapsed"]
    cat(sprintf("one call with nsim = 3 at 5e5 steps: %.2f s (compare 3 separate calls above)\n", t_multi)); add("rmh", "nsim=3 single call", t_multi)
    
    cat("\n==== D. GEYER FIT vs QUADRATURE (pattern from the 5e5-step chain, n =", npoints(Xg), ") ====\n")
    rows <- list()
    for (eps in c(2, 1, 0.5)) {
      t <- system.time(f <- ppm(Xg, trend = formula_E3, interaction = Geyer(r = 4, sat = 4), covariates = covs, rbord = 10, quad.args = list(eps = eps)))["elapsed"]
      cf <- coef(f); ii <- grep("Interaction", names(cf))
      rows[[length(rows) + 1]] <- data.frame(eps_km = eps, sec = round(t, 2), gamma_hat = round(exp(cf[ii[1]]), 3), logPL = round(as.numeric(logLik(f)), 2))
    }
    print(do.call(rbind, rows), row.names = FALSE)
    for (r in rows) add("Geyer fit", paste0("eps ", r$eps_km, " km, r4 sat4 rbord10"), r$sec, r$gamma_hat)
  })
  
  sec("F pipeline", {
    cat("\n==== F. PIPELINE TIMING AND PARALLEL SPEED-UP ====\n")
    patterns <- lapply(1:8, function(i) { set.seed(SEEDS$master + SEEDS$offsets["bootstrap"] + 100 + i); to_poly(simulate(fit, nsim = 1, drop = TRUE, verbose = FALSE)) })
    s_loo <- mk_stats(TRUE); s_pl <- mk_stats(FALSE)
    t1 <- system.time(r_loo <- lapply(patterns[1:2], s_loo))["elapsed"] / 2
    t2 <- system.time(r_pl  <- lapply(patterns[1:2], s_pl))["elapsed"] / 2
    cat(sprintf("per pattern, full pipeline (fit + intensity + L + g): leave-one-out %.2f s | plain %.2f s\n", t1, t2))
    cat(sprintf("curve difference loo vs plain over 2 bootstrap patterns: max |dL| %.4f | max |dg| %.4f\n",
                max(sapply(1:2, function(i) maxdiff(r_loo[[i]]$L, r_pl[[i]]$L))), max(sapply(1:2, function(i) maxdiff(r_loo[[i]]$g, r_pl[[i]]$g)))))
    add("pipeline", "loo per pattern", t1); add("pipeline", "plain per pattern", t2)
    ncl <- min(4, parallel::detectCores(logical = FALSE))
    ts <- system.time(ser <- lapply(patterns, s_pl))["elapsed"]
    cl <- parallel::makeCluster(ncl)
    tp <- tryCatch({
      parallel::clusterExport(cl, c("covs", "W_km", "formula_E3", "QUAD_KM", "RGRID", "rL", "rG", "col_of", "G_ARGS", "L_ARGS", "mk_stats"), envir = environment())
      parallel::clusterEvalQ(cl, suppressPackageStartupMessages({ library(spatstat.geom); library(spatstat.explore); library(spatstat.model); library(spatstat.random) }))
      t0 <- system.time(pr <- parallel::parLapply(cl, patterns, mk_stats(FALSE)))["elapsed"]
      list(t = t0, d = max(sapply(seq_along(ser), function(i) maxdiff(ser[[i]]$g, pr[[i]]$g))))
    }, error = function(e) { cat("parallel FAILED:", conditionMessage(e), "\n"); NULL }, finally = parallel::stopCluster(cl))
    cat(sprintf("8 patterns, plain pipeline: serial %.1f s | %d workers", ts, ncl))
    if (!is.null(tp)) cat(sprintf(" %.1f s | speed-up x%.2f | max difference vs serial %.2g\n", tp$t, ts / tp$t, tp$d)) else cat("\n")
    add("8 patterns", "serial plain", ts); if (!is.null(tp)) add("8 patterns", paste0("parallel x", ncl), tp$t)
  })
  
  write.csv(do.call(rbind, tab), file.path(out, "stats_profile_v2.csv"), row.names = FALSE)
}, error = function(e) cat("\nSCRIPT ERROR:", conditionMessage(e), "\n"), finally = sink())

writeLines(capture.output(sessionInfo()), file.path(out, "sessionInfo_18.txt"))
cat("\nDONE. Outputs in:", out, "\n")