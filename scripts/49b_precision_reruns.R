# =====================================================================
# 49b_precision_reruns.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW; writes R/robust_generic.R)
# Purpose : precision rule (declared 2026-11-12): robustness settings with p in [0.04, 0.06]
#           at nsim 999 (script 49: R02_datasetC 0.045, R12_quad05 0.045) are recomputed at
#           nsim 10000 with a distinct seed; the 10000-simulation p is used; threshold 0.05.
#           Simulation chunks (10 x 1000 per setting) run on 4 workers; observed curves in the
#           main session; same generic code as script 49. A worker-vs-main identity check
#           (50 simulations, same seed) must pass first.
# Outputs : outputs/49b_report.txt, 49b_precision_result.csv, 49b_chunk_<setting>_<NN>.rds,
#           49b_sessionInfo.txt; R/robust_generic.R
# Notes   : RUN ONCE (resumable by chunk). Rule: p(10000) < 0.05 = deficiency retained;
#           otherwise not retained at nsim 10000 (R02, n = 97: inconclusive-low-n caveat).
#           Restart R before running. Workers load the cache and the generic file only (no terra).
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))
suppressPackageStartupMessages(library(parallel))

OUT_DIR <- dir_of("11_generative_gof", "outputs")
f_res <- file.path(OUT_DIR, "49b_precision_result.csv")
if (file.exists(f_res)) stop("49b_precision_result.csv exists: this analysis runs once. Stop; do not rerun.")
GENERIC_FILE <- file.path(REPO_DIR, "R", "robust_generic.R")
N_WORK <- 4; NCH <- 10; CH <- 1000
SEED_P <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 60000
GMIN_FLOOR <- GINHOM_LOWER_KM[["A"]]
RUN <- c("R12_quad05", "R02_datasetC")

# ---- open inputs BEFORE sink ----
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
SITES_G <- read.csv(f_sites, stringsAsFactors = FALSE)
f_s49 <- file.path(OUT_DIR, "49_robustness_summary.csv"); stopifnot(file.exists(f_s49))
S49 <- read.csv(f_s49, stringsAsFactors = FALSE)
flagged <- S49$setting[S49$flag == "PRECISION"]
stopifnot(setequal(flagged, RUN))

# ---- generic functions (dumped for workers) ----
mk <- function(name, X, W, cv, dum, Lmin = 1, corrL = "isotropic", corrG = "Ripley") {
  n <- npoints(X); a <- area.owin(W); half <- PCF_STOYAN_C / sqrt(n / a); gmin <- max(Lmin, GMIN_FLOOR, ceiling(half * 2) / 2)
  list(name = name, X = X, W = W, cov = cv, n = n, area = a, bw = half / sqrt(5),
       rL = RGRID[RGRID >= Lmin & RGRID <= R_MAX_DECISION_KM], rG = RGRID[RGRID >= gmin & RGRID <= R_MAX_DECISION_KM],
       corrL = corrL, corrG = corrG, dum = dum, covW = lapply(cv, function(im) im[W, drop = FALSE]))
}
cov_na <- function(P, S) {
  if (npoints(P) == 0) return(logical(0))
  Z <- do.call(cbind, lapply(S$cov, function(im) lookup.im(im, P$x, P$y, naok = TRUE)))
  rowSums(is.na(Z)) > 0
}
rob_clean <- function(P, S) { b <- cov_na(P, S); if (any(b)) P[!b] else P }
rob_poly <- function(P, S) { Q <- P[S$W]; if (Window(Q)$type != "polygonal") stop("not polygonal"); Q }
rob_fit <- function(Xp, S) {
  Q <- quadscheme(Xp, dummy = S$dum, method = "grid")
  w <- w.quad(Q); isd <- is.data(Q); xq <- x.quad(Q); yq <- y.quad(Q)
  Z <- vapply(S$cov, function(im) lookup.im(im, xq, yq, naok = TRUE), numeric(length(xq)))
  ok <- stats::complete.cases(Z) & w > 0
  X1 <- cbind(`(Intercept)` = 1, Z)
  g <- suppressWarnings(stats::glm.fit(X1[ok, , drop = FALSE], (isd / w)[ok], weights = w[ok], family = poisson(), intercept = FALSE))
  Xd <- X1[isd & ok, , drop = FALSE]
  if (nrow(Xd) != npoints(Xp)) stop("data points dropped in rob_fit: ", nrow(Xd), " of ", npoints(Xp))
  list(beta = g$coefficients, lam_data = as.numeric(exp(Xd %*% g$coefficients)))
}
rob_curves <- function(Xp, ft, S) {
  Lh <- do.call(Linhom, c(list(Xp, lambda = ft$lam_data, r = RGRID, correction = S$corrL), L_ARGS))
  gh <- do.call(pcfinhom, c(list(Xp, lambda = ft$lam_data, r = RGRID, correction = S$corrG, bw = S$bw), G_ARGS_NS))
  list(L = (col_of(Lh, c("iso", "trans", "border")) - Lh$r)[match(round(S$rL, 6), round(Lh$r, 6))],
       g = col_of(gh, c("iso", "trans"))[match(round(S$rG, 6), round(gh$r, 6))])
}
make_setting <- function(nm) {
  if (nm == "R02_datasetC") return(mk(nm, X_obs_km[SITES_G$dataset_C], W_km, covs, DUMMY_FIXED, Lmin = 2))
  if (nm == "R12_quad05") return(mk(nm, X_obs_km, W_km, covs, default.dummy(X_obs_km, eps = 0.5)))
  stop("unknown setting ", nm)
}
chunk_run <- function(S, nsim, seed) {
  set.seed(seed)
  ft <- rob_fit(S$X, S)
  lam <- exp(Reduce(`+`, Map(function(a, b) a * b, S$covW, ft$beta[-1])) + ft$beta[1])
  sL <- matrix(NA_real_, length(S$rL), nsim); sG <- matrix(NA_real_, length(S$rG), nsim)
  for (k in seq_len(nsim)) {
    res <- tryCatch({ Xs <- rob_poly(rpoispp(lam), S); Xc <- rob_clean(Xs, S); fs <- rob_fit(Xc, S); rob_curves(Xc, fs, S) }, error = function(e) NULL)
    if (!is.null(res)) { sL[, k] <- res$L; sG[, k] <- res$g }
  }
  list(sL = sL, sG = sG)
}
chunk_task <- function(job) {
  t0 <- unname(proc.time()["elapsed"])
  S <- make_setting(job$setting); r <- chunk_run(S, job$nsim, job$seed)
  c(list(setting = job$setting, ch = job$ch, sec = unname(proc.time()["elapsed"]) - t0), r)
}
init_worker3 <- function(repo, cache, generic, sites_path) {
  assign("REPO_DIR", repo, envir = globalenv()); Sys.setenv(RECOVERY_CACHE = cache)
  source(file.path(repo, "R", "recovery_worker.R"), local = globalenv())
  source(file.path(repo, "R", "recovery_functions.R"), local = globalenv())
  source(generic, local = globalenv())
  assign("SITES_G", read.csv(sites_path, stringsAsFactors = FALSE), envir = globalenv())
  TRUE
}
dump(c("GMIN_FLOOR", "mk", "cov_na", "rob_clean", "rob_poly", "rob_fit", "rob_curves", "make_setting", "chunk_run", "chunk_task"), file = GENERIC_FILE)

REPORT <- file.path(OUT_DIR, "49b_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("49b_precision_reruns | ", format(Sys.time()), "\n", sep = "")
  cat(sprintf("settings %s | nsim %d (%d chunks of %d) | workers %d | seed base %d\n", paste(RUN, collapse = ", "), NCH * CH, NCH, CH, N_WORK, SEED_P))
  cat("flagged in script 49:", paste(flagged, collapse = ", "), "| p(999):", paste(sprintf("%s %.4f", S49$setting[S49$flag == "PRECISION"], S49$p[S49$flag == "PRECISION"]), collapse = ", "), "\n")
  cl <- makeCluster(N_WORK)
  ini <- clusterCall(cl, init_worker3, REPO_DIR, CACHE_FILE, GENERIC_FILE, f_sites)
  cat("workers initialised:", all(unlist(ini)), "\n")
  ex <- clusterCall(cl, function() c(exists("chunk_task"), exists("SITES_G"), "terra" %in% loadedNamespaces()))
  cat("worker check (chunk_task, SITES_G present; terra loaded):", paste(ex[[1]], collapse = " "), "\n")
  
  cat("\n== Worker-vs-main identity check (R02, 50 simulations, same seed) ==\n")
  S02 <- make_setting("R02_datasetC")
  m <- chunk_run(S02, 50, SEED_P + 999)
  w <- parLapply(cl[1], list(list(setting = "R02_datasetC", ch = 0, nsim = 50, seed = SEED_P + 999)), chunk_task)[[1]]
  d <- max(abs(m$sL - w$sL), abs(m$sG - w$sG), na.rm = TRUE); same_na <- identical(is.na(m$sL), is.na(w$sL)) && identical(is.na(m$sG), is.na(w$sG))
  cat(sprintf("max |difference| %.2e | NA pattern identical %s | worker time %.0f s\n", d, same_na, w$sec))
  if (!(d < 1e-10 && same_na)) stop("worker path differs from main path; stop and paste output")
  
  cat("\n== Chunks ==\n")
  jobs <- list()
  for (nm in RUN) for (ch in seq_len(NCH)) {
    f_c <- file.path(OUT_DIR, sprintf("49b_chunk_%s_%02d.rds", nm, ch))
    if (!file.exists(f_c)) jobs[[length(jobs) + 1]] <- list(setting = nm, ch = ch, nsim = CH, seed = SEED_P + 1000 * match(nm, RUN) + ch, file = f_c)
  }
  cat("chunks to run:", length(jobs), "\n")
  i <- 1
  while (i <= length(jobs)) {
    grp <- jobs[i:min(i + N_WORK - 1, length(jobs))]; t0 <- unname(proc.time()["elapsed"])
    res <- parLapplyLB(cl, grp, chunk_task)
    for (k in seq_along(res)) saveRDS(list(sL = res[[k]]$sL, sG = res[[k]]$sG), grp[[k]]$file)
    cat(sprintf("  round done: %s | %.0f s wall [%s]\n", paste(sprintf("%s#%d", vapply(grp, function(j) j$setting, ""), vapply(grp, function(j) j$ch, 1)), collapse = " "),
                unname(proc.time()["elapsed"]) - t0, format(Sys.time(), "%H:%M")))
    i <- i + N_WORK
  }
  
  cat("\n== Results ==\n")
  OUT <- data.frame()
  for (nm in RUN) {
    S <- make_setting(nm); ft <- rob_fit(S$X, S); ob <- rob_curves(S$X, ft, S)
    CK <- lapply(sprintf("49b_chunk_%s_%02d.rds", nm, seq_len(NCH)), function(f) readRDS(file.path(OUT_DIR, f)))
    sL <- do.call(cbind, lapply(CK, function(x) x$sL)); sG <- do.call(cbind, lapply(CK, function(x) x$sG))
    okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
    csL <- GET::create_curve_set(list(r = S$rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
    csG <- GET::create_curve_set(list(r = S$rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
    p <- erl_p(csL, csG, 2); N <- sum(okc); se <- sqrt(p * (1 - p) / N)
    dec <- if (is.na(p)) "UNDETERMINED" else if (p < TESTS$alpha) "DEFICIENCY retained (nsim 10000)" else
      if (S$n <= 114) "NOT retained at nsim 10000 (n <= 114: inconclusive)" else "NOT retained at nsim 10000"
    p999 <- S49$p[S49$setting == nm]
    cat(sprintf("  %-14s n %3d | usable %d of %d | p(10000) %.5f (SE %.4f; 95%% %.4f-%.4f) | p(999) %.4f | %s\n",
                nm, S$n, N, NCH * CH, p, se, p - 1.96 * se, p + 1.96 * se, p999, dec))
    OUT <- rbind(OUT, data.frame(setting = nm, n = S$n, nsim = N, p_10000 = p, mc_se = se, p_999 = p999, decision = dec, date = format(Sys.time())))
  }
  write.csv(OUT, f_res, row.names = FALSE)
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink(); if (exists("cl")) try(stopCluster(cl), silent = TRUE) })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "49b_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")