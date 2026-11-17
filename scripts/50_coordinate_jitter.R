# =====================================================================
# 50_coordinate_jitter.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW; writes R/jitter_tasks.R)
# Purpose : coordinate-uncertainty robustness of the frozen M1 deficiency. Primary scenario:
#           Gaussian jitter of the best-supported coordinates by source class (FIELD 10 m,
#           AUTHOR_REPORTED 25 m, GIS_MAP 250 m), 200 jitter replicates plus 100 nominal-
#           coordinate baseline replicates (distinct seeds), full M1 test per replicate at
#           nsim 199 (E3, joint ERL L_inhom + g_inhom, nstep 2), 4 workers, batches of 20.
#           Draws outside the window or in NA-covariate pixels are redrawn (<= 50 tries;
#           fallbacks to the nominal position are counted).
# Outputs : outputs/50_report.txt, 50_batch_NN.csv, 50_jitter_replicates.csv,
#           50_jitter_summary.csv, 50_sessionInfo.txt; R/jitter_tasks.R
# Notes   : RUN ONCE (resumable by batch). Retention = share of replicates with p < 0.05;
#           frozen criterion >= 0.80, read against the baseline retention at the same nsim
#           (Monte Carlo noise alone keeps retention below 100%). Precision rule does not
#           apply to replicates. Restart R before running. Workers load the cache and the
#           generic files only (no terra).
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
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
f_res <- file.path(OUT_DIR, "50_jitter_summary.csv")
if (file.exists(f_res)) stop("50_jitter_summary.csv exists: this analysis runs once. Stop; do not rerun.")
GENERIC_FILE <- file.path(REPO_DIR, "R", "robust_generic.R"); stopifnot(file.exists(GENERIC_FILE))
JIT_FILE <- file.path(REPO_DIR, "R", "jitter_tasks.R")
N_J <- 200; N_B <- 100; NSIM <- TESTS$nsim_recovery; BATCH <- 20; N_WORK <- 4
SEED_J <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 70000

# ---- open inputs BEFORE sink ----
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
SITES_G <- read.csv(f_sites, stringsAsFactors = FALSE)
stopifnot(nrow(SITES_G) == 197, max(abs(X_obs_km$x - SITES_G$x_m / 1000), abs(X_obs_km$y - SITES_G$y_m / 1000)) < 1e-3,
          all(SITES_G$coordinate_source %in% names(JITTER_SIGMA_M)))
source(GENERIC_FILE)

# ---- jitter functions (dumped for workers) ----
pairs_at <- function(P, r) sum(closepaircounts(P, r)) / 2
jitter_obs <- function(seed) {
  set.seed(seed)
  sig <- unname(JITTER_SIGMA_M[SITES_G$coordinate_source]) / 1000
  x <- X_obs_km$x; y <- X_obs_km$y; fb <- 0L; disp <- numeric(length(x))
  for (i in seq_along(x)) {
    ok <- FALSE
    for (tr in seq_len(50)) {
      xn <- x[i] + rnorm(1, 0, sig[i]); yn <- y[i] + rnorm(1, 0, sig[i])
      if (inside.owin(xn, yn, W_km)) {
        z <- vapply(covs, function(im) lookup.im(im, xn, yn, naok = TRUE), numeric(1))
        if (!anyNA(z)) { disp[i] <- sqrt((xn - x[i])^2 + (yn - y[i])^2); x[i] <- xn; y[i] <- yn; ok <- TRUE; break }
      }
    }
    if (!ok) fb <- fb + 1L
  }
  list(X = ppp(x, y, window = W_km, check = FALSE), fallbacks = fb, disp = disp)
}
jit_test <- function(X, nsim, seed) {
  S <- mk("J", X, W_km, covs, DUMMY_FIXED)
  ft <- rob_fit(X, S); ob <- rob_curves(X, ft, S)
  lam <- exp(Reduce(`+`, Map(function(a, b) a * b, S$covW, ft$beta[-1])) + ft$beta[1])
  set.seed(seed)
  sL <- matrix(NA_real_, length(S$rL), nsim); sG <- matrix(NA_real_, length(S$rG), nsim)
  pr <- matrix(NA_real_, 3, nsim); nn <- rep(NA_real_, nsim)
  for (k in seq_len(nsim)) {
    res <- tryCatch({
      Xs <- rob_poly(rpoispp(lam), S); Xc <- rob_clean(Xs, S); fs <- rob_fit(Xc, S)
      list(cs = rob_curves(Xc, fs, S), n = npoints(Xc), pr = vapply(c(2, 3, 5), function(r) pairs_at(Xc, r), numeric(1)))
    }, error = function(e) NULL)
    if (!is.null(res)) { sL[, k] <- res$cs$L; sG[, k] <- res$cs$g; nn[k] <- res$n; pr[, k] <- res$pr }
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = S$rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = S$rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc), n_sim = mean(nn[okc]),
       pobs = vapply(c(2, 3, 5), function(r) pairs_at(X, r), numeric(1)), psim = rowMeans(pr[, okc, drop = FALSE]))
}
jit_task <- function(job) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch({
    if (job$type == "jitter") { jj <- jitter_obs(job$seed_j); X <- jj$X; fb <- jj$fallbacks; md <- median(jj$disp) } else { X <- X_obs_km; fb <- 0L; md <- 0 }
    r <- jit_test(X, job$nsim, job$seed_t)
    list(n = npoints(X), p = r$p, usable = r$usable, n_sim = r$n_sim, po2 = r$pobs[1], po3 = r$pobs[2], po5 = r$pobs[3],
         ps2 = r$psim[1], ps3 = r$psim[2], ps5 = r$psim[3], fallbacks = fb, med_disp_km = md, err = NA_character_)
  }, error = function(e) list(err = conditionMessage(e)))
  c(list(type = job$type, i = job$i, sec = unname(proc.time()["elapsed"]) - t0), out)
}
jrow_df <- function(L) do.call(rbind, lapply(L, function(x) {
  g <- function(f) if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  data.frame(type = x$type, i = x$i, sec = x$sec, n = g("n"), p = g("p"), usable = g("usable"), n_sim = g("n_sim"),
             po2 = g("po2"), po3 = g("po3"), po5 = g("po5"), ps2 = g("ps2"), ps3 = g("ps3"), ps5 = g("ps5"),
             fallbacks = g("fallbacks"), med_disp_km = g("med_disp_km"),
             err = if (is.null(x$err)) NA_character_ else as.character(x$err), stringsAsFactors = FALSE)
}))
init_worker4 <- function(repo, cache, generic, jit, sites_path) {
  assign("REPO_DIR", repo, envir = globalenv()); Sys.setenv(RECOVERY_CACHE = cache)
  source(file.path(repo, "R", "recovery_worker.R"), local = globalenv())
  source(file.path(repo, "R", "recovery_functions.R"), local = globalenv())
  assign("SITES_G", read.csv(sites_path, stringsAsFactors = FALSE), envir = globalenv())
  source(generic, local = globalenv())
  source(jit, local = globalenv())
  TRUE
}
dump(c("pairs_at", "jitter_obs", "jit_test", "jit_task", "jrow_df"), file = JIT_FILE)
wilson <- function(k, n) as.numeric(suppressWarnings(prop.test(k, n, correct = FALSE)$conf.int))

REPORT <- file.path(OUT_DIR, "50_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("50_coordinate_jitter | ", format(Sys.time()), "\n", sep = "")
  cat(sprintf("replicates: jitter %d + baseline %d | nsim %d | batch %d | workers %d | seed base %d\n", N_J, N_B, NSIM, BATCH, N_WORK, SEED_J))
  cat("sigma (m): ", paste(sprintf("%s %g", names(JITTER_SIGMA_M), JITTER_SIGMA_M), collapse = ", "), " | classes: ",
      paste(sprintf("%s %d", names(table(SITES_G$coordinate_source)), as.integer(table(SITES_G$coordinate_source))), collapse = ", "), "\n", sep = "")
  
  cat("\n== A. Jitter geometry (200 main-session draws) ==\n")
  D <- lapply(seq_len(200), function(i) jitter_obs(SEED_J + i))
  disp <- do.call(rbind, lapply(D, function(d) d$disp)); cls <- SITES_G$coordinate_source
  GE <- do.call(rbind, lapply(names(JITTER_SIGMA_M), function(cl) data.frame(class = cl, sites = sum(cls == cl), sigma_m = JITTER_SIGMA_M[[cl]],
                                                                             expected_median_m = 1.1774 * JITTER_SIGMA_M[[cl]], observed_median_m = 1000 * median(disp[, cls == cl]),
                                                                             observed_p95_m = 1000 * quantile(disp[, cls == cl], 0.95))))
  print(GE, digits = 4, row.names = FALSE)
  cat(sprintf("fallbacks to the nominal position: %d sites in %d of 200 draws\n", sum(vapply(D, function(d) d$fallbacks, 1L)),
              sum(vapply(D, function(d) d$fallbacks > 0, logical(1)))))
  
  cl <- makeCluster(N_WORK)
  ini <- clusterCall(cl, init_worker4, REPO_DIR, CACHE_FILE, GENERIC_FILE, JIT_FILE, f_sites)
  cat("\nworkers initialised:", all(unlist(ini)), "\n")
  ex <- clusterCall(cl, function() c(exists("jit_task"), exists("SITES_G"), "terra" %in% loadedNamespaces()))
  cat("worker check (jit_task, SITES_G present; terra loaded):", paste(ex[[1]], collapse = " "), "\n")
  
  cat("\n== B. Worker smoke (2 jitter + 2 baseline, nsim 19; completion check only) ==\n")
  sm <- parLapply(cl, list(list(type = "jitter", i = 901, nsim = 19, seed_j = SEED_J + 901, seed_t = SEED_J + 5901),
                           list(type = "jitter", i = 902, nsim = 19, seed_j = SEED_J + 902, seed_t = SEED_J + 5902),
                           list(type = "base", i = 901, nsim = 19, seed_j = NA, seed_t = SEED_J + 6901),
                           list(type = "base", i = 902, nsim = 19, seed_j = NA, seed_t = SEED_J + 6902)), jit_task)
  smd <- jrow_df(sm); print(smd[, c("type", "i", "n", "usable", "fallbacks", "med_disp_km", "sec", "err")], digits = 3, row.names = FALSE)
  if (any(!is.na(smd$err)) || any(is.na(smd$p))) stop("worker smoke failed; paste output")
  
  cat("\n== C. Batches ==\n")
  key <- c((seq_len(N_J)) / N_J, (seq_len(N_B)) / N_B)
  typ <- c(rep("jitter", N_J), rep("base", N_B)); idx <- c(seq_len(N_J), seq_len(N_B))
  ord <- order(key, typ)
  JOBS <- lapply(ord, function(o) list(type = typ[o], i = idx[o], nsim = NSIM,
                                       seed_j = SEED_J + 1000 + idx[o],
                                       seed_t = SEED_J + (if (typ[o] == "jitter") 500000 else 700000) + idx[o]))
  nb <- ceiling(length(JOBS) / BATCH)
  for (b in seq_len(nb)) {
    f_b <- file.path(OUT_DIR, sprintf("50_batch_%02d.csv", b))
    if (file.exists(f_b)) { cat(sprintf("  batch %2d: found, skipped\n", b)); next }
    grp <- JOBS[((b - 1) * BATCH + 1):min(b * BATCH, length(JOBS))]; t0 <- unname(proc.time()["elapsed"])
    res <- parLapplyLB(cl, grp, jit_task); db <- jrow_df(res); db$batch <- b
    write.csv(db, f_b, row.names = FALSE)
    cat(sprintf("  batch %2d of %d: %d replicates, %.0f s wall, failed %d [%s]\n", b, nb, nrow(db), unname(proc.time()["elapsed"]) - t0,
                sum(!is.na(db$err)), format(Sys.time(), "%H:%M")))
  }
  
  cat("\n== D. Retention ==\n")
  fs <- file.path(OUT_DIR, sprintf("50_batch_%02d.csv", seq_len(nb)))
  ALL <- do.call(rbind, lapply(fs[file.exists(fs)], read.csv, stringsAsFactors = FALSE))
  write.csv(ALL, file.path(OUT_DIR, "50_jitter_replicates.csv"), row.names = FALSE)
  SUMR <- data.frame()
  for (ty in c("jitter", "base")) {
    d <- ALL[ALL$type == ty & !is.na(ALL$p), ]; k <- sum(d$p < TESTS$alpha); n <- nrow(d); ci <- wilson(k, n)
    q <- quantile(d$p, c(0.1, 0.25, 0.5, 0.75, 0.9))
    SUMR <- rbind(SUMR, data.frame(type = ty, reps = n, failed = sum(ALL$type == ty) - n, retained = k, retention = k / n, ci_lo = ci[1], ci_hi = ci[2],
                                   p_q10 = q[1], p_q25 = q[2], p_median = q[3], p_q75 = q[4], p_q90 = q[5], min_usable = min(d$usable),
                                   pair_ratio_2 = median(d$po2 / d$ps2), pair_ratio_3 = median(d$po3 / d$ps3), pair_ratio_5 = median(d$po5 / d$ps5),
                                   mean_fallbacks = mean(d$fallbacks), mean_median_disp_m = 1000 * mean(d$med_disp_km)))
  }
  print(SUMR, digits = 3, row.names = FALSE)
  rj <- SUMR[SUMR$type == "jitter", ]; rb <- SUMR[SUMR$type == "base", ]
  dif <- rj$retention - rb$retention; se <- sqrt(rj$retention * (1 - rj$retention) / rj$reps + rb$retention * (1 - rb$retention) / rb$reps)
  cat(sprintf("\njitter retention %.3f (Wilson %.3f-%.3f) | baseline retention %.3f (Wilson %.3f-%.3f) | difference %.3f (approx 95%% %.3f to %.3f)\n",
              rj$retention, rj$ci_lo, rj$ci_hi, rb$retention, rb$ci_lo, rb$ci_hi, dif, dif - 1.96 * se, dif + 1.96 * se))
  verdict <- if (rj$retention >= 0.80) "frozen criterion MET (jitter retention >= 0.80)" else "frozen criterion NOT MET (jitter retention < 0.80)"
  cat("verdict:", verdict, "| baseline retention", sprintf("%.3f", rb$retention), "shows the ceiling at nsim", NSIM, "\n")
  SUMR$verdict <- verdict; SUMR$date <- format(Sys.time())
  write.csv(SUMR, f_res, row.names = FALSE)
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink(); if (exists("cl")) try(stopCluster(cl), silent = TRUE) })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "50_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")