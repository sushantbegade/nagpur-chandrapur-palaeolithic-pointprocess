# =====================================================================
# 56_block_validation.R
# Stage   : 12_validation
# Folder  : 12_validation/scripts/   (NEW; creates 12_validation/outputs)
# Purpose : buffered spatial block validation of the FIRST-ORDER environmental model
#           (frozen: blocks 60 km, 5 folds, buffer 10 km). Held-out Poisson log-likelihood
#           gain per test site over a homogeneous intensity fitted on the same training region,
#           for E1, E2, E3 (primary), E4, over 50 random block->fold assignments.
#           Pixel-level Poisson regression on the 500 m grid (exposure 0.25 km2 per pixel);
#           buffer = distance to the nearest test-block rectangle; predicted intensity floored at
#           1% of the homogeneous intensity (unfloored value also reported); all-zero training
#           columns dropped and flagged; separation flag |coef| > 8.
# Outputs : outputs/56_report.txt, 56_blockcv.csv, 56_sessionInfo.txt,
#           protocol_deviation_log_before_56.csv; docs/protocol_deviation_log.csv (1 row)
# Notes   : Validation answers ONLY whether first-order environmental structure transfers across
#           space; it is not evidence about second-order structure, Gibbs interaction or
#           the M1 deficiency. Reading rules (set before running): TRANSFERS if the pooled
#           floored gain per held-out site is > 0 in >= 90% of assignments for E3; a median
#           calibration ratio (expected / observed test sites) outside 0.7-1.4 is reported as a
#           scale bias; separation or dropped columns in any fold are reported.
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

OUT_DIR <- dir_of("12_validation", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
BLOCK <- VALIDATION$block_km; NF <- VALIDATION$folds; BUF <- VALIDATION$buffer_km
NREP <- 50; FLOOR_FRAC <- 0.01; SEP_THR <- 8
SEED_V <- SEEDS$master + unname(SEEDS$offsets["validation"]) + 1
CONT <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5")
LITH <- c("lith_basalt", "lith_sandstone", "lith_carbonate", "lith_other")
MODELS <- list(E1 = CONT, E2 = c(CONT, LITH), E3 = names(covs), E4 = c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", LITH))

f_log <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv"); stopifnot(file.exists(f_log))
LOG0 <- read.csv(f_log, stringsAsFactors = FALSE)
LOG_COLS <- c("date", "item", "reason", "effect_on_inference", "approved_by")
ROW <- data.frame(
  date = "2026-11-14",
  item = "Clarification: block validation metric and implementation (held-out Poisson log-likelihood gain over a homogeneous intensity; pixel-level fit)",
  reason = "The freeze fixes blocks (60 km), folds (5) and buffer (10 km) but not the metric. Chosen before any validation result: held-out Poisson log-likelihood gain per test site over a homogeneous intensity fitted on the same training region; pixel-level Poisson regression on the 500 m grid (exposure 0.25 km2 per pixel) approximating the Berman-Turner fit; buffer = distance to the nearest test-block rectangle (slightly conservative where blocks overhang the window); predicted intensity floored at 1% of the homogeneous intensity (unfloored also reported); 50 random block-to-fold assignments; E1, E2, E3, E4.",
  effect_on_inference = "Validation addresses first-order transfer only and is not evidence about the second-order deficiency or about interaction. Gate 6 (FAIL) is unaffected.",
  approved_by = "author-delegated (assistant decision; clarification, no frozen item changed)", stringsAsFactors = FALSE)
if (identical(names(LOG0), LOG_COLS)) {
  file.copy(f_log, file.path(OUT_DIR, "protocol_deviation_log_before_56.csv"), overwrite = TRUE)
  if (!(ROW$item %in% LOG0$item)) write.csv(rbind(LOG0, ROW), f_log, row.names = FALSE)
  LOG_MSG <- "Log row appended (backup in outputs)."
} else { write.csv(ROW, file.path(OUT_DIR, "56_log_row_for_manual_append.csv"), row.names = FALSE); LOG_MSG <- "Log columns differ; row saved for manual append." }

# ---- pixel table ----
im1 <- covs[[1]]; nr <- nrow(im1$v); nc <- ncol(im1$v)
stopifnot(all(diff(im1$yrow) > 0), all(diff(im1$xcol) > 0))
gx <- rep(im1$xcol, each = nr); gy <- rep(im1$yrow, times = nc)
A_PIX <- im1$xstep * im1$ystep
Zall <- vapply(covs, function(im) as.vector(im$v), numeric(nr * nc))
ok <- inside.owin(gx, gy, W_km) & stats::complete.cases(Zall)
ix <- round((X_obs_km$x - im1$xcol[1]) / im1$xstep) + 1; iy <- round((X_obs_km$y - im1$yrow[1]) / im1$ystep) + 1
lin <- (ix - 1) * nr + iy
counts <- tabulate(lin, nbins = nr * nc)
fr <- Frame(W_km); xmin <- fr$xrange[1]; ymin <- fr$yrange[1]
bx <- floor((gx - xmin) / BLOCK); by <- floor((gy - ymin) / BLOCK); bid <- bx * 100 + by
blocks <- sort(unique(bid[ok]))
rect_dist <- function(b) {
  xlo <- xmin + floor(b / 100) * BLOCK; ylo <- ymin + (b %% 100) * BLOCK
  dx <- pmax(xlo - gx, 0, gx - (xlo + BLOCK)); dy <- pmax(ylo - gy, 0, gy - (ylo + BLOCK)); sqrt(dx^2 + dy^2)
}
fit_glm <- function(cols, tr) {
  X <- cbind(1, Zall[tr, cols, drop = FALSE]); y <- counts[tr]
  keep <- c(TRUE, colSums(abs(X[, -1, drop = FALSE])) > 0)
  g <- suppressWarnings(stats::glm.fit(X[, keep, drop = FALSE], y, family = poisson(), offset = rep(log(A_PIX), nrow(X)), intercept = FALSE))
  beta <- rep(0, ncol(X)); beta[keep] <- g$coefficients; beta[is.na(beta)] <- 0
  list(beta = beta, dropped = sum(!keep), sep = max(abs(beta[-1])) > SEP_THR)
}
eval_set <- function(beta, cols, te, lam0) {
  lam <- as.vector(exp(cbind(1, Zall[te, cols, drop = FALSE]) %*% beta)); y <- counts[te]
  ll <- function(l) sum(y * log(l * A_PIX) - l * A_PIX)
  list(ll = ll(lam), llf = ll(pmax(lam, FLOOR_FRAC * lam0)), ll0 = ll(rep(lam0, length(y))), n = sum(y), exp = sum(lam) * A_PIX)
}

REPORT <- file.path(OUT_DIR, "56_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("56_block_validation | ", format(Sys.time()), "\n", sep = "")
  cat(LOG_MSG, "\n")
  cat(sprintf("blocks %g km, folds %d, buffer %g km | assignments %d | seed %d | floor %.0f%% of homogeneous intensity\n", BLOCK, NF, BUF, NREP, SEED_V, 100 * FLOOR_FRAC))
  
  cat("\n== 0. Grid and block checks ==\n")
  cat(sprintf("pixels %d x %d | valid pixels %d (area %.1f km2) | sites mapped %d of %d | max site-to-pixel-centre distance %.3f km (limit %.3f)\n",
              nr, nc, sum(ok), sum(ok) * A_PIX, sum(counts[ok]), npoints(X_obs_km), max(sqrt((X_obs_km$x - gx[lin])^2 + (X_obs_km$y - gy[lin])^2)), im1$xstep * sqrt(2) / 2))
  stopifnot(sum(counts[ok]) == npoints(X_obs_km))
  cat(sprintf("blocks with valid pixels: %d | sites per block: median %d, range %d-%d\n", length(blocks),
              median(vapply(blocks, function(b) sum(counts[ok & bid == b]), 1)), min(vapply(blocks, function(b) sum(counts[ok & bid == b]), 1)),
              max(vapply(blocks, function(b) sum(counts[ok & bid == b]), 1))))
  fa <- fit_glm(names(covs), ok); f0 <- fast_fit(X_obs_km)
  cat(sprintf("pixel-level E3 fit vs primary Berman-Turner fit: max |coefficient difference| %.3f (primary coefficients range %.2f to %.2f)\n",
              max(abs(fa$beta - f0$beta[names(c("(Intercept)" = 0, covs))])), min(f0$beta), max(f0$beta)))
  
  cat("\n== 1. Cross-validation ==\n")
  set.seed(SEED_V); ROWS <- data.frame(); t0 <- unname(proc.time()["elapsed"])
  for (rep in seq_len(NREP)) {
    fold_of <- setNames(rep_len(seq_len(NF), length(blocks)), sample(blocks))
    for (f in seq_len(NF)) {
      tb <- as.numeric(names(fold_of)[fold_of == f])
      te <- ok & bid %in% tb
      dmin <- Reduce(pmin, lapply(tb, rect_dist))
      tr <- ok & !(bid %in% tb) & dmin >= BUF
      n_te <- sum(counts[te]); n_tr <- sum(counts[tr])
      if (n_te == 0 || n_tr < 20) next
      lam0 <- n_tr / (sum(tr) * A_PIX)
      for (m in names(MODELS)) {
        fg <- fit_glm(MODELS[[m]], tr); ev <- eval_set(fg$beta, MODELS[[m]], te, lam0)
        ROWS <- rbind(ROWS, data.frame(rep = rep, fold = f, model = m, n_train_pix = sum(tr), n_train_sites = n_tr, n_test_pix = sum(te), n_test_sites = n_te,
                                       gain_site = (ev$ll - ev$ll0) / n_te, gain_site_floor = (ev$llf - ev$ll0) / n_te, calib = ev$exp / n_te, dropped = fg$dropped, sep = fg$sep))
      }
    }
    if (rep %% 10 == 0) cat(sprintf("  assignment %d of %d done [%.0f s]\n", rep, NREP, unname(proc.time()["elapsed"]) - t0))
  }
  write.csv(ROWS, file.path(OUT_DIR, "56_blockcv.csv"), row.names = FALSE)
  
  cat("\n== 2. Pooled results over assignments (floored gain per held-out site; log-likelihood units) ==\n")
  SUMR <- data.frame()
  for (m in names(MODELS)) {
    d <- ROWS[ROWS$model == m, ]
    pg <- vapply(split(d, d$rep), function(x) sum(x$gain_site_floor * x$n_test_sites) / sum(x$n_test_sites), 1)
    pu <- vapply(split(d, d$rep), function(x) sum(x$gain_site * x$n_test_sites) / sum(x$n_test_sites), 1)
    cb <- vapply(split(d, d$rep), function(x) sum(x$calib * x$n_test_sites) / sum(x$n_test_sites), 1)
    SUMR <- rbind(SUMR, data.frame(model = m, df = length(MODELS[[m]]), pooled_gain_median = median(pg), q25 = quantile(pg, .25), q75 = quantile(pg, .75),
                                   share_assign_gain_gt0 = mean(pg > 0), share_folds_gain_gt0 = mean(d$gain_site_floor > 0), unfloored_median = median(pu),
                                   calib_median = median(cb), folds_with_separation = mean(d$sep), folds_with_dropped_cols = mean(d$dropped > 0),
                                   train_sites_median = median(d$n_train_sites), test_sites_median = median(d$n_test_sites), folds = nrow(d)))
  }
  print(SUMR, digits = 3, row.names = FALSE)
  insamp <- vapply(names(MODELS), function(m) { fg <- fit_glm(MODELS[[m]], ok); ev <- eval_set(fg$beta, MODELS[[m]], ok, sum(counts[ok]) / (sum(ok) * A_PIX)); (ev$ll - ev$ll0) / ev$n }, 1)
  cat("\nin-sample gain per site (same models, all data):", paste(sprintf("%s %.3f", names(insamp), insamp), collapse = " | "), "\n")
  e3 <- SUMR[SUMR$model == "E3", ]
  cat(sprintf("\nverdict (pre-set rule) E3: share of assignments with pooled gain > 0 = %.2f (>= 0.90: %s) | median calibration ratio %.2f (0.7-1.4: %s)\n",
              e3$share_assign_gain_gt0, e3$share_assign_gain_gt0 >= 0.90, e3$calib_median, e3$calib_median >= 0.7 && e3$calib_median <= 1.4))
  cat("E3 verdict:", if (e3$share_assign_gain_gt0 >= 0.90) "environmental model TRANSFERS (better than a homogeneous intensity out of block)" else "transfer NOT established at the pre-set level", "\n")
  cat("\n-- first assignment, E3, by fold --\n")
  print(ROWS[ROWS$rep == 1 & ROWS$model == "E3", c("fold", "n_train_sites", "n_test_sites", "n_train_pix", "n_test_pix", "gain_site", "gain_site_floor", "calib", "dropped", "sep")], digits = 3, row.names = FALSE)
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "56_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")