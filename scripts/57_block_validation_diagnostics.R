# =====================================================================
# 57_block_validation_diagnostics.R
# Stage   : 12_validation
# Folder  : 12_validation/scripts/   (NEW)
# Purpose : EXPLORATORY diagnostics of the negative out-of-block gain in script 56:
#           (A) from 56_blockcv.csv: pooled gain excluding folds flagged for separation; fold gain
#               vs test-site count tercile and vs |log calibration|;
#           (B) rank skill (mean percentile of held-out sites among all test pixels by predicted
#               intensity; 0.5 = chance) and floored log-likelihood gain for E1 and E3 at block sizes
#               30, 60 and 100 km (5 folds, buffer 10 km, 50 assignments each);
#           (C) sign stability of the E3 coefficients over the 250 training fits at 60 km.
# Outputs : outputs/57_report.txt, 57_blockcv_exploratory.csv, 57_coef_stability.csv, 57_sessionInfo.txt
# Notes   : Exploratory; changes no reported result and selects no model. Reading rules (set before
#           running): RANK TRANSFER at 60 km if E3 pooled percentile > 0.5 in >= 90% of assignments;
#           SEPARATION FOLDS explain the negative gain only if the E3 pooled gain excluding them is > 0
#           in >= 90% of assignments. Same pixel grid and fit as script 56.
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
f_cv <- file.path(OUT_DIR, "56_blockcv.csv"); stopifnot(file.exists(f_cv))
CV56 <- read.csv(f_cv, stringsAsFactors = FALSE)
NF <- VALIDATION$folds; BUF <- VALIDATION$buffer_km; NREP <- 50; FLOOR_FRAC <- 0.01; SEP_THR <- 8
BLOCKS <- c(30, 60, 100)
SEED_V <- SEEDS$master + unname(SEEDS$offsets["validation"]) + 2
CONT <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5")
MODELS <- list(E1 = CONT, E3 = names(covs))

im1 <- covs[[1]]; nr <- nrow(im1$v); nc <- ncol(im1$v)
gx <- rep(im1$xcol, each = nr); gy <- rep(im1$yrow, times = nc); A_PIX <- im1$xstep * im1$ystep
Zall <- vapply(covs, function(im) as.vector(im$v), numeric(nr * nc))
ok <- inside.owin(gx, gy, W_km) & stats::complete.cases(Zall)
ix <- round((X_obs_km$x - im1$xcol[1]) / im1$xstep) + 1; iy <- round((X_obs_km$y - im1$yrow[1]) / im1$ystep) + 1
counts <- tabulate((ix - 1) * nr + iy, nbins = nr * nc)
fr <- Frame(W_km); xmin <- fr$xrange[1]; ymin <- fr$yrange[1]
bid_for <- function(blk) floor((gx - xmin) / blk) * 1000 + floor((gy - ymin) / blk)
rect_dist <- function(b, blk) {
  xlo <- xmin + floor(b / 1000) * blk; ylo <- ymin + (b %% 1000) * blk
  dx <- pmax(xlo - gx, 0, gx - (xlo + blk)); dy <- pmax(ylo - gy, 0, gy - (ylo + blk)); sqrt(dx^2 + dy^2)
}
fit_glm <- function(cols, tr) {
  X <- cbind(1, Zall[tr, cols, drop = FALSE]); y <- counts[tr]
  keep <- c(TRUE, colSums(abs(X[, -1, drop = FALSE])) > 0)
  g <- suppressWarnings(stats::glm.fit(X[, keep, drop = FALSE], y, family = poisson(), offset = rep(log(A_PIX), nrow(X)), intercept = FALSE))
  beta <- rep(0, ncol(X)); beta[keep] <- g$coefficients; beta[is.na(beta)] <- 0
  list(beta = beta, sep = max(abs(beta[-1])) > SEP_THR)
}
eval_set <- function(beta, cols, te, lam0) {
  lam <- as.vector(exp(cbind(1, Zall[te, cols, drop = FALSE]) %*% beta)); y <- counts[te]
  ll <- function(l) sum(y * log(l * A_PIX) - l * A_PIX)
  r <- rank(lam, ties.method = "average")
  list(gain = (ll(pmax(lam, FLOOR_FRAC * lam0)) - ll(rep(lam0, length(y)))) / sum(y), pct = sum(y * (r - 0.5) / length(lam)) / sum(y), n = sum(y))
}

REPORT <- file.path(OUT_DIR, "57_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("57_block_validation_diagnostics (EXPLORATORY) | ", format(Sys.time()), "\n", sep = "")
  
  cat("\n== A. From script 56: separation folds and fold-level structure ==\n")
  for (m in c("E1", "E2", "E3", "E4")) {
    d <- CV56[CV56$model == m, ]
    pool <- function(x) if (nrow(x)) sum(x$gain_site_floor * x$n_test_sites) / sum(x$n_test_sites) else NA_real_
    all_g <- vapply(split(d, d$rep), pool, 1); ex_g <- vapply(split(d, d$rep), function(x) pool(x[!x$sep, ]), 1)
    cat(sprintf("  %s: pooled gain all folds median %.3f | excluding separation folds median %.3f (share of assignments > 0: %.2f) | folds with separation %.3f\n",
                m, median(all_g, na.rm = TRUE), median(ex_g, na.rm = TRUE), mean(ex_g > 0, na.rm = TRUE), mean(d$sep)))
  }
  d3 <- CV56[CV56$model == "E3", ]
  tert <- cut(d3$n_test_sites, quantile(d3$n_test_sites, c(0, 1/3, 2/3, 1)), include.lowest = TRUE)
  cat("  E3 fold gain by test-site tercile:\n"); print(aggregate(gain_site_floor ~ tert, data = d3, FUN = function(x) c(mean = mean(x), median = median(x), n = length(x))), digits = 3)
  cat(sprintf("  E3 Spearman(fold gain, |log calibration ratio|) = %.3f | Spearman(fold gain, test sites) = %.3f\n",
              cor(d3$gain_site_floor, abs(log(d3$calib)), method = "spearman"), cor(d3$gain_site_floor, d3$n_test_sites, method = "spearman")))
  
  cat("\n== B/C. Rank skill and gain by block size (E1, E3; 5 folds; buffer ", BUF, " km; ", NREP, " assignments) ==\n", sep = "")
  set.seed(SEED_V); ROWS <- data.frame(); BETAS <- NULL
  for (blk in BLOCKS) {
    bid <- bid_for(blk); blocks <- sort(unique(bid[ok])); t0 <- unname(proc.time()["elapsed"])
    cat(sprintf("  block %g km: %d blocks with valid pixels\n", blk, length(blocks)))
    for (rep in seq_len(NREP)) {
      fold_of <- setNames(rep_len(seq_len(NF), length(blocks)), sample(blocks))
      for (f in seq_len(NF)) {
        tb <- as.numeric(names(fold_of)[fold_of == f]); te <- ok & bid %in% tb
        tr <- ok & !(bid %in% tb) & Reduce(pmin, lapply(tb, rect_dist, blk = blk)) >= BUF
        if (sum(counts[te]) == 0 || sum(counts[tr]) < 20) next
        lam0 <- sum(counts[tr]) / (sum(tr) * A_PIX)
        for (m in names(MODELS)) {
          fg <- fit_glm(MODELS[[m]], tr); ev <- eval_set(fg$beta, MODELS[[m]], te, lam0)
          ROWS <- rbind(ROWS, data.frame(block_km = blk, rep = rep, fold = f, model = m, n_test = ev$n, gain = ev$gain, pct = ev$pct, sep = fg$sep))
          if (blk == 60 && m == "E3") BETAS <- rbind(BETAS, fg$beta[-1])
        }
      }
    }
    cat(sprintf("    done in %.0f s\n", unname(proc.time()["elapsed"]) - t0))
  }
  write.csv(ROWS, file.path(OUT_DIR, "57_blockcv_exploratory.csv"), row.names = FALSE)
  pool2 <- function(d, v) sum(d[[v]] * d$n_test) / sum(d$n_test)
  SM <- do.call(rbind, lapply(split(ROWS, list(ROWS$block_km, ROWS$model)), function(d) {
    pg <- vapply(split(d, d$rep), pool2, 1, v = "gain"); pp <- vapply(split(d, d$rep), pool2, 1, v = "pct")
    data.frame(block_km = d$block_km[1], model = d$model[1], folds = nrow(d), gain_median = median(pg), share_gain_gt0 = mean(pg > 0),
               pct_median = median(pp), pct_q25 = quantile(pp, .25), pct_q75 = quantile(pp, .75), share_pct_gt_half = mean(pp > 0.5),
               fold_share_pct_gt_half = mean(d$pct > 0.5), sep_share = mean(d$sep))
  }))
  print(SM, digits = 3, row.names = FALSE)
  
  cat("\n-- C. E3 coefficient sign stability across training fits at 60 km (", nrow(BETAS), " fits) --\n", sep = "")
  b_full <- fit_glm(MODELS$E3, ok)$beta[-1]; nm <- MODELS$E3
  STAB <- data.frame(term = nm, full_data = b_full, same_sign_share = vapply(seq_along(nm), function(j) mean(sign(BETAS[, j]) == sign(b_full[j])), 1),
                     fold_q05 = apply(BETAS, 2, quantile, 0.05), fold_median = apply(BETAS, 2, median), fold_q95 = apply(BETAS, 2, quantile, 0.95))
  print(STAB, digits = 3, row.names = FALSE)
  write.csv(STAB, file.path(OUT_DIR, "57_coef_stability.csv"), row.names = FALSE)
  
  cat("\n== Verdict (pre-set rules) ==\n")
  e3 <- SM[SM$block_km == 60 & SM$model == "E3", ]
  cat(sprintf("E3 at 60 km: share of assignments with pooled percentile > 0.5 = %.2f -> %s\n", e3$share_pct_gt_half,
              if (e3$share_pct_gt_half >= 0.90) "RANK TRANSFER (order is predicted better than chance, level is not)" else "rank transfer NOT established at the pre-set level"))
  ex3 <- vapply(split(CV56[CV56$model == "E3" & !CV56$sep, ], CV56$rep[CV56$model == "E3" & !CV56$sep]), function(x) sum(x$gain_site_floor * x$n_test_sites) / sum(x$n_test_sites), 1)
  cat(sprintf("E3 gain excluding separation folds: share of assignments > 0 = %.2f -> %s\n", mean(ex3 > 0),
              if (mean(ex3 > 0) >= 0.90) "separation folds explain the negative gain" else "separation folds do NOT explain the negative gain"))
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "57_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")