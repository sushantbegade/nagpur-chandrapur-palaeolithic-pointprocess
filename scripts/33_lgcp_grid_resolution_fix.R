# =====================================================================
# 33_lgcp_grid_resolution_fix.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : test whether the excess pair correlation of rLGCP patterns
#           (scripts 31-32) is caused by rLGCP resampling mu onto a coarse
#           default grid (~1.7 km). Variants: default (control), eps = 0.5,
#           dimyx = dim(mu), 500 m mask window. For each: Lambda pixel step,
#           mean n at 4x, pcfinhom with TRUE lambda vs exp(var*exp(-r/scale)),
#           clik2 sigma2 / alpha at 4x.
# Outputs : outputs/33_report.txt, 33_lgcp_grid.csv, 33_lgcp_grid_summary.csv,
#           33_sessionInfo.txt
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           Passing eps / dimyx through rLGCP '...' is UNVERIFIED; failures
#           print and continue. Acceptance rule (set before running): Lambda
#           step <= 0.6 km AND V1 ratio within 0.85-1.15 at r = 2,3,5,8,12 km
#           for BOTH cells. Whole run under withCallingHandlers (warnings tallied).
#           Restart R before running. Main session only (no workers).
# =====================================================================
rm(list = ls()); while (sink.number() > 0) sink()

CONFIG_PATH <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"  # <- only edit-line
stopifnot(file.exists(CONFIG_PATH)); source(CONFIG_PATH)

CACHE_FILE <- "recovery_cache_simpl.rds"
stopifnot(file.exists(file.path(REPO_DIR, "data_derived", CACHE_FILE)))
Sys.setenv(RECOVERY_CACHE = CACHE_FILE)
source(file.path(REPO_DIR, "R", "recovery_worker.R"))
source(file.path(REPO_DIR, "R", "recovery_functions.R"))

OUT_DIR <- dir_of("10_simulation_recovery", "outputs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEED6    <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 21000
HI_MULT  <- 4
N_X4     <- 10
N_FIT    <- 6
RG       <- seq(0, 15, by = 0.5)
R_SHOW   <- c(1, 2, 3, 5, 8, 12)
R_RULE   <- c(2, 3, 5, 8, 12)
STEP_MAX <- 0.6
RATIO_LO <- 0.85
RATIO_HI <- 1.15
VARIANTS <- c("default", "eps05", "dimyx", "maskwin")
PCF_ARGS <- list(kernel = "epanechnikov", zerocor = "convolution", divisor = "r", renormalise = FALSE)
LG       <- list(L1 = c(var = 0.5, scale = 3), L4 = c(var = 1.5, scale = 8))

# ---- helpers ----
WARN <- new.env(); WARN$msgs <- character(0)
note_w  <- function(w) { WARN$msgs <- c(WARN$msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
try_ <- function(label, expr) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch(expr, error = function(e) { cat(sprintf("[%s] FAILED: %s\n", label, conditionMessage(e))); NULL })
  cat(sprintf("[%s] %.1f s\n", label, unname(proc.time()["elapsed"]) - t0))
  out
}
quiet   <- function(expr) tryCatch(expr, error = function(e) NULL)
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
fill_im <- function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
gen_var <- function(variant, v, s, shift, nsim, seed) {
  set.seed(seed)
  mu <- mu_fill + shift - v / 2
  switch(variant,
         default = rLGCP(model = "exp", mu = mu, param = list(var = v, scale = s), win = W_km,
                         nsim = nsim, saveLambda = TRUE, drop = TRUE),
         eps05   = rLGCP(model = "exp", mu = mu, param = list(var = v, scale = s), win = W_km,
                         nsim = nsim, saveLambda = TRUE, drop = TRUE, eps = 0.5),
         dimyx   = rLGCP(model = "exp", mu = mu, param = list(var = v, scale = s), win = W_km,
                         nsim = nsim, saveLambda = TRUE, drop = TRUE, dimyx = dim(mu_fill)),
         maskwin = rLGCP(model = "exp", mu = mu, param = list(var = v, scale = s), win = WMASK,
                         nsim = nsim, saveLambda = TRUE, drop = TRUE))
}
pcf_at <- function(Q, lam) tryCatch({
  gh <- do.call(pcfinhom, c(list(Q, lambda = lam, r = RG, correction = "Ripley", bw = 0.5), PCF_ARGS))
  col_of(gh, c("iso", "trans", "border"))[match(round(R_SHOW, 6), round(gh$r, 6))]
}, error = function(e) rep(NA_real_, length(R_SHOW)))
fit2 <- function(X)
  quiet(kppm(X, trend = formula_E3, clusters = "LGCP", model = "exponential", covariates = covs, method = "clik2"))
pars_of <- function(f) if (is.null(f)) c(NA_real_, NA_real_) else unname(c(f$par["sigma2"], f$par["alpha"]))

REPORT <- file.path(OUT_DIR, "33_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("33_lgcp_grid_resolution_fix | ", format(Sys.time()), "\n", sep = "")
    mu0 <- Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta)) + truth$b0
    mu_fill <- fill_im(mu0)
    lam4 <- exp(mu_fill + log(HI_MULT))
    cat(sprintf("mu_fill grid: %d x %d pixels, step x %.3f km, step y %.3f km\n",
                nrow(mu_fill$v), ncol(mu_fill$v), mu_fill$xstep, mu_fill$ystep))
    WMASK <- as.mask(W_km, eps = 0.5)
    cat(sprintf("Mask window: %d x %d pixels, step %.3f km\n", WMASK$dim[1], WMASK$dim[2], WMASK$xstep))
    RES <- data.frame(); SUM <- data.frame()
    
    for (nm in names(LG)) {
      v <- unname(LG[[nm]]["var"]); s <- unname(LG[[nm]]["scale"])
      th <- exp(v * exp(-R_SHOW / s))
      for (vr in VARIANTS) {
        cat(sprintf("\n-- %s (var %.1f, scale %.0f km), variant %s --\n", nm, v, s, vr))
        P <- try_(sprintf("rLGCP %s %s x%d, %d sims", nm, vr, HI_MULT, N_X4),
                  gen_var(vr, v, s, log(HI_MULT), N_X4, SEED6 + 10))
        if (is.null(P)) { SUM <- rbind(SUM, data.frame(cell = nm, variant = vr, ok = FALSE, step_km = NA_real_,
                                                       mean_n = NA_real_, s2_med = NA_real_, al_med = NA_real_, min_ratio = NA_real_,
                                                       max_ratio = NA_real_, pass = FALSE)); next }
        Pl <- as_list(P)
        Lam1 <- attr(Pl[[1]], "Lambda")
        step <- if (is.null(Lam1)) NA_real_ else Lam1$xstep
        Pq <- lapply(Pl, to_poly); nn <- vapply(Pq, npoints, numeric(1))
        cat(sprintf("  Lambda pixel step %.3f km | mean n %.1f (expect %d)\n", step, mean(nn), HI_MULT * RECOVERY_ENV$expected_n))
        G <- t(vapply(Pq, function(Q) {
          lt <- lookup.im(lam4, Q$x, Q$y, naok = TRUE); ok <- !is.na(lt) & lt > 0
          pcf_at(Q[ok], lt[ok]) }, numeric(length(R_SHOW))))
        m <- colMeans(G, na.rm = TRUE); se <- apply(G, 2, sd, na.rm = TRUE) / sqrt(nrow(G))
        print(data.frame(r_km = R_SHOW, theory = th, V1_mean = m, se = se, ratio = m / th), digits = 3, row.names = FALSE)
        fits <- lapply(Pq[seq_len(min(N_FIT, length(Pq)))], fit2)
        pp <- vapply(fits, pars_of, numeric(2))
        cat(sprintf("  clik2 x%d: ok %d/%d | sigma2 median %.3f (truth %.2f) | alpha median %.2f (truth %.0f)\n",
                    HI_MULT, sum(!is.na(pp[1, ])), ncol(pp), median(pp[1, ], na.rm = TRUE), v,
                    median(pp[2, ], na.rm = TRUE), s))
        ratio <- (m / th)[R_SHOW %in% R_RULE]
        pass <- is.finite(step) && step <= STEP_MAX && all(ratio >= RATIO_LO & ratio <= RATIO_HI)
        RES <- rbind(RES, data.frame(cell = nm, variant = vr, r_km = R_SHOW, theory = th, mean = m, se = se))
        SUM <- rbind(SUM, data.frame(cell = nm, variant = vr, ok = TRUE, step_km = step, mean_n = mean(nn),
                                     s2_med = median(pp[1, ], na.rm = TRUE), al_med = median(pp[2, ], na.rm = TRUE),
                                     min_ratio = min(ratio), max_ratio = max(ratio), pass = pass))
        cat(sprintf("  rule: step <= %.1f km and ratio in [%.2f, %.2f] at r >= 2 -> %s (ratio range %.2f-%.2f)\n",
                    STEP_MAX, RATIO_LO, RATIO_HI, if (pass) "PASS" else "FAIL", min(ratio), max(ratio)))
      }
    }
    if (nrow(RES)) write.csv(RES, file.path(OUT_DIR, "33_lgcp_grid.csv"), row.names = FALSE)
    if (nrow(SUM)) write.csv(SUM, file.path(OUT_DIR, "33_lgcp_grid_summary.csv"), row.names = FALSE)
    
    cat("\n== SUMMARY ==\n"); print(SUM, digits = 3, row.names = FALSE)
    both <- vapply(VARIANTS, function(vr) all(SUM$pass[SUM$variant == vr]) && sum(SUM$variant == vr) == length(LG), logical(1))
    cat("\nVariants passing for BOTH cells:", if (any(both)) paste(VARIANTS[both], collapse = ", ") else "none", "\n")
    
    cat("\n== Warning tally (whole run) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "33_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")