# =====================================================================
# 32_lgcp_generator_decomposition.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : locate the excess pair correlation of rLGCP patterns (script 31):
#           (A) field second moment from the saved Lambda images:
#               sum(L_i L_j) / sum(E_i E_j) and mean(Z_i Z_j) over pixel pairs
#               at r in R_SHOW, vs theory exp(var*exp(-r/scale)) and
#               var*exp(-r/scale);
#           (B) at 4x intensity, pcfinhom variants on the same patterns:
#               V1 true lambda (bw 0.5, Ripley), V2 true lambda (bw 1.0,
#               translate), V3 realised Lambda at the points (expect 1),
#               V4 points redrawn with rpoispp(Lambda), scored with true lambda.
# Outputs : outputs/32_report.txt, 32_field_g.csv, 32_pcf_checks.csv,
#           32_sessionInfo.txt
# Notes   : Synthetic patterns only. NO fit to the archaeological data.
#           Whole run under withCallingHandlers: warnings are tallied.
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

SEED5    <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 19000
HI_MULT  <- 4
N_FIELD  <- 10
N_X4     <- 10
RG       <- seq(0, 15, by = 0.5)
R_SHOW   <- c(1, 2, 3, 5, 8, 12)
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
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
fill_im <- function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
gen_lgcp <- function(v, s, shift, nsim, seed, save) {
  set.seed(seed)
  rLGCP(model = "exp", mu = mu_fill + shift - v / 2, param = list(var = v, scale = s),
        win = W_km, nsim = nsim, saveLambda = save, drop = TRUE)
}
shift_x <- function(A, k) list(a = A[, seq_len(ncol(A) - k), drop = FALSE], b = A[, seq_len(ncol(A) - k) + k, drop = FALSE])
shift_y <- function(A, k) list(a = A[seq_len(nrow(A) - k), , drop = FALSE], b = A[seq_len(nrow(A) - k) + k, , drop = FALSE])
sums_dir <- function(sh, L, E, Z, k) {
  l <- sh(L, k); e <- sh(E, k); z <- sh(Z, k)
  ok <- is.finite(l$a * l$b) & is.finite(e$a * e$b) & is.finite(z$a * z$b)
  c(sum((l$a * l$b)[ok]), sum((e$a * e$b)[ok]), sum((z$a * z$b)[ok]), sum(ok))
}
field_pairs <- function(P, v) tryCatch({
  Lam <- attr(P, "Lambda"); M <- Lam$v; nr <- nrow(M); nc <- ncol(M)
  gx <- rep(Lam$xcol, each = nr); gy <- rep(Lam$yrow, times = nc)
  ins <- matrix(inside.owin(gx, gy, W_km), nrow = nr)
  Mu  <- matrix(lookup.im(mu_fill, gx, gy, naok = TRUE), nrow = nr)
  ok  <- ins & is.finite(M) & is.finite(Mu) & (M > 0)
  L <- ifelse(ok, M, NA_real_); E <- ifelse(ok, exp(Mu), NA_real_); Z <- ifelse(ok, log(M) - (Mu - v / 2), NA_real_)
  vapply(R_SHOW, function(r) {
    kx <- max(1L, as.integer(round(r / Lam$xstep))); ky <- max(1L, as.integer(round(r / Lam$ystep)))
    s <- sums_dir(shift_x, L, E, Z, kx) + sums_dir(shift_y, L, E, Z, ky)
    c(s[1] / s[2], s[3] / s[4])
  }, numeric(2))
}, error = function(e) { cat("  field_pairs FAILED:", conditionMessage(e), "\n"); matrix(NA_real_, 2, length(R_SHOW)) })
pcf_at <- function(Q, lam, bw, corr) tryCatch({
  gh <- do.call(pcfinhom, c(list(Q, lambda = lam, r = RG, correction = corr, bw = bw), PCF_ARGS))
  col_of(gh, c("iso", "trans", "border"))[match(round(R_SHOW, 6), round(gh$r, 6))]
}, error = function(e) rep(NA_real_, length(R_SHOW)))
msd <- function(M) rbind(mean = colMeans(M, na.rm = TRUE), se = apply(M, 2, sd, na.rm = TRUE) / sqrt(nrow(M)))

REPORT <- file.path(OUT_DIR, "32_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  withCallingHandlers({
    cat("32_lgcp_generator_decomposition | ", format(Sys.time()), "\n\n", sep = "")
    mu0 <- Reduce(`+`, Map(function(a, b) a * b, covs, truth$beta)) + truth$b0
    mu_fill <- fill_im(mu0)
    lam4 <- exp(mu_fill + log(HI_MULT))
    FG <- data.frame(); PC <- data.frame()
    
    cat("== A. Field second moment from saved Lambda (x1, ", N_FIELD, " sims per cell) ==\n", sep = "")
    for (nm in names(LG)) {
      v <- unname(LG[[nm]]["var"]); s <- unname(LG[[nm]]["scale"])
      P1 <- try_(sprintf("rLGCP %s x1 saveLambda", nm), gen_lgcp(v, s, 0, N_FIELD, SEED5 + 10, TRUE))
      if (is.null(P1)) next
      FP <- lapply(as_list(P1), field_pairs, v = v)
      G <- t(vapply(FP, function(m) m[1, ], numeric(length(R_SHOW))))
      C <- t(vapply(FP, function(m) m[2, ], numeric(length(R_SHOW))))
      thg <- exp(v * exp(-R_SHOW / s)); thc <- v * exp(-R_SHOW / s)
      mg <- msd(G); mc <- msd(C)
      cat(sprintf("\n-- %s: var %.1f, scale %.0f km --\n", nm, v, s))
      print(data.frame(r_km = R_SHOW, theory_g = thg, field_g = mg["mean", ], se_g = mg["se", ],
                       ratio_g = mg["mean", ] / thg, theory_cov = thc, field_cov = mc["mean", ], se_cov = mc["se", ]),
            digits = 3, row.names = FALSE)
      FG <- rbind(FG, data.frame(cell = nm, r_km = R_SHOW, theory_g = thg, field_g = mg["mean", ], se_g = mg["se", ],
                                 theory_cov = thc, field_cov = mc["mean", ], se_cov = mc["se", ]))
    }
    if (nrow(FG)) write.csv(FG, file.path(OUT_DIR, "32_field_g.csv"), row.names = FALSE)
    
    cat("\n== B. pcfinhom variants at x", HI_MULT, " (", N_X4, " sims per cell) ==\n", sep = "")
    for (nm in names(LG)) {
      v <- unname(LG[[nm]]["var"]); s <- unname(LG[[nm]]["scale"])
      P4 <- try_(sprintf("rLGCP %s x%d saveLambda", nm, HI_MULT), gen_lgcp(v, s, log(HI_MULT), N_X4, SEED5 + 20, TRUE))
      if (is.null(P4)) next
      P4l <- as_list(P4)
      V <- list(V1 = NULL, V2 = NULL, V3 = NULL, V4 = NULL); nn <- numeric(0); nna <- 0
      set.seed(SEED5 + 30)
      for (P in P4l) {
        Lam <- attr(P, "Lambda"); Q <- to_poly(P); nn <- c(nn, npoints(Q))
        lt <- lookup.im(lam4, Q$x, Q$y, naok = TRUE)
        lr <- lookup.im(Lam, Q$x, Q$y, naok = TRUE)
        bad <- is.na(lt) | is.na(lr) | lr <= 0; nna <- nna + sum(bad)
        Qk <- Q[!bad]; lt <- lt[!bad]; lr <- lr[!bad]
        V$V1 <- rbind(V$V1, pcf_at(Qk, lt, 0.5, "Ripley"))
        V$V2 <- rbind(V$V2, pcf_at(Qk, lt, 1.0, "translate"))
        V$V3 <- rbind(V$V3, pcf_at(Qk, lr, 0.5, "Ripley"))
        Pown <- tryCatch(to_poly(rpoispp(Lam)), error = function(e) NULL)
        if (is.null(Pown)) { V$V4 <- rbind(V$V4, rep(NA_real_, length(R_SHOW))); next }
        lo <- lookup.im(lam4, Pown$x, Pown$y, naok = TRUE); Pown <- Pown[!is.na(lo)]; lo <- lo[!is.na(lo)]
        V$V4 <- rbind(V$V4, pcf_at(Pown, lo, 0.5, "Ripley"))
      }
      cat(sprintf("\n-- %s: var %.1f, scale %.0f km | mean n %.1f (expect %d) | points dropped for NA/zero lambda: %d --\n",
                  nm, v, s, mean(nn), HI_MULT * RECOVERY_ENV$expected_n, nna))
      th <- exp(v * exp(-R_SHOW / s))
      labs <- c(V1 = "true lambda, bw 0.5, Ripley", V2 = "true lambda, bw 1.0, translate",
                V3 = "realised Lambda at points (expect 1)", V4 = "own rpoispp(Lambda), true lambda")
      tab <- data.frame(r_km = R_SHOW, theory = th)
      for (k in names(V)) { m <- msd(V[[k]]); tab[[paste0(k, "_mean")]] <- m["mean", ]; tab[[paste0(k, "_se")]] <- m["se", ]
      PC <- rbind(PC, data.frame(cell = nm, variant = labs[[k]], r_km = R_SHOW,
                                 theory = if (k == "V3") 1 else th, mean = m["mean", ], se = m["se", ])) }
      print(tab, digits = 3, row.names = FALSE)
      cat("  variants: V1", labs[["V1"]], "| V2", labs[["V2"]], "| V3", labs[["V3"]], "| V4", labs[["V4"]], "\n")
    }
    if (nrow(PC)) write.csv(PC, file.path(OUT_DIR, "32_pcf_checks.csv"), row.names = FALSE)
    
    cat("\n== Warning tally (whole run) ==\n")
    if (length(WARN$msgs)) print(head(sort(table(substr(WARN$msgs, 1, 160)), decreasing = TRUE), 10)) else cat("none\n")
    cat("\n== END ==\n")
  }, warning = note_w)
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "32_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")