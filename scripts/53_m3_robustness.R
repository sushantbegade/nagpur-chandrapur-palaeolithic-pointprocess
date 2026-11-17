# =====================================================================
# 53_m3_robustness.R
# Stage   : 11_generative_gof
# Folder  : 11_generative_gof/scripts/   (NEW; writes R/robust_m3.R)
# Purpose : M3 (E + LGCP: kppm exponential clik2; conditional simulation rLGCP eps 0.5 in chunks
#           of 100 + rpoint(n | Lambda); nsim 999; joint ERL L_inhom + g_inhom, nstep 2) under
#           the robustness settings R01 redeposited-excluded (195), R06 eroded window 5 km (169),
#           R08 E2 (10 df), R09 E4 (8 df). 4 workers, one setting per worker. Generic curve
#           code identical to script 49. Worker-vs-main identity check first (R08, nsim 19).
# Outputs : outputs/53_report.txt, 53_m3_robustness_summary.csv, 53_sessionInfo.txt;
#           R/robust_m3.R
# Notes   : RUN ONCE (stops if the summary exists). Reading rules (set before running):
#           p >= 0.06 -> M3 adequate (outcome retained); p in [0.04, 0.06] -> precision rule
#           (nsim 10000 rerun, script 53b); p < 0.04 -> M3 deficient in that setting (then M4
#           would be descriptive only; no interaction claim, Gate 6). sigma2_hat / alpha_hat are
#           estimator outputs. Restart R before running. Workers load the cache and generic files
#           only (no terra).
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
f_res <- file.path(OUT_DIR, "53_m3_robustness_summary.csv")
if (file.exists(f_res)) stop("53_m3_robustness_summary.csv exists: this analysis runs once. Stop; do not rerun.")
GENERIC_FILE <- file.path(REPO_DIR, "R", "robust_generic.R"); stopifnot(file.exists(GENERIC_FILE))
M3_FILE <- file.path(REPO_DIR, "R", "robust_m3.R")
NSIM <- TESTS$nsim_empirical; N_WORK <- 4
SEED_M <- SEEDS$master + unname(SEEDS$offsets["empirical"]) + 100000
RUN <- c("R01_redep", "R06_eroded5", "R08_E2", "R09_E4")
CONT <- c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", "z_ldist_riv3", "z_ldist_riv5")
LITH <- c("lith_basalt", "lith_sandstone", "lith_carbonate", "lith_other")

# ---- open inputs BEFORE sink ----
f_man <- file.path(REPO_DIR, "docs", "freeze_v1_manifest.csv"); stopifnot(file.exists(f_man))
MAN <- read.csv(f_man, stringsAsFactors = FALSE)
f_sites <- MAN$path[basename(MAN$path) == "sites_NA.csv"]; stopifnot(length(f_sites) == 1, file.exists(f_sites))
SITES_G <- read.csv(f_sites, stringsAsFactors = FALSE)
f_m3 <- file.path(OUT_DIR, "46_m3_result.csv"); f_r49 <- file.path(OUT_DIR, "49_robustness_summary.csv")
stopifnot(file.exists(f_m3), file.exists(f_r49))
R46 <- read.csv(f_m3, stringsAsFactors = FALSE); R49 <- read.csv(f_r49, stringsAsFactors = FALSE)
source(GENERIC_FILE)

# ---- M3 functions (dumped for workers) ----
as_list <- function(x) if (is.list(x) && !is.ppp(x)) x else list(x)
fill_im <- function(m) im(ifelse(is.na(m$v), median(m$v, na.rm = TRUE), m$v), xcol = m$xcol, yrow = m$yrow)
make_setting_m3 <- function(nm) {
  if (nm == "R01_redep") return(mk(nm, X_obs_km[!SITES_G$sens_redeposited], W_km, covs, DUMMY_FIXED))
  if (nm == "R06_eroded5") { Wer <- erosion(W_km, 5); Xe <- X_obs_km[Wer]; return(mk(nm, Xe, Wer, covs, default.dummy(Xe, eps = QUAD_KM))) }
  if (nm == "R08_E2") return(mk(nm, X_obs_km, W_km, covs[c(CONT, LITH)], DUMMY_FIXED))
  if (nm == "R09_E4") return(mk(nm, X_obs_km, W_km, covs[c("z_elev_p1", "z_elev_p2", "z_lslope", "z_tpi2km", LITH)], DUMMY_FIXED))
  stop("unknown setting ", nm)
}
m3_run <- function(S, nsim, seed, chunk = 100) {
  set.seed(seed)
  fm <- as.formula(paste("~", paste(names(S$cov), collapse = " + ")))
  fit <- kppm(S$X, trend = fm, clusters = "LGCP", model = "exponential", covariates = S$cov, method = "clik2")
  s2 <- unname(fit$par["sigma2"]); al <- unname(fit$par["alpha"]); b <- coef(fit)
  eta <- Reduce(`+`, Map(function(a, nmc) a * unname(b[nmc]), S$cov, names(S$cov))) + unname(b["(Intercept)"])
  mu <- fill_im(eta); n <- npoints(S$X)
  ft <- rob_fit(S$X, S); ob <- rob_curves(S$X, ft, S)
  sL <- matrix(NA_real_, length(S$rL), nsim); sG <- matrix(NA_real_, length(S$rG), nsim); nn <- rep(NA_real_, nsim)
  k0 <- 0
  while (k0 < nsim) {
    m <- min(chunk, nsim - k0)
    Sl <- as_list(rLGCP(model = "exp", mu = mu - s2 / 2, param = list(var = s2, scale = al), win = S$W, eps = 0.5,
                        saveLambda = TRUE, nsim = m, drop = TRUE))
    for (j in seq_len(m)) {
      res <- tryCatch({
        Xs <- rob_poly(rpoint(n, f = attr(Sl[[j]], "Lambda"), win = S$W, forcewin = TRUE), S)
        Xc <- rob_clean(Xs, S); fs <- rob_fit(Xc, S); list(cs = rob_curves(Xc, fs, S), n = npoints(Xc))
      }, error = function(e) NULL)
      if (!is.null(res)) { sL[, k0 + j] <- res$cs$L; sG[, k0 + j] <- res$cs$g; nn[k0 + j] <- res$n }
    }
    k0 <- k0 + m
  }
  okc <- colSums(is.na(sL)) == 0 & colSums(is.na(sG)) == 0
  csL <- GET::create_curve_set(list(r = S$rL, obs = ob$L, sim_m = sL[, okc, drop = FALSE]))
  csG <- GET::create_curve_set(list(r = S$rG, obs = ob$g, sim_m = sG[, okc, drop = FALSE]))
  list(p = erl_p(csL, csG, 2), usable = sum(okc), sigma2 = s2, alpha = al, n_mean = mean(nn[okc]), n = n)
}
m3_task <- function(job) {
  t0 <- unname(proc.time()["elapsed"])
  S <- make_setting_m3(job$setting)
  r <- tryCatch(m3_run(S, job$nsim, job$seed), error = function(e) list(err = conditionMessage(e)))
  c(list(setting = job$setting, sec = unname(proc.time()["elapsed"]) - t0), r)
}
m3_row <- function(x) {
  g <- function(f) if (is.null(x[[f]])) NA_real_ else as.numeric(x[[f]])
  data.frame(setting = x$setting, n = g("n"), nsim = g("usable"), p = g("p"), sigma2 = g("sigma2"), alpha = g("alpha"),
             n_sim_mean = g("n_mean"), sec = g("sec"), err = if (is.null(x$err)) NA_character_ else as.character(x$err), stringsAsFactors = FALSE)
}
init_worker5 <- function(repo, cache, generic, m3file, sites_path) {
  assign("REPO_DIR", repo, envir = globalenv()); Sys.setenv(RECOVERY_CACHE = cache)
  source(file.path(repo, "R", "recovery_worker.R"), local = globalenv())
  source(file.path(repo, "R", "recovery_functions.R"), local = globalenv())
  assign("SITES_G", read.csv(sites_path, stringsAsFactors = FALSE), envir = globalenv())
  source(generic, local = globalenv()); source(m3file, local = globalenv())
  TRUE
}
dump(c("CONT", "LITH", "as_list", "fill_im", "make_setting_m3", "m3_run", "m3_task", "m3_row"), file = M3_FILE)

REPORT <- file.path(OUT_DIR, "53_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("53_m3_robustness | ", format(Sys.time()), "\n", sep = "")
  cat(sprintf("settings %s | nsim %d | workers %d | seed base %d\n", paste(RUN, collapse = ", "), NSIM, N_WORK, SEED_M))
  cat(sprintf("primary M3 (script 46): p %.4f, sigma2_hat %.4f, alpha_hat %.3f km\n", R46$p_erl_nstep2, R46$sigma2_hat, R46$alpha_hat))
  cl <- makeCluster(N_WORK)
  ini <- clusterCall(cl, init_worker5, REPO_DIR, CACHE_FILE, GENERIC_FILE, M3_FILE, f_sites)
  cat("workers initialised:", all(unlist(ini)), "\n")
  ex <- clusterCall(cl, function() c(exists("m3_task"), exists("SITES_G"), "terra" %in% loadedNamespaces()))
  cat("worker check (m3_task, SITES_G present; terra loaded):", paste(ex[[1]], collapse = " "), "\n")
  
  cat("\n== Worker-vs-main identity check (R08, nsim 19, same seed) ==\n")
  m <- m3_run(make_setting_m3("R08_E2"), 19, SEED_M + 999)
  w <- parLapply(cl[1], list(list(setting = "R08_E2", nsim = 19, seed = SEED_M + 999)), m3_task)[[1]]
  cat(sprintf("main: p %.4f sigma2 %.6f alpha %.5f | worker: p %.4f sigma2 %.6f alpha %.5f | worker time %.0f s\n",
              m$p, m$sigma2, m$alpha, w$p, w$sigma2, w$alpha, w$sec))
  if (!(isTRUE(all.equal(m$p, w$p)) && isTRUE(all.equal(m$sigma2, w$sigma2)) && isTRUE(all.equal(m$alpha, w$alpha)))) stop("worker path differs from main path; stop and paste output")
  
  cat("\n== M3 runs (one setting per worker) ==\n")
  jobs <- lapply(seq_along(RUN), function(k) list(setting = RUN[k], nsim = NSIM, seed = SEED_M + k))
  t0 <- unname(proc.time()["elapsed"])
  res <- parLapplyLB(cl, jobs, m3_task)
  cat(sprintf("wall %.0f s\n", unname(proc.time()["elapsed"]) - t0))
  OUT <- do.call(rbind, lapply(res, m3_row))
  OUT$mc_se <- sqrt(OUT$p * (1 - OUT$p) / OUT$nsim)
  OUT$flag <- ifelse(is.na(OUT$p), "FAILED", ifelse(OUT$p >= 0.04 & OUT$p <= 0.06, "PRECISION", ifelse(OUT$p < 0.04, "M3 DEFICIENT", "M3 ADEQUATE")))
  OUT$m1_p_script49 <- R49$p[match(OUT$setting, R49$setting)]
  OUT$date <- format(Sys.time())
  write.csv(OUT, f_res, row.names = FALSE)
  cat("\n== Results ==\n")
  print(OUT[, c("setting", "n", "nsim", "p", "mc_se", "flag", "sigma2", "alpha", "m1_p_script49", "sec", "err")], digits = 3, row.names = FALSE)
  cat("\nprecision-rule settings (script 53b):", paste(OUT$setting[OUT$flag == "PRECISION"], collapse = " "), "\n")
  cat("\n== END ==\n")
}, finally = { while (sink.number() > 0) sink(); if (exists("cl")) try(stopCluster(cl), silent = TRUE) })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "53_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")