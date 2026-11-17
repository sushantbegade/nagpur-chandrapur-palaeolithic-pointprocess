# =====================================================================
# 25_signature_probe.R
# Stage   : 10_simulation_recovery
# Folder  : 10_simulation_recovery/scripts/   (NEW file, overwrites nothing)
# Purpose : Probe installed signatures and behaviour of kppm / vcov /
#           simulate (LGCP) and ppm / vcov / simulate with new.coef
#           (Poisson, Geyer) BEFORE writing the M3/M4 test functions.
#           Synthetic patterns only. NO fit to the archaeological data.
# Outputs : outputs/25_report.txt, outputs/25_sessionInfo.txt
# Notes   : rmhcontrol(nrep = 1e4) is used for probe simulations ONLY
#           (timing / behaviour, not valid for production).
#           All calls go through try_(); failures print and continue.
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
SEED0      <- SEEDS$master + unname(SEEDS$offsets["recovery"]) + 5000
NREP_PROBE <- 1e4

# ---- open inputs BEFORE sink ----
f_m <- file.path(OUT_DIR, "geyer_multipliers_simpl.csv"); stopifnot(file.exists(f_m))
gm <- read.csv(f_m, stringsAsFactors = FALSE)
MULT_G3 <- gm$mult[gm$cell == "G3"]; stopifnot(length(MULT_G3) == 1)

# ---- helpers (no assignments inside) ----
try_ <- function(label, expr) {
  t0 <- unname(proc.time()["elapsed"])
  out <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      cat(sprintf("  [%s] warning: %s\n", label, conditionMessage(w))); invokeRestart("muffleWarning") }),
    error = function(e) { cat(sprintf("[%s] FAILED: %s\n", label, conditionMessage(e))); NULL })
  cat(sprintf("[%s] %.1f s\n", label, unname(proc.time()["elapsed"]) - t0))
  out
}
show_formals <- function(generic, cls) {
  m <- utils::getS3method(generic, cls, optional = TRUE)
  if (is.null(m)) cat(sprintf("%s.%s: NOT FOUND\n", generic, cls))
  else cat(sprintf("%s.%s formals: %s\n", generic, cls, paste(names(formals(m)), collapse = ", ")))
}
filter_args <- function(generic, cls, args) {
  m <- utils::getS3method(generic, cls, optional = TRUE)
  if (is.null(m)) return(args)
  args[names(args) %in% names(formals(m))]
}
describe_sims <- function(sims, label) {
  L <- if (is.ppp(sims)) list(sims) else as.list(sims)
  for (k in seq_along(L)) {
    P <- L[[k]]
    cat(sprintf("  %s sim %d: class %s, window %s, n = %d, inside W_km = %d\n", label, k, class(P)[1],
                Window(P)$type, npoints(P), sum(inside.owin(P$x, P$y, W_km))))
  }
  invisible(NULL)
}
describe_vcov <- function(V, b, label) {
  cat(sprintf("  %s: dim %s, length(coef) %d, anyNA %s, names match coef %s\n", label,
              paste(dim(V), collapse = "x"), length(b), anyNA(V), identical(rownames(V), names(b))))
  cat("  sqrt(diag) range:", format(range(sqrt(diag(V)), na.rm = TRUE), digits = 4), "\n")
  invisible(NULL)
}
draw_coef <- function(b, V, seed) { set.seed(seed); MASS::mvrnorm(1, mu = b, Sigma = V) }
sim_geyer <- function(mult, gamma, r, sat, n_start, nsim, seed, nrep = 5e5) {
  set.seed(seed)
  rmh(rmhmodel(cif = "geyer", par = list(beta = 1, gamma = gamma, r = r, sat = sat),
               w = W_km, trend = truth$lamW * mult),
      start = rmhstart(n.start = n_start), control = rmhcontrol(nrep = nrep),
      nsim = nsim, drop = TRUE, verbose = FALSE)
}

REPORT <- file.path(OUT_DIR, "25_report.txt")
tryCatch({
  sink(REPORT, split = TRUE)
  cat("25_signature_probe | ", format(Sys.time()), "\n", sep = "")
  cat("spatstat.model", as.character(packageVersion("spatstat.model")),
      "| spatstat.random", as.character(packageVersion("spatstat.random")),
      "| MASS available:", requireNamespace("MASS", quietly = TRUE), "\n\n")
  
  cat("== 0. Formals (S3 methods) ==\n")
  for (g in list(c("kppm", "ppp"), c("kppm", "formula"), c("vcov", "kppm"), c("simulate", "kppm"),
                 c("coef", "kppm"), c("vcov", "ppm"), c("simulate", "ppm"), c("coef", "ppm")))
    show_formals(g[1], g[2])
  cat("kppm generic formals:", paste(names(formals(kppm)), collapse = ", "), "\n")
  cat("new.coef in simulate.ppm:", "new.coef" %in% names(formals(utils::getS3method("simulate", "ppm"))), "\n")
  
  cat("\n== A. LGCP (kppm, clik2, exponential) on a synthetic Poisson pattern ==\n")
  set.seed(SEED0 + 1)
  XP <- to_poly(rpoispp(truth$lamW)); cat("Synthetic Poisson pattern n =", npoints(XP), "\n")
  fitK <- try_("kppm LGCP clik2", kppm(XP, trend = formula_E3, clusters = "LGCP", model = "exponential",
                                       covariates = covs, method = "clik2"))
  if (!is.null(fitK)) {
    cat("class:", paste(class(fitK), collapse = ","), "\n")
    cat("names:", paste(names(fitK), collapse = ", "), "\n")
    bK <- try_("coef(kppm)", coef(fitK)); if (!is.null(bK)) { cat("length(coef):", length(bK), "\n"); print(round(bK, 3)) }
    for (nm in c("clusters", "clustpar", "clustargs", "modelpar", "par", "modelname")) {
      cat("fit[['", nm, "']]: ", sep = ""); print(fitK[[nm]]) }
    Vk <- try_("vcov(kppm)", vcov(fitK))
    if (!is.null(Vk) && !is.null(bK)) describe_vcov(Vk, bK, "vcov(kppm)")
    ak <- filter_args("simulate", "kppm", list(nsim = 2, drop = TRUE, verbose = FALSE))
    cat("simulate.kppm args passed:", paste(names(ak), collapse = ", "), "\n")
    Sk <- try_("simulate(kppm)", do.call(simulate, c(list(fitK), ak)))
    if (!is.null(Sk)) try_("describe kppm sims", describe_sims(Sk, "LGCP"))
  }
  
  cat("\n== B. Poisson ppm (fixed quadrature) on the same pattern ==\n")
  fitP <- try_("ppm Poisson", ppm(make_quad(XP), trend = formula_E3, covariates = covs))
  if (!is.null(fitP)) {
    bP <- coef(fitP); cat("length(coef):", length(bP), "\n")
    Vp <- try_("vcov(ppm Poisson)", vcov(fitP))
    if (!is.null(Vp)) describe_vcov(Vp, bP, "vcov(ppm Poisson)")
    ap <- filter_args("simulate", "ppm", list(nsim = 2, drop = TRUE, verbose = FALSE))
    cat("simulate.ppm args passed:", paste(names(ap), collapse = ", "), "\n")
    Sp <- try_("simulate(ppm Poisson)", do.call(simulate, c(list(fitP), ap)))
    if (!is.null(Sp)) try_("describe Poisson sims", describe_sims(Sp, "Poisson"))
  }
  
  cat("\n== C. Geyer (G3 truth, fit at r = 3, sat = 3) ==\n")
  XG <- try_("simulate G3 truth", sim_geyer(MULT_G3, 2.5, 3, 3, 197, 1, SEED0 + 3))
  if (!is.null(XG)) {
    cat("Geyer pattern n =", npoints(XG), "\n")
    fitG <- try_("geyer_fit r=3 sat=3", geyer_fit(XG, r = 3, sat = 3))
    if (!is.null(fitG)) {
      bG <- coef(fitG); cat("length(coef):", length(bG), "| last names:", paste(tail(names(bG), 3), collapse = ", "), "\n")
      cat("Interaction coef (log gamma):", format(tail(bG, 1), digits = 4), " exp:", format(exp(tail(bG, 1)), digits = 4), "\n")
      Vg <- try_("vcov(ppm Geyer)", vcov(fitG))
      if (!is.null(Vg)) describe_vcov(Vg, bG, "vcov(ppm Geyer)")
      drawG <- if (!is.null(Vg) && requireNamespace("MASS", quietly = TRUE)) try_("draw new.coef", draw_coef(bG, Vg, SEED0 + 4)) else NULL
      if (!is.null(drawG)) cat("draw vs coef, max abs diff:", format(max(abs(drawG - bG)), digits = 4), "\n")
      ag <- filter_args("simulate", "ppm", list(nsim = 2, drop = TRUE, verbose = FALSE,
                                                control = rmhcontrol(nrep = NREP_PROBE)))
      cat("simulate.ppm args passed (no new.coef):", paste(names(ag), collapse = ", "), "\n")
      Sg0 <- try_("simulate(ppm Geyer) point est", do.call(simulate, c(list(fitG), ag)))
      if (!is.null(Sg0)) try_("describe Geyer sims (point)", describe_sims(Sg0, "Geyer point"))
      if (!is.null(drawG)) {
        ag1 <- filter_args("simulate", "ppm", list(nsim = 2, drop = TRUE, verbose = FALSE, new.coef = drawG,
                                                   control = rmhcontrol(nrep = NREP_PROBE)))
        cat("simulate.ppm args passed (new.coef):", paste(names(ag1), collapse = ", "), "\n")
        Sg1 <- try_("simulate(ppm Geyer) new.coef", do.call(simulate, c(list(fitG), ag1)))
        if (!is.null(Sg1)) try_("describe Geyer sims (draw)", describe_sims(Sg1, "Geyer draw"))
      }
    }
  }
  cat("\n== END PROBE ==\n")
}, finally = { while (sink.number() > 0) sink() })

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "25_sessionInfo.txt"))
cat("Done. Report:", REPORT, "\n")