# =============================================================================
# 23_log_window_fallback.R
# Purpose: log the activation of the pre-registered fallback window (freeze v1) and the
#   extension of the required-settings list. Appends one row to docs/protocol_deviation_log.csv.
# Not executed by assistant. Run once.
# =============================================================================
REPO_CONFIG <- "E:/Projects & Researches/NC_PointProcess_Repo/R/config.R"
while (sink.number() > 0) sink()
source(REPO_CONFIG)
dl_f <- file.path(REPO_DIR, "docs", "protocol_deviation_log.csv")
dl <- read.csv(dl_f, stringsAsFactors = FALSE)
new <- data.frame(
  date = format(Sys.Date()),
  item = "ACTIVATION of pre-registered fallback window (100 m simplified polygon) for ALL simulation-based tests (recovery and empirical confirmatory); exact window becomes a REQUIRED sensitivity setting",
  reason = "Script 22: per-pattern summary-statistic cost 0.655 s (exact, 18,926 vertices) vs 0.275 s (simplified, 1,378 vertices), x2.38; curve step 0.415 s vs 0.045 s. Script 21: simplified window area change -0.0009%, symmetric difference 28.9 km2, all 197 sites inside. Script 08: K(r) change below 7e-5 for r >= 1 km.",
  effect_on_inference = "Negligible geometric change (area -0.0009%). Interaction claim must also hold on the exact window (added to DECISION$required_settings: exact_window). Makes the claim rule stricter, not looser.",
  approved_by = "author-delegated (assistant)", stringsAsFactors = FALSE)
if (!new$item %in% dl$item) { write.csv(rbind(dl[, names(new)], new), dl_f, row.names = FALSE); cat("Appended fallback activation to docs/protocol_deviation_log.csv\n") } else cat("Row already present\n")