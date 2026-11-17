# config.R -- sourced by every script
RAW_DIR  <- "E:/Projects & Researches/Nagpur-Chandrapur Enhanced Dataset"
REPO_DIR <- "E:/Projects & Researches/NC_PointProcess_Repo"
CRS_PROJ <- 32644        # EPSG:32644 WGS84/UTM 44N  -- PROPOSED, confirm at freeze
MASTER_SEED <- 20260101  # PROPOSED, confirm at freeze; batch seeds derived from this
set.seed(MASTER_SEED)
dir_of <- function(stage, sub = 'outputs') file.path(REPO_DIR, stage, sub)
