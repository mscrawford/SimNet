
# -------------------------------------------------------------------------
# Directories

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../")
base_dir          <- getwd()
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")
model_data_dir    <- paste0(base_dir, "/data/brms_models")


# -------------------------------------------------------------------------
# Options

READ_CACHE <- FALSE
SAVE_CACHE <- TRUE


# -------------------------------------------------------------------------
# Focal biodiversity metrics

# MEASURES <- c("Shannon", "Richness", "FDis")
MEASURES <- c("Shannon")

# -------------------------------------------------------------------------
# Run scripts

# Derive dataset and plots for biomass
Y_VAL = "biomass"
source(paste0(scripts_dir, "/SimNet_ReadModels.R"))
source(paste0(scripts_dir, "/SimNet_BRMS.R"))
source(paste0(scripts_dir, "/SimNet_BRMS_derivePlots.R"))
source(paste0(scripts_dir, "/SimNet_Correlations_vs_Estimates.R"))

# Derive dataset and plots for productivity
Y_VAL = "productivity"
source(paste0(scripts_dir, "/SimNet_ReadModels.R"))
source(paste0(scripts_dir, "/SimNet_BRMS.R"))
source(paste0(scripts_dir, "/SimNet_BRMS_derivePlots.R"))
source(paste0(scripts_dir, "/SimNet_Correlations_vs_Estimates.R"))
