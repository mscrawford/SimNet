library(tidyverse)

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
# Focal metrics

# X_VALS <- c("Shannon", "Richness", "FDis")
X_VALS <- c("Shannon")

# Y_VALS <- c("biomass", "productivity")
Y_VALS <- c("biomass")


# -------------------------------------------------------------------------
# Run analysis

purrr::walk(.x = Y_VALS,
            .f = ~ {
                assign(x = "y_val", value = .x, envir = .GlobalEnv)
                source(paste0(scripts_dir, "/SimNet_ReadModels.R"))
                source(paste0(scripts_dir, "/SimNet_BRMS.R"))
                source(paste0(scripts_dir, "/SimNet_BRMS_derivePlots.R"))
                source(paste0(scripts_dir, "/SimNet_FunctionDominanceCorr.R"))
            })
