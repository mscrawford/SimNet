library(cowplot)
set.seed(1987)
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

source(paste0(scripts_dir, "/to_test/fx_cforest_party.R"))

### Dryland (Bjoern)
source(paste0(scripts_dir, "/to_test/readBjoern.R"))

bjoern <- models$bjoern %>% 
 select(-Model, -SeedRain)

#Biomass
model<-bjoern
model<-subset(model, select=-pRoot)
modelName = "Dryland"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
C1 <- fx_cforest_single_condition(modelName,model,1,meta)
C2 <- fx_cforest_single_condition(modelName,model,1,iso)
C3 <- fx_cforest_single_condition(modelName,model,32,meta)
C4 <- fx_cforest_single_condition(modelName,model,32,iso)

rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
#rf_df <- readRDS(fileName)
p1 <- fx_plot(rf_df,modelName)

#Productivity
modelName = "Dryland_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
C1 <- fx_cforest_single_condition(modelName,model,1,meta)
C2 <- fx_cforest_single_condition(modelName,model,1,iso)
C4 <- fx_cforest_single_condition(modelName,model,32,iso)
C3 <- fx_cforest_single_condition(modelName,model,32,meta)

rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
#rf_df <- readRDS(fileName)
p2 <- fx_plot(rf_df,modelName)
