#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

source(paste0(scripts_dir, "/to_test/randomForest_4conditions.R"))

### Grass3 (IBC-grass)
source(paste0(scripts_dir, "/to_test/readIBC.R"))

IBC_grass <- models$Grass3 %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))%>%
  select(-Model, -SeedRain) 

G3C1 <- randomForest_4conditions("G3C1",IBC_grass,1,meta)
G3C2 <- randomForest_4conditions("G3C2",IBC_grass,1,iso)
G3C3 <- randomForest_4conditions("G3C3",IBC_grass,32,meta)
G3C4 <- randomForest_4conditions("G3C4",IBC_grass,32,iso)
