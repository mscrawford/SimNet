#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

source(paste0(scripts_dir, "/to_test/randomForest_4conditions.R"))

### Forest1 (PPA)
source(paste0(scripts_dir, "/to_test/readPPA.R"))

PPA <- models$Forest1 %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))%>%
  select(-Model, -SeedRain)

F1C1 <- randomForest_4conditions("F1C1",PPA,1,meta)
F1C2 <- randomForest_4conditions("F1C2",PPA,1,iso)
F1C3 <- randomForest_4conditions("F1C3",PPA,32,meta)
F1C4 <- randomForest_4conditions("F1C4",PPA,32,iso)
