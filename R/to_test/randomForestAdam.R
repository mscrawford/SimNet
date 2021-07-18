#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

source(paste0(scripts_dir, "/to_test/randomForest_4conditions.R"))

### Grass1 (Adam's model)
source(paste0(scripts_dir, "/to_test/readAdam.R"))

adam <- models$Grass1 %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))%>%
  select(-Model, -SeedRain) %>%
  mutate(SpeciesID = as.factor(SpeciesID))

G1C1 <- randomForest_4conditions("G1C1",adam,1,meta)
G1C2 <- randomForest_4conditions("G1C2",adam,1,iso)
G1C3 <- randomForest_4conditions("G1C3",adam,32,meta)
G1C4 <- randomForest_4conditions("G1C4",adam,32,iso)
