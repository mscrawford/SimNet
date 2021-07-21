#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

source(paste0(scripts_dir, "/to_test/randomForest_4conditions.R"))

### Grass2 (Lindsay's model)
source(paste0(scripts_dir, "/to_test/readLindsay.R"))

lindsay <- models$Grass2%>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))%>%
    select(-Model, -SeedRain) 

G2C1 <- randomForest_4conditions("G2C1",lindsay,1,meta)
G2C2 <- randomForest_4conditions("G2C2",lindsay,1,iso)
G2C3 <- randomForest_4conditions("G2C3",lindsay,32,meta)
G2C4 <- randomForest_4conditions("G2C4",lindsay,32,iso)
