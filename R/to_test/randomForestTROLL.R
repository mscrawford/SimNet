#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

source(paste0(scripts_dir, "/to_test/randomForest_4conditions.R"))

### Forest2 (TROLL)
source(paste0(scripts_dir, "/to_test/readTROLL.R"))

troll <- models$Forest2 %>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))%>%
    select(-Model, -SeedRain) %>%
    mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
    select(-hmax, -ah)

F2C1 <- randomForest_4conditions("F2C1",troll,1,meta)
F2C2 <- randomForest_4conditions("F2C2",troll,1,iso)
F2C3 <- randomForest_4conditions("F2C3",troll,32,meta)
F2C4 <- randomForest_4conditions("F2C4",troll,32,iso)
