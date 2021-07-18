#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

#source(paste0(scripts_dir, "/readModels.R"))
source(paste0(scripts_dir, "/to_test/randomForest_4conditions.R"))

### Grass1 (Adam's model)
source(paste0(scripts_dir, "/to_test/readAdam.R"))

adam <- models$Grass1 %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))%>%
  select(-Model, -SeedRain) 

G1C1 <- randomForest_4conditions("G1C1",adam,1,meta)
G1C2 <- randomForest_4conditions("G1C2",adam,1,iso)
G1C3 <- randomForest_4conditions("G1C3",adam,32,meta)
G1C4 <- randomForest_4conditions("G1C4",adam,32,iso)

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

### Grass3 (IBC-grass)
source(paste0(scripts_dir, "/to_test/readIBC.R"))

IBC_grass <- models$Grass3 %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))%>%
  select(-Model, -SeedRain) %>%
  select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)
 

G3C1 <- randomForest_4conditions("G3C1",IBC_grass,1,meta)
G3C2 <- randomForest_4conditions("G3C2",IBC_grass,1,iso)
G3C3 <- randomForest_4conditions("G3C3",IBC_grass,32,meta)
G3C4 <- randomForest_4conditions("G3C4",IBC_grass,32,iso)

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

### Dryland (Bjoern)
source(paste0(scripts_dir, "/to_test/readBjoern.R"))

bjoern <- models$bjoern %>% 
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))%>%
    select(-Model, -SeedRain) 

DC1 <- randomForest_4conditions("DC1",bjoern,1,meta)
DC2 <- randomForest_4conditions("DC2",bjoern,1,iso)
DC3 <- randomForest_4conditions("DC3",bjoern,32,meta)
DC4 <- randomForest_4conditions("DC4",bjoern,32,iso)

