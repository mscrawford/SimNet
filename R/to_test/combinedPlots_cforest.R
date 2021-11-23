set.seed(1987)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

#source(paste0(scripts_dir, "/readModels.R"))
source(paste0(scripts_dir, "/to_test/fx_cforest_party.R"))

### Grass1 (Adam's model)
source(paste0(scripts_dir, "/to_test/readAdam.R"))

adam <- models$Grass1 %>%
 ungroup() %>% # There shouldn't be groups anyways
 mutate(Biomass = scales::rescale(Biomass,
                                  to = c(0, 100)))%>%
 #select(-Model, -SeedRain, -abmi) #Removed abmi because it is highly correlated to biomass
 select(-Model, -SeedRain)

fx_cforest(adam,"Grass1")

### Grass2 (Lindsay's model)
source(paste0(scripts_dir, "/to_test/readLindsay.R"))

lindsay <- models$Grass2%>%
   ungroup() %>% # There shouldn't be groups anyways
   mutate(Biomass = scales::rescale(Biomass,
                                    to = c(0, 100)))%>%
   select(-Model, -SeedRain)

fx_cforest(lindsay,"Grass2")

### Grass3 (IBC-grass)
source(paste0(scripts_dir, "/to_test/readIBC.R"))

IBC_grass <- models$Grass3 %>%
 ungroup() %>% # There shouldn't be groups anyways
 mutate(Biomass = scales::rescale(Biomass,
                                  to = c(0, 100)))%>%
 select(-Model, -SeedRain) %>%
 select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)

fx_cforest(IBC_grass,"Grass3")

### Forest1 (PPA)
source(paste0(scripts_dir, "/to_test/readPPA.R"))

PPA <- models$Forest1 %>%
 ungroup() %>% # There shouldn't be groups anyways
 mutate(Biomass = scales::rescale(Biomass,
                                  to = c(0, 100)))%>%
 select(-Model, -SeedRain)

fx_cforest(PPA,"Forest1")

### Forest2 (TROLL)
source(paste0(scripts_dir, "/to_test/readTROLL.R"))

troll <- models$Forest2 %>%
   ungroup() %>% # There shouldn't be groups anyways
   mutate(Biomass = scales::rescale(Biomass,
                                    to = c(0, 100)))%>%
   select(-Model, -SeedRain)

fx_cforest(troll,"Forest2")

### Forest2 (TROLL) h_realmax
#source(paste0(scripts_dir, "/to_test/readTROLL.R"))

troll <- models$Forest2 %>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))%>%
    select(-Model, -SeedRain) %>%
    mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
    select(-hmax, -ah, -dmax)

fx_cforest(troll,"Forest2_hrealmax")

### Dryland (Bjoern)
source(paste0(scripts_dir, "/to_test/readBjoern.R"))

bjoern <- models$bjoern %>% 
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))%>%
    select(-Model, -SeedRain) 

fx_cforest(bjoern,"Dryland")
