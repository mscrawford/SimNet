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

# fx_cforest_party("t1t",adam,1,meta)
# fx_cforest_party("t1",adam,1,meta)
# fx_cforest_party("t1m",adam,32,meta)
G1C1 <- fx_cforest_party("G1C1",adam,1,meta)
G1C2 <- fx_cforest_party("G1C2",adam,1,iso)
G1C3 <- fx_cforest_party("G1C3",adam,32,meta)
G1C4 <- fx_cforest_party("G1C4",adam,32,iso)

### Grass2 (Lindsay's model)
source(paste0(scripts_dir, "/to_test/readLindsay.R"))

lindsay <- models$Grass2%>%
   ungroup() %>% # There shouldn't be groups anyways
   mutate(Biomass = scales::rescale(Biomass,
                                    to = c(0, 100)))%>%
   select(-Model, -SeedRain)

# fx_cforest_party("t2t",lindsay,1,meta)
# fx_cforest_party("t2",lindsay,1,meta)
# fx_cforest_party("t2m",lindsay,32,meta)
G2C1 <- fx_cforest_party("G2C1",lindsay,1,meta)
G2C2 <- fx_cforest_party("G2C2",lindsay,1,iso)
G2C3 <- fx_cforest_party("G2C3",lindsay,32,meta)
G2C4 <- fx_cforest_party("G2C4",lindsay,32,iso)

### Grass3 (IBC-grass)
source(paste0(scripts_dir, "/to_test/readIBC.R"))

IBC_grass <- models$Grass3 %>%
 ungroup() %>% # There shouldn't be groups anyways
 mutate(Biomass = scales::rescale(Biomass,
                                  to = c(0, 100)))%>%
 select(-Model, -SeedRain) %>%
 select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)

# fx_cforest_party("t3t",IBC_grass,1,meta)
# fx_cforest_party("t3",IBC_grass,1,meta)
# fx_cforest_party("t3m",IBC_grass,32,meta)
G3C1 <- fx_cforest_party("G3C1",IBC_grass,1,meta)
G3C2 <- fx_cforest_party("G3C2",IBC_grass,1,iso)
G3C3 <- fx_cforest_party("G3C3",IBC_grass,32,meta)
G3C4 <- fx_cforest_party("G3C4",IBC_grass,32,iso)

### Forest1 (PPA)
source(paste0(scripts_dir, "/to_test/readPPA.R"))

PPA <- models$Forest1 %>%
 ungroup() %>% # There shouldn't be groups anyways
 mutate(Biomass = scales::rescale(Biomass,
                                  to = c(0, 100)))%>%
 select(-Model, -SeedRain)

# fx_cforest_party("t4t",PPA,1,meta)
# fx_cforest_party("t4",PPA,1,meta)
# fx_cforest_party("t4m",PPA,32,meta)
F1C1 <- fx_cforest_party("F1C1",PPA,1,meta)
F1C2 <- fx_cforest_party("F1C2",PPA,1,iso)
F1C3 <- fx_cforest_party("F1C3",PPA,32,meta)
F1C4 <- fx_cforest_party("F1C4",PPA,32,iso)

### Forest2 (TROLL)
source(paste0(scripts_dir, "/to_test/readTROLL.R"))

troll <- models$Forest2 %>%
   ungroup() %>% # There shouldn't be groups anyways
   mutate(Biomass = scales::rescale(Biomass,
                                    to = c(0, 100)))%>%
   select(-Model, -SeedRain)

# fx_cforest_party("t5t",troll,1,meta)
# fx_cforest_party("t5",troll,1,meta)
# fx_cforest_party("t5m",troll,32,meta)
F2C1 <- fx_cforest_party("F2C1",troll,1,meta)
F2C2 <- fx_cforest_party("F2C2",troll,1,iso)
F2C3 <- fx_cforest_party("F2C3",troll,32,meta)
F2C4 <- fx_cforest_party("F2C4",troll,32,iso)

### Forest2 (TROLL) h_realmax
#source(paste0(scripts_dir, "/to_test/readTROLL.R"))

troll <- models$Forest2 %>%
   ungroup() %>% # There shouldn't be groups anyways
   mutate(Biomass = scales::rescale(Biomass,
                                    to = c(0, 100)))%>%
   select(-Model, -SeedRain) %>%
   mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
   select(-hmax, -ah, -dmax)

# fx_cforest_party("t5ht",troll,1,meta)
# fx_cforest_party("t5h",troll,1,meta)
# fx_cforest_party("t5mh",troll,32,meta)
F2C1 <- fx_cforest_party("F2C1_hrm",troll,1,meta)
F2C2 <- fx_cforest_party("F2C2_hrm",troll,1,iso)
F2C3 <- fx_cforest_party("F2C3_hrm",troll,32,meta)
F2C4 <- fx_cforest_party("F2C4_hrm",troll,32,iso)

### Dryland (Bjoern)
source(paste0(scripts_dir, "/to_test/readBjoern.R"))

bjoern <- models$bjoern %>% 
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))%>%
    select(-Model, -SeedRain) 


# fx_cforest_party("t6t",bjoern,1,meta)
# fx_cforest_party("t6",bjoern,1,meta)
# fx_cforest_party("t6m",bjoern,32,meta)
DC1 <- fx_cforest_party("DC1",bjoern,1,meta)
DC2 <- fx_cforest_party("DC2",bjoern,1,iso)
DC3 <- fx_cforest_party("DC3",bjoern,32,meta)
DC4 <- fx_cforest_party("DC4",bjoern,32,iso)
