#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

#source(paste0(scripts_dir, "/readModels.R"))
source(paste0(scripts_dir, "/to_test/fx_traits_vs_biomass.R"))

### Grass1 (Adam's model)
source(paste0(scripts_dir, "/to_test/readAdam.R"))

adam <- models$Grass1 %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))%>%
  #select(-Model, -SeedRain, -abmi) #Removed abmi because it is highly correlated to biomass
  select(-Model, -SeedRain)

g_by <- c("SpeciesID", "Year", "Stage","abmi", "pNi", "no3i")
lab1 <- "no3i (nitrogen R*)"
lab2 <- "pNi (aboveground N concentration)"

G1C1_a <- fx_traits_vs_biomass_jitter("G1C1_a",adam,1,meta,g_by,"abmi","no3i","abmi",lab1)
G1C2_a <- fx_traits_vs_biomass_jitter("G1C2_a",adam,1,iso,g_by,"abmi","no3i","abmi",lab1)
G1C3_a <- fx_traits_vs_biomass_jitter("G1C3_a",adam,32,meta,g_by,"abmi","no3i","abmi",lab1)
G1C4_a <- fx_traits_vs_biomass_jitter("G1C4_a",adam,32,iso,g_by,"abmi","no3i","abmi",lab1)
G1C1_b <- fx_traits_vs_biomass_jitter("G1C1_b",adam,1,meta,g_by,"abmi","pNi","abmi",lab2)
G1C2_b <- fx_traits_vs_biomass_jitter("G1C2_b",adam,1,iso,g_by,"abmi","pNi","abmi",lab2)
G1C3_b <- fx_traits_vs_biomass_jitter("G1C3_b",adam,32,meta,g_by,"abmi","pNi","abmi",lab2)
G1C4_b <- fx_traits_vs_biomass_jitter("G1C4_b",adam,32,iso,g_by,"abmi","pNi","abmi",lab2)
#abmi2 <- fx_traits_vs_biomass_jitter("abmi2",adam,1,meta,g_by,"abmi","pNi","abmi",lab2)

G1C1 <- fx_traits_vs_biomass_jitter("G1C1",adam,1,meta,g_by,"no3i","pNi",lab1,lab2)
G1C2 <- fx_traits_vs_biomass_jitter("G1C2",adam,1,iso,g_by,"no3i","pNi",lab1,lab2)
G1C3 <- fx_traits_vs_biomass_jitter("G1C3",adam,32,meta,g_by,"no3i","pNi",lab1,lab2)
G1C4 <- fx_traits_vs_biomass_jitter("G1C4",adam,32,iso,g_by,"no3i","pNi",lab1,lab2)

#### Grass2 (Lindsay's model)
#source(paste0(scripts_dir, "/to_test/readLindsay.R"))
#
#lindsay <- models$Grass2%>%
#    ungroup() %>% # There shouldn't be groups anyways
#    mutate(Biomass = scales::rescale(Biomass,
#                                     to = c(0, 100)))%>%
#    select(-Model, -SeedRain) 
#
#g_by <- c("SpeciesID", "Year", "Stage", "Vi", "thetai")
#lab1 <- "Vi (volume of soil \n accessible to species i)"
#lab2 <- "thetai (Nitrogen uptake \n rate per unit plant biomass)"
#
#G2C1 <- fx_traits_vs_biomass_jitter("G2C1",lindsay,1,meta,g_by,"Vi","thetai",lab1,lab2)
#G2C2 <- fx_traits_vs_biomass_jitter("G2C2",lindsay,1,iso,g_by,"Vi","thetai",lab1,lab2)
#G2C3 <- fx_traits_vs_biomass_jitter("G2C3",lindsay,32,meta,g_by,"Vi","thetai",lab1,lab2)
#G2C4 <- fx_traits_vs_biomass_jitter("G2C4",lindsay,32,iso,g_by,"Vi","thetai",lab1,lab2)
#
#### Grass3 (IBC-grass)
#source(paste0(scripts_dir, "/to_test/readIBC.R"))
#
#IBC_grass <- models$Grass3 %>%
#  ungroup() %>% # There shouldn't be groups anyways
#  mutate(Biomass = scales::rescale(Biomass,
#                                   to = c(0, 100)))%>%
#  select(-Model, -SeedRain) %>%
#  select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)
# 
#g_by <- c("SpeciesID", "Year", "Stage", "LMR", "MaxMass", "Gmax", "SLA", "meanSpacerLength")
#lab1 <- "Gmax (maximum resource \n utilization per time step)"
#lab2 <- "MaxMass (Plant's \n maximum size)"
#lab3 <- "LMR (leaf \n to mass ratio)"
#lab4 <- "MeanSpacerLength"
#
#G3C1 <- fx_traits_vs_biomass_jitter("G3C1",IBC_grass,1,meta,g_by,"Gmax","MaxMass",lab1,lab2)
#G3C2 <- fx_traits_vs_biomass_jitter("G3C2",IBC_grass,1,iso,g_by,"Gmax","MaxMass",lab1,lab2)
#G3C3 <- fx_traits_vs_biomass_jitter("G3C3",IBC_grass,32,meta,g_by,"Gmax","MaxMass",lab1,lab2)
#G3C4 <- fx_traits_vs_biomass_jitter("G3C4",IBC_grass,32,iso,g_by,"Gmax","meanSpacerLength",lab1,lab4)
#
#### Forest1 (PPA)
#source(paste0(scripts_dir, "/to_test/readPPA.R"))
#
#PPA <- models$Forest1 %>%
#  ungroup() %>% # There shouldn't be groups anyways
#  mutate(Biomass = scales::rescale(Biomass,
#                                   to = c(0, 100)))%>%
#  select(-Model, -SeedRain)
#
#g_by <- c("SpeciesID", "Year", "Stage", "PC1score", "PC2score")
#lab1 <- "PC1score (associated \n with fast-slow lifecycle)"#plant height)"
#lab2 <- "PC2score (associated \n with tree stature)"# LMA -leaf mass per area)"
#
#F1C1 <- fx_traits_vs_biomass_jitter("F1C1",PPA,1,meta,g_by,"PC2score","PC1score",lab2,lab1)
#F1C2 <- fx_traits_vs_biomass_jitter("F1C2",PPA,1,iso,g_by,"PC2score","PC1score",lab2,lab1)
#F1C3 <- fx_traits_vs_biomass_jitter("F1C3",PPA,32,meta,g_by,"PC2score","PC1score",lab2,lab1)
#F1C4 <- fx_traits_vs_biomass_jitter("F1C4",PPA,32,iso,g_by,"PC2score","PC1score",lab2,lab1)
#
#### Forest2 (TROLL)
#source(paste0(scripts_dir, "/to_test/readTROLL.R"))
#
#troll <- models$Forest2 %>%
#    ungroup() %>% # There shouldn't be groups anyways
#    mutate(Biomass = scales::rescale(Biomass,
#                                     to = c(0, 100)))%>%
#    select(-Model, -SeedRain)
#
#g_by <- c("SpeciesID", "Year","Stage", "dmax", "wsg", "lma", "pmass", "nmass", "hmax", "ah")
#lab1 <- "lma (leaf \n mass per area)"
#lab2 <- "dmax (diameter at \n breast height threshold)"
#lab3 <- "wsg (wood \n specific gravity)"
#lab4 <- "P - pmass (leaf P \n content per dry mass)"
#
#F2C1 <- fx_traits_vs_biomass_jitter("F2C1",troll,1,meta,g_by,"dmax","wsg",lab2,lab3)
#F2C2 <- fx_traits_vs_biomass_jitter("F2C2",troll,1,iso,g_by,"dmax","wsg",lab2,lab3)
#F2C3 <- fx_traits_vs_biomass_jitter("F2C3",troll,32,meta,g_by,"lma","dmax",lab1,lab2)
#F2C4 <- fx_traits_vs_biomass_jitter("F2C4",troll,32,iso,g_by,"pmass","lma",lab4,lab1)
#
#### Forest2 (TROLL) h_realmax
#source(paste0(scripts_dir, "/to_test/readTROLL.R"))
#
#troll <- models$Forest2 %>%
#    ungroup() %>% # There shouldn't be groups anyways
#    mutate(Biomass = scales::rescale(Biomass,
#                                     to = c(0, 100)))%>%
#    select(-Model, -SeedRain) %>%
#    mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
#    select(-hmax, -ah, -dmax)
#
#g_by <- c("SpeciesID", "Year","Stage", "wsg", "lma", "pmass", "nmass", "h_realmax")
#lab1 <- "lma (leaf \n mass per area)"
#lab2 <- "h_realmax = \n hmax * dmax / (dmax + ah)"
#lab3 <- "wsg (wood \n specific gravity)"
#lab4 <- "P - pmass (leaf P \n content per dry mass)"
#
#F2C1 <- fx_traits_vs_biomass_jitter("F2C1_hrm",troll,1,meta,g_by,"h_realmax","wsg",lab2,lab3)
#F2C2 <- fx_traits_vs_biomass_jitter("F2C2_hrm",troll,1,iso,g_by,"h_realmax","wsg",lab2,lab3)
#F2C3 <- fx_traits_vs_biomass_jitter("F2C3_hrm",troll,32,meta,g_by,"lma","h_realmax",lab1,lab2)
#F2C4 <- fx_traits_vs_biomass_jitter("F2C4_hrm",troll,32,iso,g_by,"lma","pmass",lab1,lab4)
#
#### Dryland (Bjoern)
#source(paste0(scripts_dir, "/to_test/readBjoern.R"))
#
#bjoern <- models$bjoern %>% 
#    ungroup() %>% # There shouldn't be groups anyways
#    mutate(Biomass = scales::rescale(Biomass,
#                                     to = c(0, 100)))%>%
#    select(-Model, -SeedRain) 
#
#g_by <- c("SpeciesID", "Year","Stage", "maxSize", "pLeaf", "pRoot", "pStorage")
#lab1 <- "maxSize (maximum \n size/size at maturity) [gC]"
#lab2 <- "pLeaf (allocation \n to leaf) [gC/gC]"
#lab3 <- "pStorage (allocation \n to storage) [gC/gC]"
#lab4 <- "pRoot (allocation \n to root) [gC/gC]"
#
#DC1 <- fx_traits_vs_biomass_jitter("DC1",bjoern,1,meta,g_by,"maxSize","pRoot",lab1,lab4)
#DC2 <- fx_traits_vs_biomass_jitter("DC2",bjoern,1,iso,g_by,"maxSize","pStorage",lab1,lab3)
#DC3 <- fx_traits_vs_biomass_jitter("DC3",bjoern,32,meta,g_by,"maxSize","pLeaf",lab1,lab2)
#DC4 <- fx_traits_vs_biomass_jitter("DC4",bjoern,32,iso,g_by,"maxSize","pLeaf",lab1,lab2)
