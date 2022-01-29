library(cowplot)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

#source(paste0(scripts_dir, "/readModels.R"))
source(paste0(scripts_dir, "/to_test/fx_traits_vs_biomass.R"))

### Grass3 (IBC-grass)
source(paste0(scripts_dir, "/to_test/readIBC.R"))

IBC_grass <- models$Grass3 %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))%>%
  select(-Model, -SeedRain) %>%
  select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)
print(IBC_grass)	

g_by <- c("SpeciesID", "Year", "Stage", "LMR", "MaxMass", "Gmax", "SLA", "meanSpacerLength")
lab1 <- "Gmax (maximum resource \n utilization per time step)"
lab2 <- "MaxMass (Plant's \n maximum size)"
lab3 <- "LMR (leaf \n to mass ratio)"
lab4 <- "MeanSpacerLength"

#G3C1 <- fx_traits_vs_biomass_jitter("G3C1",IBC_grass,1,meta,g_by,"Gmax","MaxMass",lab1,lab2)
#G3C2 <- fx_traits_vs_biomass_jitter("G3C2",IBC_grass,1,iso,g_by,"Gmax","MaxMass",lab1,lab2)
#G3C3 <- fx_traits_vs_biomass_jitter("G3C3",IBC_grass,32,meta,g_by,"Gmax","MaxMass",lab1,lab2)
#G3C4 <- fx_traits_vs_biomass_jitter("G3C4",IBC_grass,32,iso,g_by,"Gmax","meanSpacerLength",lab1,lab4)
Gmax <- fx_traits_vs_biomass_jitter("Gmax_Biomass",IBC_grass,32,meta,g_by,"Gmax","Biomass",lab1,"Log mean Biomass")
MaxMass <- fx_traits_vs_biomass_jitter("MaxMass_Biomass",IBC_grass,32,meta,g_by,"MaxMass","Biomass",lab2,"Log mean Biomass")

### Grass3 (IBC-grass) noNDD
source(paste0(scripts_dir, "/to_test/readIBC.R"))

IBC_grass_noNDD <- models$IBC_grass.noNDD %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))%>%
  select(-Model, -SeedRain) %>%
  select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)
print(IBC_grass_noNDD)	

g_by <- c("SpeciesID", "Year", "Stage", "LMR", "MaxMass", "Gmax", "SLA", "meanSpacerLength")
lab1 <- "Gmax (maximum resource \n utilization per time step)"
lab2 <- "MaxMass (Plant's \n maximum size)"
lab3 <- "LMR (leaf \n to mass ratio)"
lab4 <- "MeanSpacerLength"

#G3C1 <- fx_traits_vs_biomass_jitter("G3C1_noNDD",IBC_grass_noNDD,1,meta,g_by,"Gmax","MaxMass",lab1,lab2)
#G3C2 <- fx_traits_vs_biomass_jitter("G3C2_noNDD",IBC_grass_noNDD,1,iso,g_by,"Gmax","MaxMass",lab1,lab2)
#G3C3 <- fx_traits_vs_biomass_jitter("G3C3_noNDD",IBC_grass_noNDD,32,meta,g_by,"Gmax","MaxMass",lab1,lab2)
#G3C4 <- fx_traits_vs_biomass_jitter("G3C4_noNDD",IBC_grass_noNDD,32,iso,g_by,"Gmax","meanSpacerLength",lab1,lab4)
Gmax_noNDD <- fx_traits_vs_biomass_jitter("Gmax_Biomass_noNDD",IBC_grass_noNDD,32,meta,g_by,"Gmax","Biomass",lab1,"Log mean Biomass")
MaxMass_noNDD <- fx_traits_vs_biomass_jitter("MaxMass_Biomass_noNDD",IBC_grass_noNDD,32,meta,g_by,"MaxMass","Biomass",lab2,"Log mean Biomass")

plot_grid(Gmax, MaxMass, Gmax_noNDD, MaxMass_noNDD,
         labels = c("Gmax vs. Biomass","MaxMass vs. Biomass","Gmax vs. Biomass (noNDD)","MaxMass vs. Biomass (noNDD)",base_asp = 2,
       nrow = 2,
       ncol = 2))#,
#         ncol = 2, nrow = 2, widths = 1, heights = 1,)
ggsave(paste0(tmp_dir,"/traits_vs_biomass/Gmax_MaxMass_NDD_noNDD.pdf"))
