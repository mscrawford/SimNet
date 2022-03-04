#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

source(paste0(scripts_dir, "/to_test/fx_PCA_all_models.R"))

### Grass3 (IBC-grass)
source(paste0(scripts_dir, "/to_test/readIBC.R"))

IBC_grass <- models$Grass3 %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))%>%
  select(-Model, -SeedRain) %>%
  select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)
 
G3PCA <- fx_PCA("Grass3",IBC_grass)
fviz_eig(G3PCA, addlabels = TRUE, ylim = c(0, 50))
ggsave(paste0(tmp_dir,"/PCA/ScreePlot_Grass3.pdf")) #Eigenvals organized from largest to smallest
while (!is.null(dev.list()))  dev.off()

modelName <- "Grass3"
G3PCA_corr <- fx_correlation_plot(modelName,G3PCA)

ind <- get_pca_ind(G3PCA)
df <- ind$coord

IBC_grass$PC1score <- df[,1]
IBC_grass$PC2score <- df[,2]
IBC_grass$PC3score <- df[,3]
fileName = paste0(tmp_dir,"/PCA/Grass3_PCAcoord.Rda")
saveRDS(IBC_grass,file=fileName)

#pcr_model <- pcr(Biomass~LMR+MaxMass+Gmax+SLA+meanSpacerLength, data=IBC_grass, scale =TRUE, validation="CV")
#summary(pcr_model)

### Forest2 (TROLL)
source(paste0(scripts_dir, "/to_test/readTROLL.R"))

troll <- models$Forest2 %>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))%>%
    select(-Model, -SeedRain)
modelName <- "Forest2"
F2PCA <- fx_PCA(modelName,troll)
fviz_eig(F2PCA, addlabels = TRUE, ylim = c(0, 50))
ggsave(paste0(tmp_dir,"/PCA/ScreePlot_Forest2.pdf")) #Eigenvals organized from largest to smallest
while (!is.null(dev.list()))  dev.off()

F2PCA_corr <- fx_correlation_plot(modelName,F2PCA)

F2PCR <- pcr(Productivity~hmax+dmax+ah+wsg+lma+pmass+nmass, data=troll, scale =TRUE, validation="CV")
#summary(pcr_model)
ind <- get_pca_ind(F2PCA)
df <- ind$coord

troll$PC1score <- df[,1]
troll$PC2score <- df[,2]
troll$PC3score <- df[,3]
fileName = paste0(tmp_dir,"/PCA/Forest2_PCAcoord.Rda")
saveRDS(troll,file=fileName)

#pcr_model <- pcr(Biomass~h_realmax+wsg+lma+pmass+nmass, data=troll, scale =TRUE, validation="CV")
#summary(pcr_model)

### Forest2 (TROLL)
source(paste0(scripts_dir, "/to_test/readTROLL.R"))

troll <- models$Forest2 %>%
   ungroup() %>% # There shouldn't be groups anyways
   mutate(Biomass = scales::rescale(Biomass,
                                    to = c(0, 100)))%>%
   select(-Model, -SeedRain) %>%
   mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
   select(-hmax, -ah, -dmax)

F2PCA_hrm <- fx_PCA("Forest2_hrm",troll)
fviz_eig(F2PCA, addlabels = TRUE, ylim = c(0, 50))
ggsave(paste0(tmp_dir,"/PCA/ScreePlot_Forest2_hrm.pdf")) #Eigenvals organized from largest to smallest
while (!is.null(dev.list()))  dev.off()

modelName <- "Forest2_hrm"
F2PCA_hrm_corr <- fx_correlation_plot(modelName,F2PCA_hrm)

F2PCR <- pcr(Productivity~h_realmax+wsg+lma+pmass+nmass, data=troll, scale =TRUE, validation="CV")
#summary(pcr_model)
ind <- get_pca_ind(F2PCA)
df <- ind$coord

troll$PC1score <- df[,1]
troll$PC2score <- df[,2]
troll$PC3score <- df[,3]
fileName = paste0(tmp_dir,"/PCA/Forest2_hrm_PCAcoord.Rda")
saveRDS(troll,file=fileName)

#pcr_model <- pcr(Biomass~h_realmax+wsg+lma+pmass+nmass, data=troll, scale =TRUE, validation="CV")
#summary(pcr_model)
