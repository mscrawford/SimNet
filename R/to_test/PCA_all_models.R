library("pls")
library("factoextra")
library(FactoMineR)
library(corrplot)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

source(paste0(scripts_dir, "/to_test/fx_PCA_all_models.R"))

### Grass1 (Adam's model)
#source(paste0(scripts_dir, "/to_test/readAdam.R"))
#
#adam <- models$Grass1 %>%
#  ungroup() %>% # There shouldn't be groups anyways
#  mutate(Biomass = scales::rescale(Biomass,
#                                   to = c(0, 100)))%>%
#  #select(-Model, -SeedRain, -abmi) #Removed abmi because it is highly correlated to biomass
#  select(-Model, -SeedRain)
#
#G1C1 <- fx_PCA_all_models("Grass1",adam)
#
#### Grass2 (Lindsay's model)
#source(paste0(scripts_dir, "/to_test/readLindsay.R"))
#
#lindsay <- models$Grass2%>%
#    ungroup() %>% # There shouldn't be groups anyways
#    mutate(Biomass = scales::rescale(Biomass,
#                                     to = c(0, 100)))%>%
#    select(-Model, -SeedRain) 
#
#G2C1 <- fx_PCA_all_models("Grass2",lindsay)
#
### Grass3 (IBC-grass)
source(paste0(scripts_dir, "/to_test/readIBC.R"))

#IBC_grass <- models$Grass3 %>%
#  ungroup() %>% # There shouldn't be groups anyways
#  mutate(Biomass = scales::rescale(Biomass,
#                                   to = c(0, 100)))%>%
#  select(-Model, -SeedRain) %>%
#  select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)
# 
#
#G3C1 <- fx_PCA_all_models("Grass3",IBC_grass)
##pcr_model <- pcr(Biomass~LMR+MaxMass+Gmax+SLA+meanSpacerLength, data=IBC_grass, scale =TRUE, validation="CV")
##summary(pcr_model)
#fviz_eig(G3C1, addlabels = TRUE, ylim = c(0, 50))
#ggsave(paste0(tmp_dir,"/PCA/ScreePlot_Grass3.pdf")) #Eigenvals organized from largest to smallest
#while (!is.null(dev.list()))  dev.off()
#
#var <- get_pca_var(G3C1)
#corrplot(var$contrib, is.corr=FALSE)
#recordPlot()
##ggsave(paste0(tmp_dir,"/PCA/Corr_Grass3.pdf")) #Highlight most contributing variables to each dimension
##while (!is.null(dev.list()))  dev.off()
#
##res.desc <- dimdesc(G3C1, axes = c(1,2), proba = 0.05)
### Description of dimension 1
##print(res.desc$Dim.1)
#
#ind <- get_pca_ind(G3C1)
#df <- ind$coord
#
#IBC_grass$PC1score <- df[,1]
#IBC_grass$PC2score <- df[,2]
#IBC_grass$PC3score <- df[,3]
#fileName = paste0(tmp_dir,"/PCA/Grass3_PCAcoord.Rda")
#saveRDS(IBC_grass,file=fileName)

#### Forest1 (PPA)
#source(paste0(scripts_dir, "/to_test/readPPA.R"))
#
#PPA <- models$Forest1 %>%
#  ungroup() %>% # There shouldn't be groups anyways
#  mutate(Biomass = scales::rescale(Biomass,
#                                   to = c(0, 100)))%>%
#  select(-Model, -SeedRain)
#
#F1C1 <- fx_PCA_all_models("Forest1",PPA)
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
#F2C1 <- fx_PCA_all_models("Forest2",troll)
#
### Forest2 (TROLL)
source(paste0(scripts_dir, "/to_test/readTROLL.R"))

troll <- models$Forest2 %>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))%>%
    select(-Model, -SeedRain) %>%
    mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
    select(-hmax, -ah, -dmax)

F2PCA <- fx_PCA_all_models("Forest2_hrm",troll)
#pcr_model <- pcr(Biomass~h_realmax+wsg+lma+pmass+nmass, data=troll, scale =TRUE, validation="CV")
#summary(pcr_model)
fviz_eig(F2PCA, addlabels = TRUE, ylim = c(0, 50))
ggsave(paste0(tmp_dir,"/PCA/ScreePlot_Forest2.pdf")) #Eigenvals organized from largest to smallest
while (!is.null(dev.list()))  dev.off()
var <- get_pca_var(F2PCA)
corrplot(var$contrib, is.corr=FALSE)
recordPlot()

F2PCR <- pcr(Productivity~h_realmax+wsg+lma+pmass+nmass, data=troll, scale =TRUE, validation="CV")
#summary(pcr_model)
ind <- get_pca_ind(F2PCA)
df <- ind$coord

troll$PC1score <- df[,1]
troll$PC2score <- df[,2]
troll$PC3score <- df[,3]
fileName = paste0(tmp_dir,"/PCA/Forest2_PCAcoord.Rda")
saveRDS(troll,file=fileName)

#### Dryland (Bjoern)
#source(paste0(scripts_dir, "/to_test/readBjoern.R"))
#
#bjoern <- models$bjoern %>% 
#    ungroup() %>% # There shouldn't be groups anyways
#    mutate(Biomass = scales::rescale(Biomass,
#                                     to = c(0, 100)))%>%
#    select(-Model, -SeedRain, -pRoot) 
#
#DC1 <- fx_PCA_all_models("Dryland",bjoern)
