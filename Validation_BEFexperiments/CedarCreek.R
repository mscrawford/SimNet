#### Databases from Cedar Creek:
## Biodiversity II: Effects of Plant Biodiversity on Population and Ecosystem Processes - Experiment 120
##	Plant aboveground biomass carbon and nitrogen (Years Available: 1996 - 2006): https://www.cedarcreek.umn.edu/research/data/dataset?nbe120
##	Plant aboveground biomass data (Years Available: 1996 - 2020): https://www.cedarcreek.umn.edu/research/data/dataset?ple120
##	Plant traits (Years Available: 2008): https://www.cedarcreek.umn.edu/research/data/dataset?aafe120
##	Root biomass data (Years Available: 1997 - 2015): https://www.cedarcreek.umn.edu/research/data/dataset?rbe120
#library(EML)
##eml <- read_eml("knb-lter-cdr.270.122.xml", from="xml")
#eml <- read_eml("knb-lter-cdr.414.122.xml", from="xml")
##print(head(eml))library(tidyverse)
# :vsp ~/Documents/FacilitationInBEF-main/Code/Wright et al. - MonoMix - Data Cleaning-6.7.2022.Rmd
library(readr)
library(dplyr)

set.seed(1987)
base_dir          <- setwd("../")
val_dir          <- paste0(base_dir, "/Validation_BEFexperiments/")
scripts_dir       <- paste0(base_dir, "/R/to_test/")
tmp_dir           <- paste0(val_dir, "tmp/")
cache_dir           <- paste0(tmp_dir, "cache/")
raw_data_dir      <- paste0(val_dir, "data/")

source(paste0(val_dir, "functions.R"))
#~/Documents/FacilitationInBEF-main/Code/Wright et al. - MonoMix - Data Cleaning-6.7.2022.Rmd

td <- read.delim(paste0(raw_data_dir,"e120_Plant aboveground biomass carbon and nitrogen.txt"))
td$id <- seq_along(td[,1])
#print(head((td)))
td_mono <- td[td$NumSp == 1,] %>%
    select(id, Year, pTotalCb, pTotalNb, Achmi, Agrsm, Amocan, Andge, Asctu, Elyca, Koecr, Lesca, Liaas, Luppe, Monfi, Panvi, Petpu, Poapr, Queel, Quema, Schsc, Sornu) %>%
    pivot_longer(!c(pTotalCb, pTotalNb,id, Year), names_to = "sp", values_to = "count")
td_N <- td_mono[td_mono$count == 1,] %>%
    mutate(Species=recode(sp, "Achmi" = "Achillea millefolium", "Agrsm" = "Agropyron smithii", "Amoca" = "Amorpha canescens", "Andge" = "Andropogon gerardi", "Asctu" = "Asclepias tuberosa", "Bargr" = "Bare ground", "Elyca" = "Elymus canadensis", "Koecr" = "Koeleria cristata", "Lesca" = "Lespedeza capitata", "Liaas" = "Liatris aspera", "Luppe" = "Lupinus perennis", "Monfi" = "Monarda fistulosa", "Panvi" = "Panicum virgatum", "Petca" = "Petalostemum candidum", "Petpu" = "Petalostemum purpureum", "Petvi" = "Petalostemum villosum", "Poapr" = "Poa pratensis", "Queel" = "Quercus ellipsoidalis", "Quema" = "Quercus macrocarpa", "Schsc" = "Schizachyrium scoparium", "Solri" = "Solidago rigida", "Sornu" = "Sorghastrum nutans")) %>%
  group_by(Species)%>%
  summarise(leafN=mean(pTotalNb), leafC=mean(pTotalCb))%>%
  ungroup()

#print(td$pTotalNb == 1.74661,)
#print(unique(td_N$Species))
#print(td_N)

species<-read.csv(paste0(raw_data_dir,"CDRLTERe120PlantedSpecies.csv"), sep=",")
target_spp <- unique(species$Species)
print(target_spp)
BigBio<-read.csv(paste0(raw_data_dir,"CedarCreek-BigBio.csv"), sep=";")
#print(str(BigBio))
#print(summary(BigBio))
species.names<-read.csv(paste0(raw_data_dir,"BigBio-SpeciesNames.csv"), sep=";")
print(species.names)
species.traits<-read.csv(paste0(raw_data_dir,"CedarCreek_traits.csv"), sep=",")
#print(dim(species.traits))
species.traits<-species.traits %>%
	select(-Year)

BigBio<-BigBio[BigBio$NumSp>0,]
BigBio$Biomass..g.m2.<-as.numeric(BigBio$Biomass..g.m2.)
BigBioSpMatrix<-BigBio[,c(4,18:35)]
BigBioSpMatrix<-BigBioSpMatrix[!duplicated(BigBioSpMatrix$Plot),]

BigBio$Species<-if_else(BigBio$Species=="Achillea millefolium(lanulosa)", "Achillea millefolium", BigBio$Species)

BigBio<-BigBio%>%
  filter(Species %in% species.names$Species.names)
print("############# head")
print(head(BigBio))

BigBio<-BigBio %>%
  group_by(Year, Plot, Species, NumSp) %>%
  summarise(bm = mean(Biomass..g.m2.)) %>%
  ungroup()

print(head(BigBio))
#prin()

BigBio<-merge(BigBio, BigBioSpMatrix, by="Plot")

BigBio<-BigBio%>%
  arrange(Species)%>%
  pivot_wider(names_from=Species, values_from=bm)

BigBio<-BigBio%>%
  select(-c("Petalostemum candidum","Petalostemum villosum", "Solidago rigida"))
BigBio<-as.data.frame(BigBio)
for(i in 22:39){
  BigBio[,i]<-if_else(BigBio[,i-18]==1&is.na(BigBio[,i]),0,BigBio[,i])
  BigBio[,i]<-ifelse(BigBio[,i-18]==0 & is.na(BigBio[,i]) == FALSE,NA,BigBio[,i])
}

BigBio<-BigBio%>%
  select(-c(Achmi:Sornu))%>%
  pivot_longer(cols=c(4:21), names_to="Species", values_to="bm")

count.table<-aggregate(data = BigBio,                # Applying aggregate
                          Species ~ Plot,
                          function(Species) length(unique(Species)))

######### find plots with woody species, to remove them
#woody<-BigBio[grep("Quercus",BigBio$Species),]
#wp <- unique(woody$Plot)
#BigBio <- BigBio %>%
#    filter(!(Plot %in% wp))
#########

######### Remove oak trees
non_woody<-BigBio[-grep("Quercus ellipsoidalis|Quercus macrocarpa",BigBio$Species),]

CedarSmall_mean <- non_woody %>%
    group_by(Year, Species, NumSp) %>%
    summarise(species_biomass_m2 = mean((bm * NumSp), na.rm=TRUE)) %>%
    ungroup()
CedarSmall_mean<-na.omit(CedarSmall_mean)
print("&&&&&&&& non_woody &&&&&&&&")
print(summary(CedarSmall_mean))
print("&&&&&&&&&&&&&&&&")

CedarSmall_mean <- merge(CedarSmall_mean,species.traits, by.x="Species", by.y="Species")
CedarSmall_mean <- merge(CedarSmall_mean,td_N, by.x="Species", by.y="Species")
print(summary(CedarSmall_mean))
#print(head(CedarSmall_mean))
print(dim(CedarSmall_mean))


CedarSmall <- non_woody %>%
#    group_by(Year, Plot, Species, NumSp) %>%
#    summarise(species_biomass_m2 = mean((bm * NumSp), na.rm=TRUE)) %>%
    mutate(species_biomass_m2 = bm * NumSp) %>%
    select(-bm)# %>%
#    ungroup()
CedarSmall<-na.omit(CedarSmall)
print("&&&&&&&& non_woody &&&&&&&&")
print(summary(CedarSmall))
print("&&&&&&&&&&&&&&&&")

CedarSmall <- merge(CedarSmall,species.traits, by.x="Species", by.y="Species")
CedarSmall <- merge(CedarSmall,td_N, by.x="Species", by.y="Species")
print(summary(CedarSmall))
#print(head(CedarSmall))
print(dim(CedarSmall))

print("&&&&&&&& 2007 &&&&&&&&")
cedar2007 <- non_woody %>%
    group_by(Year, Plot, Species, NumSp) %>%
    summarise(species_biomass_m2 = mean((bm * NumSp), na.rm=TRUE)) %>%
    filter(Year == 2007) %>%
    ungroup()
cedar2007<-na.omit(cedar2007)

cedar2007 <- merge(cedar2007,species.traits, by.x="Species", by.y="Species")
cedar2007 <- merge(cedar2007,td_N, by.x="Species", by.y="Species")
print(summary(cedar2007))
#print(head(cedar2007))
print(dim(cedar2007))

cedar7y <- CedarSmall %>%
    filter(Year < 2008) %>%
    ungroup()
cedar7y<-na.omit(cedar7y)
print("    %%%%%&&&&&&&&$$$$$$ cedar < 2007")
print(summary(cedar7y))

#CedarSmall <- BigBio[BigBio$height_.m. < 3,] #Check that there are no trees in dataframe
#print(dim(CedarSmall))
# Add an ')id' column to facilitate cforest analysis
CedarSmall$id <- seq_along(CedarSmall[,1])
cedar2007$id <- seq_along(cedar2007[,1])
cedar7y$id <- seq_along(cedar7y[,1])
#CedarSmall$Plot <- seq_along(CedarSmall[,1])

#print("##############################    Merged    ##########################")

CedarSmall_log <- CedarSmall %>% 
    mutate(Log_P_A.L = log(P_A.L)) %>%
    mutate(Log_P.A = log(P.A)) %>%
    mutate(Log_Seed_weight_.g. = log(Seed_weight_.g.)) %>%
    #mutate(species_biomass_m2 = log(species_biomass_m2)) %>%
    select(-c(P_A.L,P.A,Seed_weight_.g.))
#print(str(CedarSmall))

print("&&&&&&&& head 2007 &&&&&&&&")
cedar2007log <- cedar2007 %>% 
    mutate(Log_P_A.L = log(P_A.L)) %>%
    mutate(Log_P.A = log(P.A)) %>%
    mutate(Log_Seed_weight_.g. = log(Seed_weight_.g.)) %>%
    #mutate(species_biomass_m2 = log(species_biomass_m2)) %>%
    select(-c(P_A.L,P.A,Seed_weight_.g.))

#############################################################################################
############################### Traits vs. Biomass ##########################################
#############################################################################################
print("##############################    Scatter    ##########################")
fx_plot_trait_Vs_biomass(CedarSmall, "Monoculture", "ScatterCedar_mono.png")
fx_plot_trait_Vs_biomass(CedarSmall, "Mixture", "ScatterCedar_mix.png")
fx_plot_trait_Vs_biomass(cedar2007log, "Monoculture", "ScatterCedar2007_mono.png")
fx_plot_trait_Vs_biomass(cedar2007log, "Mixture", "ScatterCedar2007_mix.png")
fx_plot_trait_Vs_biomass(cedar7y, "Monoculture", "ScatterCedar_7y_mono.png")
fx_plot_trait_Vs_biomass(cedar7y, "Mixture", "ScatterCedar_7y_mix.png")

############################################
#print("################ cforest ###################")
############################################
#
##Traits:  Area_of_leaf_blade_cm2, P.A, P_A.L, SLA_cm2.g, Seed_weight_.g., height_.m.
##Year | Definition: Year sampled
##Species | Definition: Species Name
##Area_of_leaf_blade_cm2 | Definition: Leaf laminar area (cm^2) | Unit: squareCentimeters
##P/A | Definition: Perimeter per area (cm-1) | Unit: waveNumber
##P_A*L | Definition: Perimeter per area times leaf laminar length (unitless) | Unit: dimensionless
##SLA_cm2/g | Definition: Specific leaf length (cm^2/g) | Unit: centimetersSquaredPerGram
##Seed_weight_(g) | Definition: Seed weight (g) | Unit: gram
##height_(m) | Definition: Average plant height at maturity (m) | Unit: meter
## Seed_weight_.g.
## SLA_cm2.g
## P.A
## P_A.L
## height_.m.
## Area_of_leaf_blade_cm2

cedar_mono <- fx_cforest_data_sets(CedarSmall,"Monoculture")
write.csv(cedar_mono, paste0(cache_dir,"cforest_cedar_mono.csv"), row.names=FALSE)
prin()

cedar_mix <- fx_cforest_data_sets(CedarSmall,"Mixture")
write.csv(cedar_mix, paste0(cache_dir,"cforest_cedar_mix.csv"), row.names=FALSE)

cedar_mono_7y <- fx_cforest_data_sets(cedar7y,"Monoculture")
write.csv(cedar_mono_7y, paste0(cache_dir,"cforest_cedar_mono_7y.csv"), row.names=FALSE)

cedar_mix_7y <- fx_cforest_data_sets(cedar7y,"Mixture")
write.csv(cedar_mix_7y, paste0(cache_dir,"cforest_cedar_mix_7y.csv"), row.names=FALSE)

cforest_mono <- fx_cforest_data_sets(cedar2007,"Monoculture")
write.csv(cforest_mono, paste0(cache_dir,"cforest_cedar_mono_2007.csv"), row.names=FALSE)

cforest_mix <- fx_cforest_data_sets(cedar2007,"Mixture")
write.csv(cforest_mix, paste0(cache_dir,"cforest_cedar_mix_2007.csv"), row.names=FALSE)

cedar2007_2t<-cedar2007%>%
select(id, Year, Species, NumSp, species_biomass_m2, leafN, height_.m.)
print("%%%%%%%%%%%%%%%%%%%%% data  (rf) %%%%%%%%%%%%%%%%%%%%%%%%%")
print(summary(cedar2007_2t))

cforest_2t_mono <- fx_cforest_data_sets(cedar2007_2t,"Monoculture")
write.csv(cforest_2t_mono, paste0(cache_dir,"cforest_cedar_2t_mono_2007"), row.names=FALSE)

cforest_2t_mix <- fx_cforest_data_sets(cedar2007_2t,"Mixture")
write.csv(cforest_2t_mix, paste0(cache_dir,"cforest_cedar_2t_mix_2007"), row.names=FALSE)


#cforest_mono <- fx_cforest_data_sets(CedarSmall,"Monoculture")
#write.csv(cforest_mono, paste0(cache_dir,"cforest_cedar_mono_mean.csv"), row.names=FALSE)
#
#cforest_mix <- fx_cforest_data_sets(CedarSmall,"Mixture")
#write.csv(cforest_mix, paste0(cache_dir,"cforest_cedar_mix_mean.csv"), row.names=FALSE)

CedarSmall<-CedarSmall%>%
select(id, Year, Species, NumSp, species_biomass_m2, leafN, height_.m.)
#select(id, Plot, Year, Species, NumSp, species_biomass_m2, leafN, height_.m.)

cforest_2t_mono <- fx_cforest_data_sets(CedarSmall,"Monoculture")
write.csv(cforest_2t_mono, paste0(cache_dir,"cforest_cedar_2t_mono_mean.csv"), row.names=FALSE)

cforest_2t_mix <- fx_cforest_data_sets(CedarSmall,"Mixture")
write.csv(cforest_2t_mix, paste0(cache_dir,"cforest_cedar_2t_mix_mean.csv"), row.names=FALSE)
