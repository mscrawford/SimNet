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
#
library(party)
library(png)
library(readr)
library(tidyr)

set.seed(1987)
base_dir          <- setwd("../")
val_dir          <- paste0(base_dir, "/Validation_BEFexperiments/")
scripts_dir       <- paste0(base_dir, "/R/to_test/")
tmp_dir           <- paste0(val_dir, "tmp/")
cache_dir           <- paste0(tmp_dir, "cache/")
raw_data_dir      <- paste0(val_dir, "data/")

source(paste0(scripts_dir, "fx_traits_vs_biomass.R"))
source(paste0(scripts_dir, "fx_cforest_party.R"))
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

BigBio<-read.csv(paste0(raw_data_dir,"CedarCreek-BigBio.csv"), sep=";")
#print(dim(BigBio))
#print(str(BigBio))
#print(summary(BigBio))
species.names<-read.csv(paste0(raw_data_dir,"BigBio-SpeciesNames.csv"), sep=";")
species.traits<-read.csv(paste0(raw_data_dir,"CedarCreek_traits.csv"), sep=",")
#print(species.traits)
species.traits<-species.traits %>%
	select(-Year)
BigBio<-BigBio[BigBio$NumSp>0,]
BigBio$Biomass..g.m2.<-as.numeric(BigBio$Biomass..g.m2.)
BigBioSpMatrix<-BigBio[,c(4,18:35)]
BigBioSpMatrix<-BigBioSpMatrix[!duplicated(BigBioSpMatrix$Plot),]
BigBio<-BigBio %>%
  group_by(Year, Plot, Species, NumSp)%>%
  summarise(species_biomass_m2=mean(Biomass..g.m2.))%>%
  ungroup()

BigBio<-merge(BigBio, BigBioSpMatrix, by="Plot")
#print(head(BigBio))
#print(summary(BigBio))

BigBio$Species<-if_else(BigBio$Species=="Achillea millefolium(lanulosa)", "Achillea millefolium", BigBio$Species)

BigBio<-BigBio%>%
  filter(Species %in% species.names$Species.names) %>%
  mutate(bm = species_biomass_m2) %>%
  mutate(species_biomass_m2 = bm * NumSp) %>%
  select(Plot, Year, Species, NumSp, species_biomass_m2)
print("##############################    orig    ##########################")
#print(str(BigBio))

BigBio <- merge(BigBio,species.traits, by.x="Species", by.y="Species")
BigBio <- merge(BigBio,td_N, by.x="Species", by.y="Species")
#print(head(BigBio))
#### remove trees from data set
CedarSmall <- BigBio[BigBio$height_.m. < 3,]
print(head(CedarSmall))

# Add an 'id' column to facilitate cforest analysis
CedarSmall$id <- seq_along(CedarSmall[,1])

#print("##############################    Merged    ##########################")
#print(str(CedarSmall))

#bigbio.mono<-CedarSmall[CedarSmall$Plot==151,] #%>% 
bigbio.mono<-CedarSmall[CedarSmall$NumSp==1,] #%>% 

bigbio.mix<-CedarSmall[CedarSmall$NumSp>1,]
#bigbio.mix<-CedarSmall[CedarSmall$NumSp==16,]

#print("##############################    Monoculture    ##########################")
#print(bigbio.mono)
#print(str(bigbio.mono))
#print("##############################    Mix    ##########################")
#print(str(bigbio.mix))

#############################################################################################
############################### Traits vs. Biomass ##########################################
#############################################################################################
print("##############################    Scatter    ##########################")
bigbio.mono %>%
    gather(-Species, -Plot, -id, -Year, -NumSp, -species_biomass_m2, key= "var", value = "value") %>%
    ggplot(aes(x=value, y=log(species_biomass_m2))) +
    geom_point() +
    labs(y = "log biomass") +
	stat_smooth(fullrange = TRUE, color="red", method="loess", se=FALSE) +
    facet_wrap(~var, scales ="free") +
    theme_bw() +
    theme(strip.text = element_text(size=25))

ggsave(file=paste0(tmp_dir,"ScatterCedar_mono.png")
       , width=18, height=14, dpi=300
)
while (!is.null(dev.list()))  dev.off()

bigbio.mix %>%
    gather(-Species, -Plot, -id, -Year, -NumSp, -species_biomass_m2, key= "var", value = "value") %>%
    ggplot(aes(x=value, y=log(species_biomass_m2))) +
    geom_point() +
    labs(y = "log biomass") +
	stat_smooth(fullrange = TRUE, color="red", method="loess", se=FALSE) +
    facet_wrap(~var, scales ="free") +
    theme_bw() +
    theme(strip.text = element_text(size=25))

ggsave(file=paste0(tmp_dir,"ScatterCedar_mix.png")
       , width=18, height=14, dpi=300
)
while (!is.null(dev.list()))  dev.off()

### mean biomass
#df <- bigbio.mono %>% 
#    group_by(Species) %>%
#    summarise(Biomass_mono = mean(species_biomass_m2), leafC = mean(leafC), leafN = mean(leafN),    Area_of_leaf_blade_cm2 = mean(Area_of_leaf_blade_cm2), P.A = mean(P.A), P_A.L = mean(P_A.L), SLA_cm2.g = mean(SLA_cm2.g), Seed_weight_.g. = mean(Seed_weight_.g.), height_.m. = mean(height_.m.))
#print(head(df))
#
#Area_of_leaf_blade_cm2
#P.A
#P_A.L
#SLA_cm2.g
#Seed_weight_.g.
#height_.m.
#leafC
#leafN

#df_mix <- bigbio.mix %>% 
#    group_by(Species) %>%
#    summarise(Biomass_mix = mean(species_biomass_m2), leafC = mean(leafC), leafN = mean(leafN),    Area_of_leaf_blade_cm2 = mean(Area_of_leaf_blade_cm2), P.A = mean(P.A), P_A.L = mean(P_A.L), SLA_cm2.g = mean(SLA_cm2.g), Seed_weight_.g. = mean(Seed_weight_.g.), height_.m. = mean(height_.m.))
#print(head(df_mix))

##df <- aggregate(df$species_biomass_m2, list(df$Species), mean, simplify = TRUE) #Calculate mean biomass per specie
##print(head(df))
##setnames(df, old=c("Group.1", "x"), new=c("Species","Biomass_mono"))
#df_t_mono$Biomass_mono <- df_t_mono$species_biomass_m2
#df_t_mix$Biomass_mix <- df_t_mix$species_biomass_m2
#png(paste0(dir, "Scatter_Cedar_meanBM.png"), 
#    #width = 1663, height = 3061, units = "px", pointsize = 40,  res = NA,
#    width = 3326, height = 3061, units = "px", pointsize = 40,  res = NA,
#    bg = "white", type = c("cairo", "cairo-png", "Xlib", "quartz"))
#par(mar=c(0,0,2,0))
#png(paste0(dir, "Scatter_Cedar_BM.png"), 
#    #width = 1663, height = 3061, units = "px", pointsize = 40,  res = NA,
#    width = 3326, height = 3061, units = "px", pointsize = 40,  res = NA,
#    bg = "white", type = c("cairo", "cairo-png", "Xlib", "quartz"))
#par(mar=c(0,0,2,0))

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
#
#cforest_mono <- fx_cforest_data_sets(CedarSmall,"Monoculture")
#write.csv(cforest_mono, paste0(cache_dir,"cforest_cedar_mono.csv"), row.names=FALSE)
#
#cforest_mix <- fx_cforest_data_sets(CedarSmall,"Mixture")
#write.csv(cforest_mix, paste0(cache_dir,"cforest_cedar_mix.csv"), row.names=FALSE)
#
CedarSmall<-CedarSmall%>%
select(id, Plot, Year, Species, NumSp, species_biomass_m2, leafN, height_.m.)

cforest_2t_mono <- fx_cforest_data_sets(CedarSmall,"Monoculture")
write.csv(cforest_2t_mono, paste0(cache_dir,"cforest_cedar_2t_mono.csv"), row.names=FALSE)

cforest_2t_mix <- fx_cforest_data_sets(CedarSmall,"Mixture")
write.csv(cforest_2t_mix, paste0(cache_dir,"cforest_cedar_2t_mix.csv"), row.names=FALSE)
