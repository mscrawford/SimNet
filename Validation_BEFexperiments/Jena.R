library(tidyverse)
library(stringr)
library(dplyr)

set.seed(1987)
base_dir          <- setwd("../")
val_dir          <- paste0(base_dir, "/Validation_BEFexperiments/")
scripts_dir       <- paste0(base_dir, "/R/to_test/")
tmp_dir           <- paste0(val_dir, "tmp/")
cache_dir           <- paste0(tmp_dir, "cache/")
raw_data_dir      <- paste0(val_dir, "data/")

source(paste0(val_dir, "functions.R"))

JenaDF<-read.csv(paste0(raw_data_dir,"JenaTBE-BarryWeigelt.et.al.csv"), sep=",")
JenaDF <- JenaDF %>%
    mutate_at("Species", str_replace, "_"," ")

# Select only the aboveground biomass values (or should we combine them with the belowground ones?)
JenaDF <- JenaDF[JenaDF$Above.Below=="Aboveground",] #%>% 
#print(summary(JenaDF))
#print(head(JenaDF))
#print(unique(JenaDF$Species))

dataset_502_18 <- read.csv(paste0(raw_data_dir,"502_18_data.csv"), sep=",") %>%
    unite(plotcode, c(Block,Plot_code), sep = "", remove = FALSE)

speciesNames <- dataset_502_18 %>%
    select(plotcode, full_name)
#print(unique(speciesNames)) #54 Spp.

monoculture_bm <- read.csv(paste0(raw_data_dir,"176_4_data.csv"), sep=",")
mono_bm <- merge(monoculture_bm,speciesNames, by.x="plotcode", by.y="plotcode")
#print(unique(mono_bm$full_name)) #27 Spp.

species.traits <- read.csv(paste0(raw_data_dir,"traits.st.csv"), sep=",") %>%
    mutate(sp = recode(sp, "Cen.jac" = "Centaurea jacea", "Pla.lan" = "Plantago lanceolata", "Poa.pra" = "Poa pratensis", "Ave.pub" = "Avenula pubescens", "Phl.pra" = "Phleum pratensis", "Dac.glo" = "Dactylis glomerata", "Ran.acr" = "Ranunculus acris", "Ant.odo" = "Anthoxanthum odoratum", "Leu.vul" = "Leucanthemum vulgare", "Ger.pra" = "Geranium pratensis", "Hol.lan" = "Holcus lanatus", "Kna.arv" = "Knautia arvensis", "Fes.rub" = "Festuca rubra", "Rum.ace" = "Rumex acetosa","Tri.fra" = "Trifolium fragiferum", "Ono.vic" = "Onobrychis viciifolia", "Tar.off" = "Taraxacum officinale", "Dau.car" = "Daucus carota", "Ran.rep" = "Ranunculus repens", "Ver.cha" = "Veronica chamaedrys", "Tri.pra" = "Trifolium pratense", "San.off" = "Sanguisorba officinalis", "Leo.his" = "Leontodon hispidus", "Med.lup" = "Medicago lupulina", "Gle.hed" = "Glechoma hederacea", "Bel.per" = "Bellis perennis", "Alo.pra" = "Alopecurus pratensis", "Med.var" = "Medicago x varia", "Pla.med" = "Plantago media", "Cre.bie" = "Crepis biennis", "Arr.ela" = "Arrhenatherum elatius", "Gal.mol" = "Galium mollugo"))# %>%
#    select(sp, leafN, RNU, cond, RootingDepth_Target)
#select(sp, la, sd, k.rad, rtd, rbm, LCN, leafN, leafC, RootingDepth_Target)
# "Ach.mil" "Aju.rep" "" "Ant.syl" ""
# "Bro.ere" "Bro.hor" "Cam.pat" "Car.car" "Car.pra" ""
# "Cir.ole" "Cyn.cri" "" "Fes.pra" ""
# "" "Her.sph" "" "" "Lat.pra"
# "Leo.aut"  "" "Lot.cor" "Luz.cam" 
# "Pas.sat" "" "Pim.maj" ""  ""
# "Poa.tri" "Pru.vul" "" 
# "Tra.pra" "Tri.cam" "Tri.dub" "Tri.fla" "Tri.hyb"
# "Tri.rep"  "Vic.cra" "Helictotrichon pubescens"
#print("#### sp")
#print(names(species.traits))
#print(unique(species.traits$sp)) #59 Spp.

print("##############################    Filtered species-traits    ##########################")
JenaDF<-JenaDF%>%
  select(Plot, Year, Species, sown_species, Biomass.g.m2)
#print(unique(JenaDF$Species)) #13 Spp.

#JenaDF_mean<-JenaDF %>%
#  filter(Species %in% species.traits$sp) %>%
#  group_by(Year, Species, sown_species, Biomass.g.m2) %>%
#  summarise(sp_biomass=mean(Biomass..g.m2.)) %>%
#  ungroup()

#print("##############################    traits    ##########################")
### Monocultures ###
df_mono <- merge(mono_bm,species.traits, by.x="full_name", by.y="sp")
df_mono <- na.omit(df_mono) %>%
    mutate(Species = full_name, species_biomass_m2 = BM_Target, NumSp = 1, Year = year, Plot = plotcode) %>%
    select(-c(full_name, BM_Target, plotcode, year, season, BM_Weed, BM_Rest, BM_Dead, Cover_Bareground))
#print("df_mono + traits")
#print(head(df_mono))
#print(unique(df_mono$full_name))

#Remove NA
JenaDF <- na.omit(JenaDF)

log_traits <- species.traits %>%
    mutate(log_leafN = log(leafN))
#There was 1 warning in `mutate()`.
#ℹ In argument: `log_leafN = log(leafN)`.
#Caused by warning in `log()`:
#! NaNs produced

JenaDF <- merge(JenaDF,species.traits, by.x="Species", by.y="sp")
#print(unique(JenaDF$Species)) #13 Spp.
#print(head(JenaDF))

JenaDF<-JenaDF[JenaDF$sown_species>0,]
JenaDF$Biomass.g.m2<-as.numeric(JenaDF$Biomass.g.m2)

#JenaDF$NumSp <- JenaDF$sown.species.richness
JenaDF$NumSp <- JenaDF$sown_species
#print(head(JenaDF))
JenaDF$species_biomass_m2 <- (JenaDF$Biomass.g.m2 * JenaDF$NumSp) # weighted biomass

JenaDF<-JenaDF %>%
    select(-sown_species, -Biomass.g.m2)#%>%
#  #group_by(Plot, Species, sown.species.richness)%>%
#  group_by(Year, Plot, Species, NumSp)%>%
#  summarise(species_biomass_m2=mean(Biomass.g.m2))#%>%
#  ungroup()

# Add an 'id' column to facilitate cforest analysis
JenaDF$id <- seq_along(JenaDF[,1])
df_mono$id <- seq_along(df_mono[,1])

print("##############################    Scatter    ##########################")

fx_plot_trait_Vs_biomass(df_mono,"Monoculture","ScatterJenaFons_justmono.png")
prin()
fx_plot_trait_Vs_biomass(JenaDF,"Monoculture","ScatterJenaFons_mono.png")
fx_plot_trait_Vs_biomass(JenaDF,"Mixture","ScatterJenaFons_mix.png")

###########################################
############### cforest ###################
###########################################

JenaDF <- na.omit(JenaDF)
df_mono <- na.omit(df_mono)
print("##############################   cforest Monoculture    ##########################")
cforest_mono <- fx_cforest_data_sets(JenaDF,"Monoculture")
write.csv(cforest_mono, paste0(cache_dir,"cforest_Jena_Fons_mono.csv"), row.names=FALSE)

print("##############################   cforest Mix    ##########################")
cforest_mix <- fx_cforest_data_sets(JenaDF,"Mixture")
write.csv(cforest_mix, paste0(cache_dir,"cforest_Jena_Fons_mix.csv"), row.names=FALSE)

print("##############################   cforest Monoculture    ##########################")
cforest_mono <- fx_cforest_data_sets(df_mono,"Monoculture")
write.csv(cforest_mono, paste0(cache_dir,"cforest_Jena_Fons_justmono.csv"), row.names=FALSE)

#print("##############################   cforest Monoculture    ##########################")
#cforest_mono <- fx_cforest_data_sets(bigbio.mono,"Monoculture")
#write.csv(cforest_mono, paste0(cache_dir,"cforest_Jena_Fons_mono_2t.csv"), row.names=FALSE)
#
#print("##############################   cforest Mix    ##########################")
#cforest_mix <- fx_cforest_data_sets(bigbio.mix,"Mixture")
#write.csv(cforest_mix, paste0(cache_dir,"cforest_Jena_Fons_mix_2t.csv"), row.names=FALSE)


#%%%%% To do %%%%%
#Check code where the biomass comes from and to see why we only have 13 spp.


#S1.2. Trait measurements
#Supplementary Table 2: Overview of traits
#Trait ; Unit ; Description
#shoot:root ratio ; g g-1 ; Shoot mass per root mass
#shoot:root N ratio ; unitless ; Leaf nitrogen uptake / root nitrogen uptake
#plant height ; cm ; Standing height of the shoot
#leaf biomass production rate ; g day-1 ; Maximum daily leaf dry mass production
#total leaf area ; cm2 ; Total area of all leaves of plant
#leaf area ; mm2 ; Average area of a single leaf
#leaf thickness ; mm ; Leaf thickness
#specific leaf area ; mm2 g-1 ; Fresh leaf area per leaf dry mass
#leaf specific density ; g cm-3 ; Leaf dry weight per leaf fresh volume
#leaf area ratio ; cm2 g-1 ; Leaf area per shoot mass
#leaf form coefficient ; mm2 mm ; Leaf area divided by leaf perimeter
#leaf dry matter content ; g g-1 ; Leaf dry weight per leaf fresh weight
#leaf C content ; % ; Leaf carbon content
#leaf N content ; % ; Leaf nitrogen Content
#leaf conductance ; μM s-1 A-1 ; Stomatal conductance per leaf area
#leaf toughness ; N ; Leaf resistance to penetration
#stem diameter ; mm ; Diameter of stem
#stem specific density ; g cm-3 ; Stem dry weight per stem fresh volume
#erectness ; cm cm-1 ; Stretched height per standing height
#biomass fraction inflorescence ; mg mg-1 ; Inflorescence:shoot biomass fraction
#inflorescences per shoot ; nr ; Number of inflorescences per shoot
#duration flowering ; ordinal ; Duration of flowering period
#seeds projected area ; mm2 ; Total area of individual seed
#nr seedlings ; nr ; Number of plant seedlings within subplot
#seed weight ; g ; Weight of 1000 seeds
#seed width length ratio ; mm mm-1 ; Ratio of seed width to seed length
#seed dry matter content ; g g-1 ; Seed dry weight per seed fresh weight
#root area ; cm2 ; Root area
#rooting depth ; ordinal ; Depth of the root system
#root area distribution ; unitless ; Evenness of vertical root area distribution
#specific root area ; cm2 g-1 ; Root surface area per root mass
#specific root length ; cm g-1 ; Root length per root mass
#root tissue density ; g cm-3 ; Root dry weight per root volume
#root nitrogen uptake ; mg day-1 ; Nitrogen uptake into roots
#root CN ratio ; unitless ; Root total carbon:nitrogen content
#root P content ; ‰ ; P content per root dry biomass
#root K content ; ‰ ; K content per root dry biomass
#root S content ; ‰ ; S content per root dry biomass
#root Ca content ; ‰ ; Ca content per root dry biomass
#root Na content ; ‰ ; Na content per root dry biomass
#nutrient uptake efficiency ; mg g-1 ; Root nitrogen uptake:root biomass
