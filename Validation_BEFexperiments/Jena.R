library(tidyverse)
library(party)
library(stringr)

set.seed(1987)
base_dir          <- setwd("../")
val_dir          <- paste0(base_dir, "/Validation_BEFexperiments/")
scripts_dir       <- paste0(base_dir, "/R/to_test/")
tmp_dir           <- paste0(val_dir, "tmp/")
cache_dir           <- paste0(tmp_dir, "cache/")
raw_data_dir      <- paste0(val_dir, "data/")

source(paste0(scripts_dir, "fx_cforest_party.R"))

BigBio<-read.csv(paste0(raw_data_dir,"JenaTBE-BarryWeigelt.et.al.csv"), sep=",")
BigBio <- BigBio %>%
    mutate_at("Species", str_replace, "_"," ")

# Select only the aboveground biomass values (or should we combine them with the belowground ones?)
BigBio <- BigBio[BigBio$Above.Below=="Aboveground",] #%>% 
print("##############################    orig    ##########################")
#print(summary(BigBio))
#print(head(BigBio))
#print(unique(BigBio$Species))

dataset_502_18 <- read.csv(paste0(raw_data_dir,"502_18_data.csv"), sep=",") %>%
    unite(plotcode, c(Block,Plot_code), sep = "", remove = FALSE)
speciesNames <- dataset_502_18 %>%
    select(plotcode, full_name)
print(unique(speciesNames))
print("#### mono")
monoculture_bm <- read.csv(paste0(raw_data_dir,"176_4_data.csv"), sep=",")
mono_bm <- merge(monoculture_bm,speciesNames, by.x="plotcode", by.y="plotcode")
print(unique(mono_bm$full_name))
######### to do #########
# multiple regression -> check random forest results
# select species from monoculture_bm based on speciesNames

species.traits <- read.csv(paste0(raw_data_dir,"traits.st.csv"), sep=",") %>%
    mutate(sp = recode(sp, "Cen.jac" = "Centaurea jacea", "Pla.lan" = "Plantago lanceolata", "Poa.pra" = "Poa pratensis", "Ave.pub" = "Avenula pubescens", "Phl.pra" = "Phleum pratensis", "Dac.glo" = "Dactylis glomerata", "Ran.acr" = "Ranunculus acris", "Ant.odo" = "Anthoxanthum odoratum", "Leu.vul" = "Leucanthemum vulgare", "Ger.pra" = "Geranium pratensis", "Hol.lan" = "Holcus lanatus", "Kna.arv" = "Knautia arvensis", "Fes.rub" = "Festuca rubra", "Rum.ace" = "Rumex acetosa","Tri.fra" = "Trifolium fragiferum", "Ono.vic" = "Onobrychis viciifolia", "Tar.off" = "Taraxacum officinale", "Dau.car" = "Daucus carota", "Ran.rep" = "Ranunculus repens", "Ver.cha" = "Veronica chamaedrys", "Tri.pra" = "Trifolium pratense", "San.off" = "Sanguisorba officinalis", "Leo.his" = "Leontodon hispidus", "Med.lup" = "Medicago lupulina", "Gle.hed" = "Glechoma hederacea", "Bel.per" = "Bellis perennis", "Alo.pra" = "Alopecurus pratensis", "Med.var" = "Medicago x varia", "Pla.med" = "Plantago media", "Cre.bie" = "Crepis biennis", "Arr.ela" = "Arrhenatherum elatius", "Gal.mol" = "Galium mollugo")) %>%
    select(sp, leafN, RNU, cond, RootingDepth_Target)
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
print("#### sp")
print(unique(species.traits$sp))

print("##############################    Filtered species-traits    ##########################")
BigBio<-BigBio%>%
  filter(Species %in% species.traits$sp) %>%
  select(Plot, Year, Species, sown_species, Biomass.g.m2)
print(unique(BigBio$Species))

print("##############################    traits    ##########################")
### Monocultures ###
df_mono <- merge(mono_bm,species.traits, by.x="full_name", by.y="sp")
df_mono <- na.omit(df_mono) %>%
    mutate(Species = full_name, species_biomass_m2 = BM_Target ) %>%
    select(-c(full_name, BM_Target, plotcode, year, season, BM_Weed, BM_Rest, BM_Dead, Cover_Bareground))
print("df_mono + traits")
print(head(df_mono))
#print(unique(df_mono$full_name))

BigBio <- merge(BigBio,species.traits, by.x="Species", by.y="sp")
#print(str(BigBio))
BigBio<-BigBio[BigBio$sown_species>0,]
BigBio$Biomass.g.m2<-as.numeric(BigBio$Biomass.g.m2)
#print(head(BigBio))
#print(summary(BigBio))

#BigBio$NumSp <- BigBio$sown.species.richness
BigBio$NumSp <- BigBio$sown_species
#print(head(BigBio))
BigBio$species_biomass_m2 <- (BigBio$Biomass.g.m2 * BigBio$NumSp) # weighted biomass
#BigBio$Year <- BigBio$sown_species #make a Year column to avoid problems with fx_cforest_party.R

BigBio<-BigBio %>%
    select(-sown_species, -Biomass.g.m2)#%>%
#  #group_by(Plot, Species, sown.species.richness)%>%
#  group_by(Year, Plot, Species, NumSp)%>%
#  summarise(species_biomass_m2=mean(Biomass.g.m2))#%>%
#  ungroup()

# Add an 'id' column to facilitate cforest analysis
BigBio$id <- seq_along(BigBio[,1])
df_mono$id <- seq_along(df_mono[,1])
#Remove NA
BigBio <- na.omit(BigBio)

bigbio.mono <- BigBio[BigBio$NumSp==1,]
print("##############################    Monoculture    ##########################")
#print(unique(bigbio.mono$Species))
print(head(bigbio.mono))

bigbio.mix <- BigBio[BigBio$NumSp>1,]
print("##############################    Mix    ##########################")
print(summary(bigbio.mix))

print("##############################    Scatter    ##########################")
print(head(bigbio.mono))
bigbio.mono %>%
    gather(-Species, -Plot, -id, -Year, -NumSp, -species_biomass_m2, key= "var", value = "value") %>%
    ggplot(aes(x=value, y=log(species_biomass_m2))) +
    geom_point() +
    labs(y = "log biomass") +
	stat_smooth(fullrange = TRUE, color="red", method="loess", se=FALSE) +
    facet_wrap(~var, scales ="free") +
    theme_bw() +
    theme(strip.text = element_text(size=9))

ggsave(file=paste0(tmp_dir,"ScatterJenaFons_mono.png")
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
    theme(strip.text = element_text(size=9))

ggsave(file=paste0(tmp_dir,"ScatterJenaFons_mix.png")
       , width=18, height=14, dpi=300
)
while (!is.null(dev.list()))  dev.off()

###########################################
############### cforest ###################
###########################################

#print("##############################   cforest Monoculture    ##########################")
#cforest_mono <- fx_cforest_data_sets(bigbio.mono,"Monoculture")
#write.csv(cforest_mono, paste0(cache_dir,"cforest_Jena_Fons_mono.csv"), row.names=FALSE)
#
#print("##############################   cforest Mix    ##########################")
#cforest_mix <- fx_cforest_data_sets(bigbio.mix,"Mixture")
#write.csv(cforest_mix, paste0(cache_dir,"cforest_Jena_Fons_mix.csv"), row.names=FALSE)
#
print("##############################   cforest Monoculture    ##########################")
cforest_mono <- fx_cforest_data_sets(df_mono,"Monoculture")
write.csv(cforest_mono, paste0(cache_dir,"cforest_Jena_Fons_mono_2t.csv"), row.names=FALSE)

#print("##############################   cforest Monoculture    ##########################")
#cforest_mono <- fx_cforest_data_sets(bigbio.mono,"Monoculture")
#write.csv(cforest_mono, paste0(cache_dir,"cforest_Jena_Fons_mono_2t.csv"), row.names=FALSE)
#
#print("##############################   cforest Mix    ##########################")
#cforest_mix <- fx_cforest_data_sets(bigbio.mix,"Mixture")
#write.csv(cforest_mix, paste0(cache_dir,"cforest_Jena_Fons_mix_2t.csv"), row.names=FALSE)
