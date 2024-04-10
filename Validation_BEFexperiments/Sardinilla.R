library(dplyr)

set.seed(1987)
base_dir          <- setwd("../")
val_dir          <- paste0(base_dir, "/Validation_BEFexperiments/")
scripts_dir       <- paste0(base_dir, "/R/to_test/")
tmp_dir           <- paste0(val_dir, "tmp/")
cache_dir           <- paste0(tmp_dir, "cache/")
raw_data_dir      <- paste0(val_dir, "data/")

source(paste0(val_dir, "functions.R"))

data <- read.csv(paste0(raw_data_dir,"DB_Inv_Sardinilla_FS.csv"), sep=",") %>%
    mutate(Species = recode(Species, "LS"="Luehea seemannii", "CA"="Cordia alliodora", "AE"="Anacardium excelsum", "HC"="Hura crepitans", "TR"="Tabebuia rosea", "CM"="Cedrela odorata")) %>%
    mutate(NumSp = richness_s_P)
print(names(data))
#"Plot"         "Year"         "TreeID"       "richness_s_P" "Species"
#"Subplot"      "FID"          "Status"       "Harvested"    "PlotEdge"
#"ba_tree"      "Height"       "plot_size"    "X_UTM17N"     "Y_UTM17N"
#"NumSp"

db <- data[data$PlotEdge==FALSE,]
db <- db[db$Harvested==FALSE,]
db <- db[db$Status=="Alive",] %>%
    group_by(Species,Year,NumSp,Plot) %>%
    summarize(BA_sp = sum(ba_tree*(10000/plot_size)/10000, na.rm = TRUE))

print(head(db))
#db[, BA_sp := sum(ba_tree*(10000/plot_size)/10000, na.rm = TRUE), .(Plot, Year, Species)]

#print(dim(data))

td <- read.delim(paste0(raw_data_dir,"Supplement_20100505.txt"))
#GENUS$ - genus (truncated at 12 letters)
#SPECIES$ - species (truncated at 12 letters)
#FAMILY$ - family (truncated at 12 letters)
#WSG  Mean wood specific gravity (g cm-3)
#SEEDMASS  Mean seed dry mass (g)
#HEIGHT  Mean height of six largest individuals (m)
#LMA  Mean leaf mass per area determined from leaf discs for the six smallest individuals of each species (g m-2)
#RGR95SAP  95th percentile relative growth rate for saplings (cm cm-1 yr-1)
#RGR90SAP  90th percentile relative growth rate for saplings (cm cm-1 yr-1)
#RGRAVGSAP  Mean relative growth rate for saplings (cm cm-1 yr-1)
#N_RGRSAP  sample size to determine relative growth rates for saplings
#MRT25SAP  Mortality rate for the 25% of saplings with the slowest growth rates in the previous census (% 5 yr-1)
#MRT50SAP  Mortality rate for the 50% of saplings with the slowest growth rates in the previous census (% 5 yr-1)
#MRTALLSAP  Mortality rate for all saplings (% 5 yr-1)
#N_MRTSAP  sample size to determine mortality rates for saplings
#RGR95TRE  95th percentile relative growth rate for large trees (cm cm-1 yr-1)
#RGR90TRE  90th percentile relative growth rate for large trees (cm cm-1 yr-1)
#RGRAVGTRE  Mean relative growth rate for large trees (cm cm-1 yr-1)
#N_RGRTRE  sample size to determine relative growth rates for large trees
#MRT25TRE  Mortality rate for the 25% of large trees with the slowest growth rates in the previous census (% 5 yr-1)
#MRT50TRE  Mortality rate for the 50% of large trees with the slowest growth rates in the previous census (% 5 yr-1)
#MRTALLTRE  Mortality rate for all large trees (% 5 yr-1)
#N_MRTTRE  sample size to determine mortality rates for large trees
td$Species <- paste(td$GENUS.,td$SPECIES.)

td_PPA <- read.delim(paste0(raw_data_dir,"ppa_full_traits.txt")) %>%
    mutate(Species = recode(sp, "LUEHSE"="Luehea seemannii", "CORDAL"="Cordia alliodora", "ANACEX"="Anacardium excelsum", "HURACR"="Hura crepitans", "TAB1RO"="Tabebuia rosea", "CEDROD"="Cedrela odorata"))

##print(dim(td))
##print(unique(td_PPA$Species))
#
##print(unique(data$Species))
##print(unique(td$Species))
#df <- data[data$PlotEdge==FALSE,]
#print(dim(df))
#df <- df[df$Harvested==FALSE,]
#print(dim(df))
#df <- df[df$Status=="Alive",]
#print(dim(df))
#print(unique(df$Species))
#print(head(df))
#print(names(df))
#
#df_1 <- merge(td, df, by.x="Species", by.y="Species")
##print(head(df_1))
#df_1[df_1=="-99"]<-NA
##print(head(df_1))
##print(names(df_1))
##"Species"      "GENUS."       "SPECIES."     "FAMILY."      "WSG"
##"SEEDMASS"     "HEIGHT"       "LMA"          "RGR95SAP"     "RGR90SAP"
##"RGRAVGSAP"    "N_RGRSAP"     "MRT25SAP"     "MRT50SAP"     "MRTALLSAP"
##"N_MRTALLSAP"  "RGR95TRE"     "RGR90TRE"     "RGRAVGTRE"    "N_RGRTRE"
##"MRT25TRE"     "MRT50TRE"     "MRTALLTRE"    "N_MRTALLTRE"  "Plot"
##"Year"         "TreeID"       "richness_s_P" "Subplot"      "FID"
##"Status"       "Harvested"    "PlotEdge"     "ba_tree"      "Height"
##"plot_size"    "X_UTM17N"     "Y_UTM17N"     "NumSp"
#
#sardinilla_1 <- na.omit(df_1)%>% 
#    select(-c(GENUS.,SPECIES.,FAMILY.,TreeID,richness_s_P,Subplot,FID,Status,Harvested,PlotEdge,plot_size,X_UTM17N,Y_UTM17N)) %>%
#    #select(Species,WSG,HEIGHT,LMA,Year,ba_tree,Plot,NumSp) %>%
#    #select(Species,WSG,SEEDMASS,HEIGHT,Height,LMA,Year,ba_tree,NumSp) %>%
#    mutate(species_biomass_m2 = ba_tree * NumSp) %>%
#    #filter(Year > 2012) %>%
#    select(-ba_tree)
##print(head(sardinilla_1))
##print("######################### summary")
##print(unique(df_1$Species))
##print(unique(sardinilla_1$Species))
##### NOTE: after na.omit, only "Tabebuia rosea" is left.
#
#sardinilla_1$id <- seq_along(sardinilla_1[,1])
#
#df_PPA <- merge(td_PPA, df, by.x="Species", by.y="Species")
##print(names(df_PPA))
##"Species"         "sp"              "WSG_AVG"         "WSG_SD"
##"MAXHEIGHT_AVG"   "MAXHEIGHT_SD"    "SEED__DRY_MIX"   "LMADISC_AVCOMB"
##"LMADISC_SDCOMB"  "LEAFAREA_AVCOMB" "LEAFAREA_SDCOMB" "LDMC_SDCOMB"
##"LDMC_AVI"        "LDMC_AVCOMB"     "LFP_COMB"        "LFN_COMB"
##"Plot"            "Year"            "TreeID"          "richness_s_P"
##"Subplot"         "FID"             "Status"          "Harvested"
##"PlotEdge"        "ba_tree"         "Height"          "plot_size"
##"X_UTM17N"        "Y_UTM17N"        "NumSp"
#
#sardinilla <- na.omit(df_PPA)%>% 
#    select(-c(sp,WSG_SD,MAXHEIGHT_SD,LMADISC_SDCOMB,LEAFAREA_SDCOMB,LDMC_SDCOMB,TreeID,richness_s_P,Subplot,FID,Status,Harvested,PlotEdge,Height,plot_size,X_UTM17N,Y_UTM17N)) %>%
#    mutate(species_biomass_m2 = ba_tree * NumSp) %>%
#    #filter(Year > 2012) %>%
#    select(-ba_tree)
#print(head(sardinilla))
#
#sardinilla$id <- seq_along(sardinilla[,1])
#

df_PPA <- merge(td_PPA, db, by.x="Species", by.y="Species")
print(names(df_PPA))
#"Species"         "sp"              "WSG_AVG"         "WSG_SD"
#"MAXHEIGHT_AVG"   "MAXHEIGHT_SD"    "SEED__DRY_MIX"   "LMADISC_AVCOMB"
#"LMADISC_SDCOMB"  "LEAFAREA_AVCOMB" "LEAFAREA_SDCOMB" "LDMC_SDCOMB"
#"LDMC_AVI"        "LDMC_AVCOMB"     "LFP_COMB"        "LFN_COMB"
#"Year"            "NumSp"           "Plot"            "BA_sp"

sardinilla <- na.omit(df_PPA)%>% 
    select(-c(sp,WSG_SD,MAXHEIGHT_SD,LMADISC_SDCOMB,LEAFAREA_SDCOMB,LDMC_SDCOMB)) %>%
    mutate(species_biomass_m2 = BA_sp * NumSp) %>%
    #filter(Year > 2012) %>%
    select(-BA_sp)
print(head(sardinilla))

sardinilla$id <- seq_along(sardinilla[,1])

print("######################### PPA")
print(summary(sardinilla))
#print(unique(df_PPA$Species))
#print(unique(sardinilla$Species))
#print(unique(sardinilla$Year))
#### NOTE: after na.omit, all six species remain.

sardinilla2 <- sardinilla %>%
    select(-c(LFP_COMB,SEED__DRY_MIX))

sardinilla2_2012 <- sardinilla %>%
    select(-c(LFP_COMB,SEED__DRY_MIX)) %>%
    filter(Year > 2012)


print("##############################    Scatter    ##########################")
#sar_1.mono <- sardinilla_1[sardinilla_1$NumSp==1,]
#sar_1.mix <- sardinilla_1[sardinilla_1$NumSp>1,]
#fx_plot_trait_Vs_biomass(sar_1.mono, "ScatterSardinilla_1_mono.png")
#fx_plot_trait_Vs_biomass(sar_1.mix, "ScatterSardinilla_1_mix.png")

sar.mono <- sardinilla[sardinilla$NumSp==1,]
sar.mix <- sardinilla[sardinilla$NumSp>1,]
fx_plot_trait_Vs_biomass(sar.mono, "ScatterSardinilla_mono_.png")
fx_plot_trait_Vs_biomass(sar.mix, "ScatterSardinilla_mix_.png")

sar2_2012.mono <- sardinilla2_2012[sardinilla2_2012$NumSp==1,]
sar2_2012.mix <- sardinilla2_2012[sardinilla2_2012$NumSp>1,]
fx_plot_trait_Vs_biomass(sar2_2012.mono, "ScatterSardinilla2_2012_mono_.png")
fx_plot_trait_Vs_biomass(sar2_2012.mix, "ScatterSardinilla2_2012_mix_.png")

############################################
#print("################ cforest ###################")
############################################
#
#cforest_mono_1 <- fx_cforest_data_sets(sardinilla_1,"Monoculture")
#write.csv(cforest_mono_1, paste0(cache_dir,"cforest_sardinilla_1_mono.csv"), row.names=FALSE)
#
#cforest_mix_1 <- fx_cforest_data_sets(sardinilla_1,"Mixture")
#write.csv(cforest_mix_1, paste0(cache_dir,"cforest_sardinilla_1_mix.csv"), row.names=FALSE)
#
cforest_mono <- fx_cforest_data_sets(sardinilla,"Monoculture")
write.csv(cforest_mono, paste0(cache_dir,"cforest_sardinilla_mono_.csv"), row.names=FALSE)

cforest_mix <- fx_cforest_data_sets(sardinilla,"Mixture")
write.csv(cforest_mix, paste0(cache_dir,"cforest_sardinilla_mix_.csv"), row.names=FALSE)

cforest_mono <- fx_cforest_data_sets(sardinilla2,"Monoculture")
write.csv(cforest_mono, paste0(cache_dir,"cforest_sardinilla2_mono_.csv"), row.names=FALSE)

cforest_mix <- fx_cforest_data_sets(sardinilla2,"Mixture")
write.csv(cforest_mix, paste0(cache_dir,"cforest_sardinilla2_mix_.csv"), row.names=FALSE)

cforest_mono <- fx_cforest_data_sets(sardinilla2_2012,"Monoculture")
write.csv(cforest_mono, paste0(cache_dir,"cforest_sardinilla2_2012_mono_.csv"), row.names=FALSE)

cforest_mix <- fx_cforest_data_sets(sardinilla2_2012,"Mixture")
write.csv(cforest_mix, paste0(cache_dir,"cforest_sardinilla2_2012_mix_.csv"), row.names=FALSE)
### To do
#Scatter plots with color coded years
#read reviews - think of a plan
