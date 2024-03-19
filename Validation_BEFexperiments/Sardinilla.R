library(tidyverse)
library(party)
library(dplyr)

set.seed(1987)
base_dir          <- setwd("../")
val_dir          <- paste0(base_dir, "/Validation_BEFexperiments/")
scripts_dir       <- paste0(base_dir, "/R/to_test/")
tmp_dir           <- paste0(val_dir, "tmp/")
cache_dir           <- paste0(tmp_dir, "cache/")
raw_data_dir      <- paste0(val_dir, "data/")

source(paste0(val_dir, "functions.R"))
source(paste0(scripts_dir, "fx_cforest_party.R"))
data <- read.csv(paste0(raw_data_dir,"DB_Inv_Sardinilla_FS.csv"), sep=",") %>%
    mutate(Species = recode(Species, "LS"="Luehea seemannii", "CA"="Cordia alliodora", "AE"="Anacardium excelsum", "HC"="Hura crepitans", "TR"="Tabebuia rosea", "CM"="Cedrela odorata")) %>%
    mutate(NumSp = richness_s_P)
print(dim(data))
print(head(data))

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

print(dim(td))
print(head(td))

df <- merge(td, data, by.x="Species", by.y="Species")
print(dim(df))
#df1 <- merge(data, td, by.x="Species", by.y="Species")
#print(dim(df1))
print(head(df))

#print(unique(data$Species))
#print(unique(td$Species))
#print(unique(df$Species))
#print(unique(df1$Species))
df_no_edge <- df[df$PlotEdge==FALSE,]
print(dim(df_no_edge))
df_no_edge <- df_no_edge[df_no_edge$Harvested==FALSE,]
print(dim(df_no_edge))
df_no_edge <- df_no_edge[df_no_edge$Status=="Alive",]
print(dim(df_no_edge))

# "Species"      "GENUS."       "SPECIES."     "FAMILY."      "WSG"
# "SEEDMASS"     "HEIGHT"       "LMA"          "RGR95SAP"     "RGR90SAP"
# "RGRAVGSAP"    "N_RGRSAP"     "MRT25SAP"     "MRT50SAP"     "MRTALLSAP"
# "N_MRTALLSAP"  "RGR95TRE"     "RGR90TRE"     "RGRAVGTRE"    "N_RGRTRE"
# "MRT25TRE"     "MRT50TRE"     "MRTALLTRE"    "N_MRTALLTRE"  "Plot"
# "Year"         "TreeID"       "richness_s_P" "Subplot"      "FID"
# "Status"       "Harvested"    "PlotEdge"     "ba_tree"      "Height"
# "plot_size"    "X_UTM17N"     "Y_UTM17N"     "NumSp"

sardinilla <- na.omit(df_no_edge)%>% 
    select(Species,WSG,HEIGHT,LMA,Year,ba_tree,Plot,NumSp) %>%
    #select(Species,WSG,SEEDMASS,HEIGHT,Height,LMA,Year,ba_tree,NumSp)
    mutate(species_biomass_m2 = ba_tree * NumSp) %>%
    filter(Year > 2012) %>%
    select(-ba_tree)
print(head(sardinilla))

sardinilla$id <- seq_along(sardinilla[,1])

print(unique(sardinilla$Year))
print("##############################    Scatter    ##########################")
sar.mono <- sardinilla[sardinilla$NumSp==1,]
sar.mix <- sardinilla[sardinilla$NumSp>1,]
fx_plot_trait_Vs_biomass(sar.mono, "ScatterSardinilla_mono.png")
fx_plot_trait_Vs_biomass(sar.mix, "ScatterSardinilla_mix.png")
prin()

############################################
#print("################ cforest ###################")
############################################
#
cforest_mono <- fx_cforest_data_sets(sardinilla,"Monoculture")
write.csv(cforest_mono, paste0(cache_dir,"cforest_sardinilla2_mono.csv"), row.names=FALSE)

cforest_mix <- fx_cforest_data_sets(sardinilla,"Mixture")
write.csv(cforest_mix, paste0(cache_dir,"cforest_sardinilla2_mix.csv"), row.names=FALSE)


### To do
Scatter plots with color coded years
read reviews - think of a plan
