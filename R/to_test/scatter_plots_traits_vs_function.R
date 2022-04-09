base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

#source(paste0(scripts_dir, "/readModels.R"))
source(paste0(scripts_dir, "/to_test/fx_traits_vs_biomass.R"))
source(paste0(scripts_dir, "/to_test/fx_cforest_party.R"))

### Dryland (Bjoern)
bjoern <- fx_read_model("readBjoern.R","bjoern") %>%
	select(-pRoot)

lab1 <- "maxBiomass" #(maximum \n size/size at maturity) [gC]"
lab2 <- "leafAllocation" #"pLeaf (allocation \n to leaf) [gC/gC]"
lab3 <- "storageAllocation" #"pStorage (allocation \n to storage) [gC/gC]"
lab4 <- "rootsAllocation" #"pRoot (allocation \n to root) [gC/gC]"

DC1 <- fx_traits_vs_biomass_jitter("DC1",bjoern,1,meta,"maxSize","pLeaf",lab1,lab2)
DC2 <- fx_traits_vs_biomass_jitter("DC2",bjoern,1,iso,"maxSize","pLeaf",lab1,lab2)
#DC2 <- fx_traits_vs_biomass_jitter("DC2",bjoern,1,iso,"maxSize","pStorage",lab1,lab3)
DC3 <- fx_traits_vs_biomass_jitter("DC3",bjoern,32,meta,"maxSize","pLeaf",lab1,lab2)
DC4 <- fx_traits_vs_biomass_jitter("DC4",bjoern,32,iso,"maxSize","pLeaf",lab1,lab2)

DC1_P <- fx_traits_vs_biomass_jitter("DC1_P",bjoern,1,meta,"maxSize","pLeaf",lab1,lab2)
#DC1_P <- fx_traits_vs_biomass_jitter("DC1_P",bjoern,1,meta,"pLeaf","maxSize",lab2,lab1)
DC2_P <- fx_traits_vs_biomass_jitter("DC2_P",bjoern,1,iso,"maxSize","pLeaf",lab1,lab2)
DC3_P <- fx_traits_vs_biomass_jitter("DC3_P",bjoern,32,meta,"maxSize","pLeaf",lab1,lab2)
DC4_P <- fx_traits_vs_biomass_jitter("DC4_P",bjoern,32,iso,"maxSize","pLeaf",lab1,lab2)

prin()

### Grass1 (Adam's model)
adam <- fx_read_model("readAdam.R", "Grass1") %>%
	mutate(r_pNi = 1/pNi) #rpNi = reciprocal of pNi

lab1 <- "monoBiomass" #"no3i (nitrogen R*)"
lab2 <- "NUE1" #"r_pNi (reciprocal of \n aboveground N concentration)"

G1C1 <- fx_traits_vs_biomass_jitter("G1C1",adam,1,meta,"abmi","r_pNi",lab1,lab2)
G1C3 <- fx_traits_vs_biomass_jitter("G1C3",adam,32,meta,"abmi","r_pNi",lab1,lab2)
G1C2 <- fx_traits_vs_biomass_jitter("G1C2",adam,1,iso,"abmi","r_pNi",lab1,lab2)
G1C4 <- fx_traits_vs_biomass_jitter("G1C4",adam,32,iso,"abmi","r_pNi",lab1,lab2)

### Grass2 (Lindsay's model)
lindsay <- fx_read_model("readLindsay_variable.R","Grass2") 

lab1 <- "rootingVolume" #"Vi (volume of soil \n accessible to species i)"
lab2 <- "NUE2" #thetai (Nitrogen uptake \n rate per unit plant biomass)"

G2C1 <- fx_traits_vs_biomass_jitter("G2C1",lindsay,1,meta,"Vi","thetai",lab1,lab2)
G2C2 <- fx_traits_vs_biomass_jitter("G2C2",lindsay,1,iso,"Vi","thetai",lab1,lab2)
G2C3 <- fx_traits_vs_biomass_jitter("G2C3",lindsay,32,meta,"Vi","thetai",lab1,lab2)
G2C4 <- fx_traits_vs_biomass_jitter("G2C4",lindsay,32,iso,"Vi","thetai",lab1,lab2)

### Grass3 (IBC-grass)
IBC_grass <- readRDS(paste0(tmp_dir,"/PCA/Grass3_PCAcoord.Rda")) %>%
	select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, PC1score, PC2score, PC3score) %>%
	mutate(PC2score = -PC2score) %>%
	mutate(PC3score = -PC3score) %>%
	mutate(id = row_number()) %>%
	mutate_if(is.character, as.factor) %>%
	mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
	mutate(Productivity = scales::rescale(Productivity, to = c(0, 100)))

lab1 <- "LES1" #"PC1score associated with LMR \n (leaf to mass ratio)"
lab2 <- "Size/Growth" #"PC2score associated with Gmax \n (maximum resource utilization) \n and MaxMass (Plant's maximum size)"
lab3 <- "Spacing" #"PC3score associated with SLA and MeanSpacerLength"

G3C1 <- fx_traits_vs_biomass_jitter("G3C1",IBC_grass,1,meta,"PC2score","PC1score",lab2,lab1)
G3C2 <- fx_traits_vs_biomass_jitter("G3C2",IBC_grass,1,iso,"PC2score","PC1score",lab2,lab1)
G3C3 <- fx_traits_vs_biomass_jitter("G3C3",IBC_grass,32,meta,"PC2score","PC1score",lab2,lab1)
G3C4 <- fx_traits_vs_biomass_jitter("G3C4",IBC_grass,32,iso,"PC2score","PC1score",lab2,lab1)
#G3C3 <- fx_traits_vs_biomass_jitter("G3C3",IBC_grass,32,meta,"PC2score","PC3score",lab2,lab3)
#G3C4 <- fx_traits_vs_biomass_jitter("G3C4",IBC_grass,32,iso,"PC3score","PC1score",lab3,lab1)

#### Forest1 (PPA)
PPA <- fx_read_model("readPPA.R","Forest1")

lab1 <- "paceOfLife" #"PC1score (associated \n with fast-slow lifecycle)"#plant height)"
lab2 <- "MaxHeight" #"PC2score (associated \n with tree stature)"# LMA -leaf mass per area)"

F1C1 <- fx_traits_vs_biomass_jitter("F1C1",PPA,1,meta,"PC2score","PC1score",lab2,lab1)
F1C2 <- fx_traits_vs_biomass_jitter("F1C2",PPA,1,iso,"PC2score","PC1score",lab2,lab1)
F1C3 <- fx_traits_vs_biomass_jitter("F1C3",PPA,32,meta,"PC2score","PC1score",lab2,lab1)
F1C4 <- fx_traits_vs_biomass_jitter("F1C4",PPA,32,iso,"PC2score","PC1score",lab2,lab1)

F1C1_P <- fx_traits_vs_biomass_jitter("F1C1_P",PPA,1,meta,"PC2score","PC1score",lab2,lab1)
F1C2_P <- fx_traits_vs_biomass_jitter("F1C2_P",PPA,1,iso,"PC2score","PC1score",lab2,lab1)
F1C3_P <- fx_traits_vs_biomass_jitter("F1C3_P",PPA,32,meta,"PC2score","PC1score",lab2,lab1)
F1C4_P <- fx_traits_vs_biomass_jitter("F1C4_P",PPA,32,iso,"PC2score","PC1score",lab2,lab1)

### Forest2 (TROLL) h_realmax
troll <- readRDS(paste0(tmp_dir,"/PCA/Forest2_hrm_PCAcoord.Rda")) %>%
	select(c(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, PC1score, PC2score, PC3score)) %>%
	mutate(PC1score = -PC1score) %>%
	mutate(PC3score = -PC3score) %>%
	mutate(id = row_number()) %>%
	mutate_if(is.character, as.factor) %>%
	mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
	mutate(Productivity = scales::rescale(Productivity, to = c(0, 100)))

lab1 <- "LES2" #"PC1score associated with \n LMA, nmass, and pmass"
lab2 <- "MaxHeight" #"PC2score associated with \n h_realmax = hmax * dmax / (dmax + ah)"
lab3 <- "woodDensity" #"PC3score associated with \n wsg (wood specific gravity)"

#F2C1 <- fx_traits_vs_biomass_jitter("F2C1_hrm",troll,1,meta,"PC2score","PC3score",lab2,lab3)
#F2C2 <- fx_traits_vs_biomass_jitter("F2C2_hrm",troll,1,iso,"PC2score","PC3score",lab2,lab3)
#F2C3 <- fx_traits_vs_biomass_jitter("F2C3_hrm",troll,32,meta,"PC2score","PC3score",lab2,lab3)
#F2C3 <- fx_traits_vs_biomass_jitter("F2C3_hrm",troll,32,meta,"PC1score","PC2score",lab1,lab2)
F2C1 <- fx_traits_vs_biomass_jitter("F2C1_hrm",troll,1,meta,"PC2score","PC1score",lab2,lab1)
F2C2 <- fx_traits_vs_biomass_jitter("F2C2_hrm",troll,1,iso,"PC2score","PC1score",lab2,lab1)
F2C3 <- fx_traits_vs_biomass_jitter("F2C3_hrm",troll,32,meta,"PC2score","PC1score",lab2,lab1)
F2C4 <- fx_traits_vs_biomass_jitter("F2C4_hrm",troll,32,iso,"PC2score","PC1score",lab2,lab1)

#F2C1_P <- fx_traits_vs_biomass_jitter("F2C1_hrm_P",troll,1,meta,"PC2score","PC1score",lab2,lab1)
#F2C2_P <- fx_traits_vs_biomass_jitter("F2C2_hrm_P",troll,1,iso,"PC3score","PC1score",lab3,lab1)
#F2C3_P <- fx_traits_vs_biomass_jitter("F2C3_hrm_P",troll,32,meta,"PC1score","PC2score",lab1,lab2)
#F2C4_P <- fx_traits_vs_biomass_jitter("F2C4_hrm_P",troll,32,iso,"PC2score","PC1score",lab2,lab1)
F2C1_P <- fx_traits_vs_biomass_jitter("F2C1_hrm_P",troll,1,meta,"PC2score","PC1score",lab2,lab1)
F2C2_P <- fx_traits_vs_biomass_jitter("F2C2_hrm_P",troll,1,iso,"PC2score","PC1score",lab2,lab1)
F2C3_P <- fx_traits_vs_biomass_jitter("F2C3_hrm_P",troll,32,meta,"PC2score","PC1score",lab2,lab1)
F2C4_P <- fx_traits_vs_biomass_jitter("F2C4_hrm_P",troll,32,iso,"PC2score","PC1score",lab2,lab1)

### Dryland (Bjoern)
bjoern <- fx_read_model("readBjoern.R","bjoern") %>%
	select(-pRoot)

lab1 <- "maxBiomass" #(maximum \n size/size at maturity) [gC]"
lab2 <- "leafAllocation" #"pLeaf (allocation \n to leaf) [gC/gC]"
lab3 <- "storageAllocation" #"pStorage (allocation \n to storage) [gC/gC]"
lab4 <- "rootsAllocation" #"pRoot (allocation \n to root) [gC/gC]"

DC1 <- fx_traits_vs_biomass_jitter("DC1",bjoern,1,meta,"maxSize","pLeaf",lab1,lab2)
DC2 <- fx_traits_vs_biomass_jitter("DC2",bjoern,1,iso,"maxSize","pLeaf",lab1,lab2)
#DC2 <- fx_traits_vs_biomass_jitter("DC2",bjoern,1,iso,"maxSize","pStorage",lab1,lab3)
DC3 <- fx_traits_vs_biomass_jitter("DC3",bjoern,32,meta,"maxSize","pLeaf",lab1,lab2)
DC4 <- fx_traits_vs_biomass_jitter("DC4",bjoern,32,iso,"maxSize","pLeaf",lab1,lab2)

DC1_P <- fx_traits_vs_biomass_jitter("DC1_P",bjoern,1,meta,"maxSize","pLeaf",lab1,lab2)
#DC1_P <- fx_traits_vs_biomass_jitter("DC1_P",bjoern,1,meta,"pLeaf","maxSize",lab2,lab1)
DC2_P <- fx_traits_vs_biomass_jitter("DC2_P",bjoern,1,iso,"maxSize","pLeaf",lab1,lab2)
DC3_P <- fx_traits_vs_biomass_jitter("DC3_P",bjoern,32,meta,"maxSize","pLeaf",lab1,lab2)
DC4_P <- fx_traits_vs_biomass_jitter("DC4_P",bjoern,32,iso,"maxSize","pLeaf",lab1,lab2)
