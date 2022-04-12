library(gridExtra)

set.seed(1987)
#base_dir          <- setwd("~/Documents/SimNet")
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

source(paste0(scripts_dir, "/to_test/fx_cforest_party.R"))

#READ_CACHE <- FALSE
#SAVE_CACHE <- TRUE
READ_CACHE <- TRUE 
SAVE_CACHE <- FALSE

if(READ_CACHE){model <- c()}else{}

### Grass1 (Adam's model)
if(SAVE_CACHE){model <- fx_read_model("readAdam.R", "Grass1") %>%
		mutate(r_pNi = 1/pNi) %>% #rpNi = reciprocal of pNi
		select(-no3i,-pNi)
}else{}

#Biomass
modelName = "Grass1"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d1 <- fx_cforest(model,modelName,fileName)
d1[d1 == "abmi"] <- "monoBiomass"
d1[d1 == "r_pNi"] <- "NUE1" #"N use efficiency"

#Productivity
modelName = "Grass1_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d1_P <- fx_cforest(model,modelName,fileName)
d1_P[d1_P == "abmi"] <- "monoBiomass"
d1_P[d1_P == "r_pNi"] <- "NUE1" #"N use efficiency"

### Grass2 (Lindsay's model)
if(SAVE_CACHE){model <- fx_read_model("readLindsay.R","Grass2")
}else{}

#Biomass
modelName = "Grass2"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2 <- fx_cforest(model,modelName,fileName)

#Productivity
modelName = "Grass2_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_P <- fx_cforest(model,modelName,fileName)

### Grass2 (Lindsay's model) - variable
if(SAVE_CACHE){model <- fx_read_model("readLindsay_variable.R","Grass2") 
}else{}

#Biomass
modelName = "Grass2_v"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_v <- fx_cforest(model,modelName,fileName)
d2_v[d2_v == "thetai"] <- "NUE2" #"N uptake rate"
d2_v[d2_v == "Vi"] <- "rootingVolume"

#Productivity
modelName = "Grass2_v_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_v_P <- fx_cforest(model,modelName,fileName)
d2_v_P[d2_v_P == "thetai"] <- "NUE2" #"N uptake rate"
d2_v_P[d2_v_P == "Vi"] <- "rootingVolume"

### Grass3 (IBC-grass)
if(SAVE_CACHE){model <- fx_read_model("readIBC.R","Grass3")
#model <- fx_read_model("readIBC.R","IBC_grass.noNDD")
}else{}

#Biomass
modelName = "Grass3"
#modelName = "Grass3_noNDD"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3 <- fx_cforest(model,modelName,fileName)

#Productivity
modelName = "Grass3_P"
#modelName = "Grass3_P_noNDD"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3_P <- fx_cforest(model,modelName,fileName)

#READ_CACHE <- FALSE
#SAVE_CACHE <- TRUE
### Grass3 (IBC-grass) PCA - 3 components
if(SAVE_CACHE){model <- readRDS(paste0(tmp_dir,"/PCA/Grass3_PCAcoord.Rda")) %>%
	select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, PC1score, PC2score, PC3score) %>%
	mutate(id = row_number()) %>%
	mutate_if(is.character, as.factor) %>%
	mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
	mutate(Productivity = scales::rescale(Productivity, to = c(0, 100)))

}else{} 

#Biomass
modelName = "Grass3_PCA3"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3_3 <- fx_cforest(model,modelName,fileName)
d3_3[d3_3 == "PC1score"] <- "LES1" #"leaf eco. spectrum"	
d3_3[d3_3 == "PC2score"] <- "Size/Growth"	
d3_3[d3_3 == "PC3score"] <- "Spacing"

#Productivity
modelName = "Grass3_PCA3_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3_3_P <- fx_cforest(model,modelName,fileName)
d3_3_P[d3_3_P == "PC1score"] <- "LES1" #"leaf eco. spectrum"
d3_3_P[d3_3_P == "PC2score"] <- "Size/Growth"	
d3_3_P[d3_3_P == "PC3score"] <- "Spacing"
#READ_CACHE <- TRUE 
#SAVE_CACHE <- FALSE

### Forest1 (PPA)
if(SAVE_CACHE){model <- fx_read_model("readPPA.R","Forest1")
}else{}

#Biomass
modelName = "Forest1"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d4 <- fx_cforest(model,modelName,fileName)
d4[d4 == "PC1score"] <- "paceOfLife"	
d4[d4 == "PC2score"] <- "MaxHeight"	

#Productivity
modelName = "Forest1_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d4_P <- fx_cforest(model,modelName,fileName)
d4_P[d4_P == "PC1score"] <- "paceOfLife"	
d4_P[d4_P == "PC2score"] <- "MaxHeight"	

### Forest2 (TROLL)
if(SAVE_CACHE){model <- fx_read_model("readTROLL.R","Forest2")
}else{}

#Biomass
modelName = "Forest2"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5 <- fx_cforest(model,modelName,fileName)

#Productivity
modelName = "Forest2_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5_P <- fx_cforest(model,modelName,fileName)

### Forest2 (TROLL) PCA - 3 components
if(SAVE_CACHE){model <- readRDS(paste0(tmp_dir,"/PCA/Forest2_PCAcoord.Rda")) %>%
	select(c(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, PC1score, PC2score, PC3score)) %>%
	mutate(id = row_number()) %>%
	mutate_if(is.character, as.factor) %>%
	mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
	mutate(Productivity = scales::rescale(Productivity, to = c(0, 100)))

}else{} 

#Biomass
modelName = "Forest2_PCA"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_2 <- fx_cforest(model,modelName,fileName)

#Productivity
modelName = "Forest2_PCA_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_2_P <- fx_cforest(model,modelName,fileName)

### Forest2 (TROLL) h_realmax
if(SAVE_CACHE){model <- fx_read_model("readTROLL.R","Forest2") %>%
		mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
		select(-hmax, -ah, -dmax)
}else{}

#Biomass
modelName = "Forest2_hrealmax"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h <- fx_cforest(model,modelName,fileName)

#Productivity
modelName = "Forest2_hrealmax_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_P <- fx_cforest(model,modelName,fileName)
##See what the data frame looks like for only monoculture-isolation:
#t2 <- model[which(model$Stage == 'Without seed inflow' & model$Ninitial == 1), ]
#mean(t2$Biomass)
#sd(t2$Biomass)
##Compare with respect to Mixture-isolation
#t3 <- model[which(model$Stage == 'Without seed inflow' & model$Ninitial == 32), ]
#mean(t3$Biomass)
#sd(t3$Biomass)

### Forest2 (TROLL) h_realmax PCA - 3 components
if(SAVE_CACHE){model <- readRDS(paste0(tmp_dir,"/PCA/Forest2_hrm_PCAcoord.Rda")) %>%
	select(c(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, PC1score, PC2score, PC3score)) %>%
	mutate(id = row_number()) %>%
	mutate_if(is.character, as.factor) %>%
	mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
	mutate(Productivity = scales::rescale(Productivity, to = c(0, 100)))

}else{} 

#Biomass
modelName = "Forest2_hrealmax_PCA"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_3 <- fx_cforest(model,modelName,fileName)
d5h_3[d5h_3 == "PC1score"] <- "LES2" #"Leaf res. conserv."
d5h_3[d5h_3 == "PC2score"] <- "MaxHeight"	
d5h_3[d5h_3 == "PC3score"] <- "woodDensity"

#Productivity
modelName = "Forest2_hrealmax_PCA_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_3_P <- fx_cforest(model,modelName,fileName)
d5h_3_P[d5h_3_P == "PC1score"] <- "LES2" #"Leaf res. conserv."
d5h_3_P[d5h_3_P == "PC2score"] <- "MaxHeight"	
d5h_3_P[d5h_3_P == "PC3score"] <- "woodDensity"

### Dryland (Bjoern)
if(SAVE_CACHE){model <- fx_read_model("readBjoern.R","bjoern") %>%
		select(-pRoot)
}else{}

#Biomass
modelName = "Dryland"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d6 <- fx_cforest(model,modelName,fileName)
d6[d6 == "maxSize"] <- "maxBiomass"
d6[d6 == "pLeaf"] <- "leafAllocation"
d6[d6 == "pStorage"] <- "storageAllocation"

#Productivity
modelName = "Dryland_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d6_P <- fx_cforest(model,modelName,fileName)
d6_P[d6_P == "maxSize"] <- "maxBiomass"
d6_P[d6_P == "pLeaf"] <- "leafAllocation"
d6_P[d6_P == "pStorage"] <- "storageAllocation"

## All models
##Biomass
#resvar = "Biomass"
#pN1 = "cforest_grid"
#pN2 = "cforest_grid_meta"
#pN3 = "cforest_grid_2PC"
#pN4 = "cforest_grid_meta_2PC"
#pN5 = "cforest_grid_3PC"
#pN6 = "cforest_grid_meta_3PC"
#all_d <- rbind(d1,d2,d3,d4,d5h,d6)
#fx_plot_all(all_d,resvar,pN1)
#all_d <- all_d[all_d$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
#fx_plot_all(all_d,resvar,pN2)
#
### PCA 3 components for Grass 3
#all_d <- rbind(d1,d2,d3_3,d4,d5h_3,d6)
#fx_plot_all(all_d,resvar,pN5)
#all_d <- all_d[all_d$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
#fx_plot_all(all_d,resvar,pN6)
#
### PCA 3 components for Grass 3, Grass 2 variable
#all_d <- rbind(d1,d2_v,d3_3,d4,d5h_3,d6)
#fx_plot_all(all_d,resvar,pN5)
#all_d <- all_d[all_d$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
#all_d[all_d == "Mono.-Meta."] <- "Monoculture"
#all_d[all_d == "Mix.-Meta."] <- "Mixture"
#bio_df <- all_d
#fx_plot_all(all_d,resvar,paste0(pN6,"_v"))
#print(head(all_d))
#
##Productivity
#resvar = "Productivity"
#all_d_P <- rbind(d1_P,d2_P,d3_P,d4_P,d5h_P,d6_P)
#fx_plot_all(all_d_P,resvar,paste0(pN1,"_P"))
#all_d_P <- all_d_P[all_d_P$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
#fx_plot_all(all_d_P,resvar,paste0(pN2,"_P"))
#
### PCA 3 components for Grass 3
#all_d_P <- rbind(d1_P,d2_P,d3_3_P,d4_P,d5h_3_P,d6_P)
#fx_plot_all(all_d_P,resvar,paste0(pN5,"_P"))
#all_d_P <- all_d_P[all_d_P$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
#fx_plot_all(all_d_P,resvar,paste0(pN6,"_P"))
#
### PCA 3 components for Grass 3, Grass 2 variable
#all_d_P <- rbind(d1_P,d2_v_P,d3_3_P,d4_P,d5h_3_P,d6_P)
#fx_plot_all(all_d_P,resvar,paste0(pN5,"_v_P"))
#all_d_P <- all_d_P[all_d_P$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
#all_d_P[all_d_P == "Mono.-Meta."] <- "Monoculture"
#all_d_P[all_d_P == "Mix.-Meta."] <- "Mixture"
#prod_df <- all_d_P
#fx_plot_all(all_d_P,resvar,paste0(pN6,"_v_P"))
#print(head(all_d_P))

# Difference plots Biomass
d1 <- fx_diff(d1)
d2_v <- fx_diff(d2_v)
#d3_3 <- fx_diff(d3_3)
d4 <- fx_diff(d4)
d5h_3<- fx_diff(d5h_3)
d6 <- fx_diff(d6)
df <- rbind(d1,d2_v,d4,d5h_3,d6)
#df <- rbind(d1,d2_v,d3_3,d4,d5h_3,d6)
plotName <- 'Diff_import'
fx_plot_diff_mono_mix(plotName,df)

d1_P <- fx_diff(d1_P)
d2_v_P <- fx_diff(d2_v_P)
#d3_3_P <- fx_diff(d3_3_P)
d4_P <- fx_diff(d4_P)
d5h_3_P<- fx_diff(d5h_3_P)
d6_P <- fx_diff(d6_P)
df_P <- rbind(d1_P,d2_v_P,d4_P,d5h_3_P,d6_P)
print(df_P)
#df <- rbind(d1_P,d2_v_P,d3_3_P,d4_P,d5h_3_P,d6_P)
plotName <- 'Diff_import_P'
fx_plot_diff_mono_mix(plotName,df_P)
