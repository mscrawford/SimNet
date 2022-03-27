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
d1 <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Grass1_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d1_P <- fx_run_cforest(model,modelName,fileName)

### Grass2 (Lindsay's model)
if(SAVE_CACHE){model <- fx_read_model("readLindsay.R","Grass2")
}else{}

#Biomass
modelName = "Grass2"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2 <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Grass2_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_P <- fx_run_cforest(model,modelName,fileName)

### Grass2 (Lindsay's model) - variable
if(SAVE_CACHE){model <- fx_read_model("readLindsay_variable.R","Grass2") 
}else{}

#Biomass
modelName = "Grass2_v"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_v <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Grass2_v_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_v_P <- fx_run_cforest(model,modelName,fileName)

### Grass3 (IBC-grass)
if(SAVE_CACHE){model <- fx_read_model("readIBC.R","Grass3")
#model <- fx_read_model("readIBC.R","IBC_grass.noNDD")
}else{}

#Biomass
modelName = "Grass3"
#modelName = "Grass3_noNDD"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3 <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Grass3_P"
#modelName = "Grass3_P_noNDD"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3_P <- fx_run_cforest(model,modelName,fileName)

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
d3_3 <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Grass3_PCA3_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3_3_P <- fx_run_cforest(model,modelName,fileName)

### Forest1 (PPA)
if(SAVE_CACHE){model <- fx_read_model("readPPA.R","Forest1")
}else{}

#Biomass
modelName = "Forest1"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d4 <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Forest1_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d4_P <- fx_run_cforest(model,modelName,fileName)

### Forest2 (TROLL)
if(SAVE_CACHE){model <- fx_read_model("readTROLL.R","Forest2")
}else{}

#Biomass
modelName = "Forest2"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5 <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Forest2_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5_P <- fx_run_cforest(model,modelName,fileName)

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
d5h_2 <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Forest2_PCA_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_2_P <- fx_run_cforest(model,modelName,fileName)

### Forest2 (TROLL) h_realmax
if(SAVE_CACHE){model <- fx_read_model("readTROLL.R","Forest2") %>%
		mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
		select(-hmax, -ah, -dmax)
}else{}

#Biomass
modelName = "Forest2_hrealmax"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Forest2_hrealmax_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_P <- fx_run_cforest(model,modelName,fileName)
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
d5h_3 <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Forest2_hrealmax_PCA_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_3_P <- fx_run_cforest(model,modelName,fileName)

### Dryland (Bjoern)
if(SAVE_CACHE){model <- fx_read_model("readBjoern.R","bjoern") %>%
		select(-pRoot)
}else{}

#Biomass
modelName = "Dryland"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d6 <- fx_run_cforest(model,modelName,fileName)

#Productivity
modelName = "Dryland_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d6_P <- fx_run_cforest(model,modelName,fileName)

# All models
#Biomass
resvar = "Biomass"
pN1 = "cforest_grid"
pN2 = "cforest_grid_meta"
pN3 = "cforest_grid_2PC"
pN4 = "cforest_grid_meta_2PC"
pN5 = "cforest_grid_3PC"
pN6 = "cforest_grid_meta_3PC"
all_d <- rbind(d1,d2,d3,d4,d5h,d6)
all_d <- fx_edit_final_df(all_d)
fx_plot_all(all_d,resvar,pN1)
all_d <- all_d[all_d$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d,resvar,pN2)

## PCA 3 components for Grass 3
all_d <- rbind(d1,d2,d3_3,d4,d5h_3,d6)
all_d <- fx_edit_final_df(all_d)
fx_plot_all(all_d,resvar,pN5)
all_d <- all_d[all_d$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d,resvar,pN6)

## PCA 3 components for Grass 3, Grass 2 variable
all_d <- rbind(d1,d2_v,d3_3,d4,d5h_3,d6)
all_d <- fx_edit_final_df(all_d)
fx_plot_all(all_d,resvar,pN5)
all_d <- all_d[all_d$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d,resvar,paste0(pN6,"_v"))

#Productivity
resvar = "Productivity"
all_d_P <- rbind(d1_P,d2_P,d3_P,d4_P,d5h_P,d6_P)
all_d_P <- fx_edit_final_df(all_d_P)
fx_plot_all(all_d_P,resvar,paste0(pN1,"_P"))
all_d_P <- all_d_P[all_d_P$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d_P,resvar,paste0(pN2,"_P"))

## PCA 3 components for Grass 3
all_d_P <- rbind(d1_P,d2_P,d3_3_P,d4_P,d5h_3_P,d6_P)
all_d_P <- fx_edit_final_df(all_d_P)
fx_plot_all(all_d_P,resvar,paste0(pN5,"_P"))
all_d_P <- all_d_P[all_d_P$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d_P,resvar,paste0(pN6,"_P"))

## PCA 3 components for Grass 3, Grass 2 variable
all_d_P <- rbind(d1_P,d2_v_P,d3_3_P,d4_P,d5h_3_P,d6_P)
all_d_P <- fx_edit_final_df(all_d_P)
fx_plot_all(all_d_P,resvar,paste0(pN5,"_v_P"))
all_d_P <- all_d_P[all_d_P$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d_P,resvar,paste0(pN6,"_v_P"))
