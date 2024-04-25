library(ggplot2)
library(ggh4x)
library(png)
library(readr)
library(data.table)
library(tidyr)
library(dplyr)
#library(ggforce)
#library(ggbreak)
#library(patchwork)
#library(cowplot)

base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")
store_dir         <- paste0(tmp_dir,"/traits_vs_biomass/")

source(paste0(scripts_dir, "/to_test/fx_cforest_party.R"))

READ_CACHE <- TRUE 
SAVE_CACHE <- FALSE
#READ_CACHE <- FALSE
#SAVE_CACHE <- TRUE 

fx_prepare_df  <- function(modelName,data,trait1,trait2){
	dfmono <- data %>%
	filter(Ninitial == 1, Year %in% meta) %>%
	select(SpeciesID, Biomass, trait1, trait2)
    dfmono$type_bm <- rep("Monoculture",times=dim(dfmono)[1])

	dfmix <- data %>%
	filter(Ninitial == 32, Year %in% meta) %>%
	select(SpeciesID, Biomass, trait1, trait2) #%>%
	dfmix <- aggregate(. ~ SpeciesID, data = dfmix, mean) #Calculate mean response per specie
	dfmix$type_bm <- rep("Mixture",times=dim(dfmix)[1])

	df <- rbind(dfmono,dfmix) %>%
	setnames(old=c(trait1,trait2), new=c("Trait 1","Trait 2"))

	df$Model <- rep(modelName,times=dim(df)[1]) 
    print(head(df))
    return(df)
}

### Grass1 (Adam's model)
file_name <- paste0(store_dir,"adam.csv")
if(SAVE_CACHE){adam <- fx_read_model("readAdam.R", "Grass1") %>%
	       mutate(r_pNi = 1/pNi) #rpNi = reciprocal of pNi
               write.csv(adam, file_name, row.names=FALSE)
}else{adam <- read_csv(file_name)}

lab1 <- "monoBiomass" #"no3i (nitrogen R*)"
lab2 <- "NUE1" #"r_pNi (reciprocal of \n aboveground N concentration)"

d1 <- fx_prepare_df("Grass 1",adam,"abmi","r_pNi") 
D1 <- data.frame(var=c("Trait 1","Trait 2"), Lab = c(lab1,lab2))
D1$Model <- rep("Grass 1",times=2)

### Grass2 (Lindsay's model)

file_name <- paste0(store_dir,"lindsay.csv")
if(SAVE_CACHE){lindsay <- fx_read_model("readLindsay_variable.R","Grass2") %>%
	       mutate(thetai = 1/thetai) # reciprocal of thetai (NUE2)
               write.csv(lindsay, file_name, row.names=FALSE)
}else{lindsay <- read_csv(file_name)}
#formatC(lindsay["Vi"], format = "e", digits = 0)

lab1 <- "rootingVolume" #"Vi (volume of soil \n accessible to species i)"
lab2 <- "NUE2" #thetai (Nitrogen uptake \n rate per unit plant biomass)"

d2 <- fx_prepare_df("Grass 2",lindsay,"Vi","thetai")
D2 <- data.frame(var=c("Trait 1","Trait 2"), Lab = c(lab1,lab2))
D2$Model <- rep("Grass 2",times=2)

### Grass3 (IBC-grass)

file_name <- paste0(store_dir,"IBC_grass.csv")
if(SAVE_CACHE){IBC_grass <- readRDS(paste0(tmp_dir,"/PCA/Grass3_PCAcoord.Rda")) %>%
	       select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, PC1score, PC2score, PC3score) %>%
	       mutate(PC2score = -PC2score) %>%
	       mutate(PC3score = -PC3score) %>%
	       mutate(id = row_number()) %>%
	       mutate_if(is.character, as.factor)
               write.csv(IBC_grass, file_name, row.names=FALSE)
}else{IBC_grass <- read_csv(file_name)}

lab1 <- "Size/Growth" #"PC2score associated with Gmax \n (maximum resource utilization) \n and MaxMass (Plant's maximum size)"
lab2 <- "LES1" #"PC1score associated with LMR \n (leaf to mass ratio)"
lab3 <- "Spacing" #"PC3score associated with SLA and MeanSpacerLength"

d3 <- fx_prepare_df("Grass 3",IBC_grass,"PC2score","PC1score")
D3 <- data.frame(var=c("Trait 1","Trait 2"), Lab = c(lab1,lab2))
D3$Model <- rep("Grass 3",times=2)

#### Forest1 (PPA)

file_name <- paste0(store_dir,"PPA.csv")
if(SAVE_CACHE){PPA <- fx_read_model("readPPA.R","Forest1")
               write.csv(PPA, file_name, row.names=FALSE)
}else{PPA <- read_csv(file_name)}

lab1 <- "MaxHeight" #"PC2score (associated \n with tree stature)"# LMA -leaf mass per area)"
lab2 <- "GrowthSurvival" #"paceOfLife" #"PC1score (associated \n with fast-slow lifecycle)"#plant height)"

d4 <- fx_prepare_df("Forest 1",PPA,"PC2score","PC1score")
D4 <- data.frame(var=c("Trait 1","Trait 2"), Lab = c(lab1,lab2))
D4$Model <- rep("Forest 1",times=2)

### Forest2 (TROLL) h_realmax

file_name <- paste0(store_dir,"troll.csv")
if(SAVE_CACHE){troll <- readRDS(paste0(tmp_dir,"/PCA/Forest2_hrm_PCAcoord.Rda")) %>%
	       select(c(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, PC1score, PC2score, PC3score)) %>%
	       mutate(PC1score = -PC1score) %>%
	       mutate(PC3score = -PC3score) %>%
	       mutate(id = row_number()) %>%
	       mutate_if(is.character, as.factor)
               write.csv(troll, file_name, row.names=FALSE)
}else{troll <- read_csv(file_name)}
lab1 <- "MaxHeight" #"PC2score associated with \n h_realmax = hmax * dmax / (dmax + ah)"
lab2 <- "LES2" #"PC1score associated with \n LMA, nmass, and pmass"
lab3 <- "woodDensity" #"PC3score associated with \n wsg (wood specific gravity)"

d5 <- fx_prepare_df("Forest 2",troll,"PC2score","PC1score")
D5 <- data.frame(var=c("Trait 1","Trait 2"), Lab = c(lab1,lab2))
D5$Model <- rep("Forest 2",times=2)

### Dryland (Bjoern)

file_name <- paste0(store_dir,"bjoern.csv")
if(SAVE_CACHE){bjoern <- fx_read_model("readBjoern.R","bjoern") %>%
	       select(-pRoot) %>%
	       mutate(maxSize = log(maxSize))
               write.csv(bjoern, file_name, row.names=FALSE)
}else{bjoern <- read_csv(file_name)}

lab1 <- "log maxBiomass" #(maximum \n size/size at maturity) [gC]"
lab2 <- "leafAllocation" #"pLeaf (allocation \n to leaf) [gC/gC]"
lab3 <- "storageAllocation" #"pStorage (allocation \n to storage) [gC/gC]"
lab4 <- "rootsAllocation" #"pRoot (allocation \n to root) [gC/gC]"

d6 <- fx_prepare_df("Dryland",bjoern,"maxSize","pLeaf")
D6 <- data.frame(var=c("Trait 1","Trait 2"), Lab = c(lab1,lab2))
D6$Model <- rep("Dryland",times=2)

#########################################################################################
print("##############################    Scatter    ##########################")
#########################################################################################

df <- rbind(d1,d2,d3,d4,d5,d6)
labels <- rbind(D1,D2,D3,D4,D5,D6)

df$SpeciesID <- as.character(df$SpeciesID)
df$log_bm <- log(df$Biomass)

df <- df %>%
    pivot_longer(c(-SpeciesID, -Biomass, -type_bm, -Model, -log_bm), names_to= "var", values_to = "Value")
df <- merge(df,labels, by=c("Model","var"))
print(head(df))
df <- na.omit(df)

# Subset to plot loess line for all models, except Dryland
df_noD <- df[df$Model != "Dryland",]
#"geom_text() we need to assemble a data frame containing the text of the labels in one column and columns for the variables to be mapped to other aesthetics, as well as the variable(s) used for faceting." 
xlabs <- df %>%
    group_by(Model,var,type_bm,Lab) %>%
    summarise(log_bm = 1.36 * max(log_bm), Value = max(Value))

# Calculate a position below the x-axis dynamically
y_position <- min(df$log_bm) - 0.1 * diff(range(df$log_bm))

plot_name <- "Figure3.png"
level1 <- c("Grass 1", "Forest 1", "Grass 2", "Grass 3", "Forest 2", "Dryland")
level2 <- c("Monoculture","Mixture")
p <- ggplot(data = df, aes(x = Value, y = log_bm)) +
  geom_point() +
  stat_smooth(data = df_noD, aes(x = Value, y = log_bm), 
              fullrange = TRUE, color = "red", method = "loess", se = FALSE) +
 #geom_smooth(method="lm", fill=NA) +
  labs(y = "log biomass") +
  facet_nested(factor(Model, levels = level1) ~ var + factor(type_bm, levels = level2), 
               scales = "free", independent = "all", nest_line = TRUE) +
#  facet_grid2(factor(Model, levels = level1) ~ var + factor(type_bm, levels = level2), 
#              scales = "free", independent = "all", axes = "all") +
#  facet_wrap(factor(Model,levels=level1) ~ factor(type_bm,levels=level2) + var, scales = "free", ncol = 4) +
  geom_text(data = xlabs, aes(label = Lab), 
            hjust = 1, vjust = 0.5, size = 8) +  # Adjusted aesthetics
  scale_x_continuous(n.breaks = 4) +
  theme_bw() +
  theme_classic() +
  theme(text = element_text(size = 25), strip.background = element_blank(),
        strip.text = element_text(size = 25), strip.clip = "off",
        legend.text = element_text(size = 20), legend.title = element_text(size = 20),
        panel.spacing = unit(3, "lines"), strip.placement = "outside")
ggsave(file=paste0(store_dir, plot_name), width=20, height=22, dpi=300)
while (!is.null(dev.list()))  dev.off()

