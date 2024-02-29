library(png)
library(readr)

base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")
store_dir         <- paste0(tmp_dir,"/traits_vs_biomass/")

#source(paste0(scripts_dir, "/readModels.R"))
source(paste0(scripts_dir, "/to_test/fx_traits_vs_biomass.R"))
source(paste0(scripts_dir, "/to_test/fx_cforest_party.R"))

#READ_CACHE <- TRUE 
#SAVE_CACHE <- FALSE
READ_CACHE <- FALSE
SAVE_CACHE <- TRUE 

### Grass1 (Adam's model)
file_name <- paste0(store_dir,"adam.csv")
if(SAVE_CACHE){adam <- fx_read_model("readAdam.R", "Grass1") %>%
	       mutate(r_pNi = 1/pNi) #rpNi = reciprocal of pNi
               write.csv(adam, file_name, row.names=FALSE)
}else{adam <- read_csv(file_name)}

lab1 <- "monoBiomass" #"no3i (nitrogen R*)"
lab2 <- "NUE1" #"r_pNi (reciprocal of \n aboveground N concentration)"

d1 <- fx_prepare_df("Grass 1",adam,"abmi","r_pNi",lab1,lab2) 
fx_scatter_plot("mono",d1,lab1,"G1C1") 
fx_scatter_plot("mix",d1,lab1,"G1C2") 
fx_scatter_plot("mono",d1,lab2,"G1C3") 
fx_scatter_plot("mix",d1,lab2,"G1C4") 

### Grass2 (Lindsay's model)

file_name <- paste0(store_dir,"lindsay.csv")
if(SAVE_CACHE){lindsay <- fx_read_model("readLindsay_variable.R","Grass2") %>%
	       mutate(thetai = 1/thetai) # reciprocal of thetai (NUE2)
               write.csv(lindsay, file_name, row.names=FALSE)
}else{lindsay <- read_csv(file_name)}
#formatC(lindsay["Vi"], format = "e", digits = 0)

lab1 <- "rootingVolume" #"Vi (volume of soil \n accessible to species i)"
lab2 <- "NUE2" #thetai (Nitrogen uptake \n rate per unit plant biomass)"

d2 <- fx_prepare_df("Grass 2",lindsay,"Vi","thetai",lab1,lab2)
fx_scatter_plot("mono",d2,lab1,"G2C1") 
fx_scatter_plot("mix",d2,lab1,"G2C2") 
fx_scatter_plot("mono",d2,lab2,"G2C3") 
fx_scatter_plot("mix",d2,lab2,"G2C4") 

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

lab1 <- "LES1" #"PC1score associated with LMR \n (leaf to mass ratio)"
lab2 <- "Size/Growth" #"PC2score associated with Gmax \n (maximum resource utilization) \n and MaxMass (Plant's maximum size)"
lab3 <- "Spacing" #"PC3score associated with SLA and MeanSpacerLength"

d3 <- fx_prepare_df("Grass 3",IBC_grass,"PC2score","PC1score",lab2,lab1)
fx_scatter_plot("mono",d3,lab1,"G3C1") 
fx_scatter_plot("mix",d3,lab1,"G3C2") 
fx_scatter_plot("mono",d3,lab2,"G3C3") 
fx_scatter_plot("mix",d3,lab2,"G3C4") 

#### Forest1 (PPA)

file_name <- paste0(store_dir,"PPA.csv")
if(SAVE_CACHE){PPA <- fx_read_model("readPPA.R","Forest1")
               write.csv(PPA, file_name, row.names=FALSE)
}else{PPA <- read_csv(file_name)}

lab1 <- "GrowthSurvival" #"paceOfLife" #"PC1score (associated \n with fast-slow lifecycle)"#plant height)"
lab2 <- "MaxHeight" #"PC2score (associated \n with tree stature)"# LMA -leaf mass per area)"

d4 <- fx_prepare_df("Forest 1",PPA,"PC2score","PC1score",lab2,lab1)
fx_scatter_plot("mono",d4,lab1,"F1C1") 
fx_scatter_plot("mix",d4,lab1,"F1C2") 
fx_scatter_plot("mono",d4,lab2,"F1C3") 
fx_scatter_plot("mix",d4,lab2,"F1C4") 

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
lab1 <- "LES2" #"PC1score associated with \n LMA, nmass, and pmass"
lab2 <- "MaxHeight" #"PC2score associated with \n h_realmax = hmax * dmax / (dmax + ah)"
lab3 <- "woodDensity" #"PC3score associated with \n wsg (wood specific gravity)"

d5 <- fx_prepare_df("Forest 2",troll,"PC2score","PC1score",lab2,lab1)
fx_scatter_plot("mono",d5,lab1,"F2C1") 
fx_scatter_plot("mix",d5,lab1,"F2C2") 
fx_scatter_plot("mono",d5,lab2,"F2C3") 
fx_scatter_plot("mix",d5,lab2,"F2C4") 

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

d6 <- fx_prepare_df("Dryland",bjoern,"maxSize","pLeaf",lab1,lab2)
fx_scatter_plot("mono",d6,lab1,"DC1") 
fx_scatter_plot("mix",d6,lab1,"DC2") 
fx_scatter_plot("mono",d6,lab2,"DC3") 
fx_scatter_plot("mix",d6,lab2,"DC4") 

#########################################################################################
#########################################################################################
#########################################################################################

#all_d <- rbind(d1,d2,d3,d4,d5,d6)
#plot <- fx_plot_all(all_d)
#plot <- fx_single_plot("Grass 1",adam,"abmi","r_pNi",lab1,lab2) 
#plot <- fx_biomass_vs_traits(all_d)

#### Compile figure 3 in the manuscript: group all plots: 
dir <- paste0(tmp_dir,"/traits_vs_biomass/")
G1C1 <- readPNG(paste0(dir,"G1C1.png"))
G1C2 <- readPNG(paste0(dir,"G1C2.png"))
G1C3 <- readPNG(paste0(dir,"G1C3.png"))
G1C4 <- readPNG(paste0(dir,"G1C4.png"))
G2C1 <- readPNG(paste0(dir,"G2C1.png"))
G2C3 <- readPNG(paste0(dir,"G2C3.png"))
G3C1 <- readPNG(paste0(dir,"G3C1.png"))
G3C3 <- readPNG(paste0(dir,"G3C3.png"))
F1C1 <- readPNG(paste0(dir,"F1C1.png"))
F1C3 <- readPNG(paste0(dir,"F1C3.png"))
F2C1 <- readPNG(paste0(dir,"F2C1.png"))
F2C3 <- readPNG(paste0(dir,"F2C3.png"))
G2C2 <- readPNG(paste0(dir,"G2C2.png"))
G2C4 <- readPNG(paste0(dir,"G2C4.png"))
G3C2 <- readPNG(paste0(dir,"G3C2.png"))
G3C4 <- readPNG(paste0(dir,"G3C4.png"))
F1C2 <- readPNG(paste0(dir,"F1C2.png"))
F1C4 <- readPNG(paste0(dir,"F1C4.png"))
F2C2 <- readPNG(paste0(dir,"F2C2.png"))
F2C4 <- readPNG(paste0(dir,"F2C4.png"))
#F2C1 <- readPNG(paste0(dir,"F2C1_hrm.png"))
#F2C3 <- readPNG(paste0(dir,"F2C3_hrm.png"))
DC1 <- readPNG(paste0(dir,"DC1.png"))
DC3 <- readPNG(paste0(dir,"DC3.png"))
DC2 <- readPNG(paste0(dir,"DC2.png"))
DC4 <- readPNG(paste0(dir,"DC4.png"))
#F1C1_P <- readPNG(paste0(dir,"F1C1_P.png"))
#F1C3_P <- readPNG(paste0(dir,"F1C3_P.png"))
#F2C1_P <- readPNG(paste0(dir,"F2C1_hrm_P.png"))
#F2C3_P <- readPNG(paste0(dir,"F2C3_hrm_P.png"))
#DC1_P <- readPNG(paste0(dir,"DC1_P.png"))
#DC3_P <- readPNG(paste0(dir,"DC3_P.png"))

#png(paste0(dir, "all_scatter_plots_biomass.png"), 
png(paste0(dir, "all_scatter_plots.png"), 
    #width = 1663, height = 3061, units = "px", pointsize = 40,  res = NA,
    width = 3326, height = 3061, units = "px", pointsize = 40,  res = NA,
    bg = "white", type = c("cairo", "cairo-png", "Xlib", "quartz"))
par(mar=c(0,0,2,0))
#all_scatter <- plot(NA, xlim = c(0, 2), ylim = c(0, 6), type = "n", 
all_scatter <- plot(NA, xlim = c(0, 4), ylim = c(0, 6), type = "n", 
		    xaxt = "n", yaxt = "n", xlab = "", ylab = "Model",
                    #main = "Monoculture                                                                                                         Mixture")
                    main = "Monoculture                                    Mixture                                       Monoculture                                    Mixture")
               #+ theme_bw() + geom_line()

#rasterImage(image, xleft, ybottom, xright, ytop, angle = 0, interpolate = TRUE, …)

rasterImage(G1C1, 0, 5, 1, 6)
rasterImage(G1C2, 1, 5, 2, 6)
rasterImage(G1C3, 2, 5, 3, 6)
rasterImage(G1C4, 3, 5, 4, 6)
rasterImage(G2C1, 0, 4, 1, 5)
rasterImage(G2C2, 1, 4, 2, 5)
rasterImage(G2C3, 2, 4, 3, 5)
rasterImage(G2C4, 3, 4, 4, 5)
rasterImage(G3C1, 0, 3, 1, 4)
rasterImage(G3C2, 1, 3, 2, 4)
rasterImage(G3C3, 2, 3, 3, 4)
rasterImage(G3C4, 3, 3, 4, 4)
rasterImage(F1C1, 0, 2, 1, 3)
rasterImage(F1C2, 1, 2, 2, 3)
rasterImage(F1C3, 2, 2, 3, 3)
rasterImage(F1C4, 3, 2, 4, 3)
rasterImage(F2C1, 0, 1, 1, 2)
rasterImage(F2C2, 1, 1, 2, 2)
rasterImage(F2C3, 2, 1, 3, 2)
rasterImage(F2C4, 3, 1, 4, 2)
rasterImage(DC1, 0, 0, 1, 1)
rasterImage(DC2, 1, 0, 2, 1)
rasterImage(DC3, 2, 0, 3, 1)
rasterImage(DC4, 3, 0, 4, 1)
while (!is.null(dev.list()))  dev.off()

#path <- paste0(tmp_dir,"/traits_vs_biomass/")
#png(paste0(path,"all_scatter_plots_productivity.png"), 
#    width = 1663, height = 1530.5, units = "px", pointsize = 40,  res = NA,
#    bg = "white", type = c("cairo", "cairo-png", "Xlib", "quartz"))
#par(mar=c(0,0,1,0))
#all_scatter <- plot(NA, xlim = c(0, 2), ylim = c(0, 3), type = "n", 
#		    xaxt = "n", yaxt = "n", xlab = "", ylab = "",
#                    main = "Monoculture                                   Mixture")
#               #+ theme_bw() + geom_line()
#
##rasterImage(image, xleft, ybottom, xright, ytop, angle = 0, interpolate = TRUE, …)
#
#rasterImage(F1C1_P, 0, 2, 1, 3)
#rasterImage(F1C3_P, 1, 2, 2, 3)
#rasterImage(F2C1_P, 0, 1, 1, 2)
#rasterImage(F2C3_P, 1, 1, 2, 2)
#rasterImage(DC1_P, 0, 0, 1, 1)
#rasterImage(DC3_P, 1, 0, 2, 1)
#while (!is.null(dev.list()))  dev.off()
#
