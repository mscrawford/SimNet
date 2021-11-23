library(cowplot)
set.seed(1987)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

#source(paste0(scripts_dir, "/readModels.R"))
source(paste0(scripts_dir, "/to_test/fx_cforest_party.R"))

### Grass1 (Adam's model)
#source(paste0(scripts_dir, "/to_test/readAdam.R"))
#
#adam <- models$Grass1 %>%
# ungroup() %>% # There shouldn't be groups anyways
# mutate(Biomass = scales::rescale(Biomass,
#                                  to = c(0, 100)))%>%
# #select(-Model, -SeedRain, -abmi) #Removed abmi because it is highly correlated to biomass
# select(-Model, -SeedRain)
#
#model<-adam
#C3 <- fx_single_cond(model,32,meta)
#C4 <- fx_single_cond(model,32,iso)
## remove traits with negative importance values:
#model <- subset(model, select = -pNi)
#C1 <- fx_single_cond(model,1,meta)
#C2 <- fx_single_cond(model,1,iso)

modelName = "Grass1_m"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
rf_df <- readRDS(fileName)
#p1 <- fx_plot(rf_df,modelName)
d1 <- rf_df

### Grass2 (Lindsay's model)
#source(paste0(scripts_dir, "/to_test/readLindsay.R"))
#
#lindsay <- models$Grass2%>%
#   ungroup() %>% # There shouldn't be groups anyways
#   mutate(Biomass = scales::rescale(Biomass,
#                                    to = c(0, 100)))%>%
#   select(-Model, -SeedRain)

modelName = "Grass2_m"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#p2 <- fx_cforest(lindsay,modelName)
rf_df <- readRDS(fileName)
d2 <- rf_df

### Grass3 (IBC-grass)
#source(paste0(scripts_dir, "/to_test/readIBC.R"))
#
#IBC_grass <- models$Grass3 %>%
# ungroup() %>% # There shouldn't be groups anyways
# mutate(Biomass = scales::rescale(Biomass,
#                                  to = c(0, 100)))%>%
# select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)
#
#model<-IBC_grass
#C3 <- fx_single_cond(model,32,meta)
#C4 <- fx_single_cond(model,32,iso)
## remove traits with negative importance values:
#model<-subset(model, select=c(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax))
#C1 <- fx_single_cond(model,1,meta)
#C2 <- fx_single_cond(model,1,iso)

modelName = "Grass3_m"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
rf_df <- readRDS(fileName)
#p3 <- fx_plot(rf_df,modelName)
d3 <- rf_df

### Forest1 (PPA)
#source(paste0(scripts_dir, "/to_test/readPPA.R"))
#
#PPA <- models$Forest1 %>%
# ungroup() %>% # There shouldn't be groups anyways
# mutate(Biomass = scales::rescale(Biomass,
#                                  to = c(0, 100)))%>%
# select(-Model, -SeedRain)

modelName = "Forest1_m"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#p4 <- fx_cforest(PPA,modelName)
rf_df <- readRDS(fileName)
d4 <- rf_df

### Forest2 (TROLL)
#source(paste0(scripts_dir, "/to_test/readTROLL.R"))
#
#troll <- models$Forest2 %>%
#   ungroup() %>% # There shouldn't be groups anyways
#   mutate(Biomass = scales::rescale(Biomass,
#                                    to = c(0, 100)))%>%
#   select(-Model, -SeedRain)
#
#model<-troll
#C3 <- fx_single_cond(model,32,meta)
#C4 <- fx_single_cond(model,32,iso)
## remove traits with negative importance values:
#model<-subset(model, select=c(-lma,-ah))
#C2 <- fx_single_cond(model,1,iso)
#model<-subset(model, select=c(-hmax))
#C1 <- fx_single_cond(model,1,meta)

modelName = "Forest2_m"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
rf_df <- readRDS(fileName)
#p5 <- fx_plot(rf_df,modelName)
d5 <- rf_df

### Forest2 (TROLL) h_realmax
#source(paste0(scripts_dir, "/to_test/readTROLL.R"))

#troll_hrealmax <- models$Forest2 %>%
#    ungroup() %>% # There shouldn't be groups anyways
#    mutate(Biomass = scales::rescale(Biomass,
#                                     to = c(0, 100)))%>%
#    select(-Model, -SeedRain) %>%
#    mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
#    select(-hmax, -ah, -dmax)
#
#model<-troll_hrealmax
#C2 <- fx_single_cond(model,1,iso)
#C3 <- fx_single_cond(model,32,meta)
#C4 <- fx_single_cond(model,32,iso)
#model<-subset(model, select=-lma)
#C1 <- fx_single_cond(model,1,meta)

modelName = "Forest2_hrealmax_m"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
rf_df <- readRDS(fileName)
#p5h <- fx_plot(rf_df,modelName)
d5h <- rf_df

### Dryland (Bjoern)
#source(paste0(scripts_dir, "/to_test/readBjoern.R"))
#
#bjoern <- models$bjoern %>% 
#    ungroup() %>% # There shouldn't be groups anyways
#    mutate(Biomass = scales::rescale(Biomass,
#                                     to = c(0, 100)))%>%
#    select(-Model, -SeedRain) 
#
#model<-bjoern
#C1 <- fx_single_cond(model,1,meta)
#model<-subset(model, select=-pLeaf)
##model<-subset(model, select=c(-pLeaf,-pRoot))
#C2 <- fx_single_cond(model,1,iso)
#model<-subset(model, select=c(-pStorage))
#C3 <- fx_single_cond(model,32,meta)
#C4 <- fx_single_cond(model,32,iso)

modelName = "Dryland_m"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
rf_df <- readRDS(fileName)
#p6 <- fx_plot(rf_df,modelName)
d6 <- rf_df

# All models
all_d <- rbind(d1,d2,d3,d4,d5,d5h,d6) %>%
mutate(condition = recode(condition, "1 80" = "Monoculture-Metacommunity"
			  ,"1 180" = "Monoculture-Isolation"
			  ,"32 80" = "Mixture-Metacommunity"
			  ,"32 180" = "Mixture-Isolation")
       ,mName = recode(mName, "Grass1_m" = "Grass 1"
		     ,"Grass2_m" = "Grass 2"
		     ,"Grass3_m" = "Grass 3"
		     ,"Forest1_m" = "Forest 1"
		     ,"Forest2_m" = "Forest 2"
		     ,"Forest2_hrealmax_m" = "Forest 2h"
		     ,"Dryland_m" = "Dryland")
)

ggplot(all_d, aes(x=reorder(varnames, var_categ), y=sCPI, fill=type)) +
#ggplot(all_d, aes(x=varnames, y=sCPI, fill=type)) +
    geom_bar(position='dodge',stat='identity') +
            geom_text(aes(label=scientific(CPI, digits = 2),size=10)
    		  ,position = position_dodge(width = 1)
                      ,hjust=0,vjust=0.2,color="black",size=2
                      ,show.legend = FALSE,angle = 90) +
    ylab("Conditional permutation importance, scaled to correlation") +
    xlab("Traits") +
    scale_fill_brewer(palette = "Set1",name = "Trait type") +
    scale_color_brewer(palette = "Set1") +
    ylim(NA,1.25) +
    facet_grid(condition ~ mName, scales = "free_x", 
     space = "free_x") + # Let the width of facets vary and force all bars to have the same width.
    theme_bw() +
    theme(text = element_text(size = 8),legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file=paste0(tmp_dir,"/randomForest/all_models_cforest_grid.pdf")
#       , width=15, height=7, dpi=300
)
while (!is.null(dev.list()))  dev.off()
