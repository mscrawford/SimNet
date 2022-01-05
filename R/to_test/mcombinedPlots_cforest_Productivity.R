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
#ungroup() %>% # There shouldn't be groups anyways
#mutate(Biomass = scales::rescale(Biomass,
#                                 to = c(0, 100)))%>%
#select(-Model, -SeedRain)
#adam_B <- adam %>% select(-Productivity)
#adam_P <- adam %>% select(-Biomass)
#setnames(adam_P, "Productivity", "Biomass")
#
##Biomass
#model<-adam_B
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

##Productivity
#model<-adam_P
#C3 <- fx_single_cond(model,32,meta)
#C4 <- fx_single_cond(model,32,iso)
#model <- subset(model, select = -pNi)
#C1 <- fx_single_cond(model,1,meta)
#C2 <- fx_single_cond(model,1,iso)

modelName = "Grass1_m_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
rf_df <- readRDS(fileName)
#p1 <- fx_plot(rf_df,modelName)
d1_P <- rf_df

### Grass2 (Lindsay's model)
#source(paste0(scripts_dir, "/to_test/readLindsay.R"))
#
#lindsay <- models$Grass2%>%
#   ungroup() %>% # There shouldn't be groups anyways
#   mutate(Biomass = scales::rescale(Biomass,
#                                    to = c(0, 100)))%>%
#   select(-Model, -SeedRain)
# lindsay_B <- lindsay %>% select(-Productivity)
# lindsay_P <- lindsay %>% select(-Biomass)
# setnames(lindsay_P, "Productivity", "Biomass")

#Biomass
modelName = "Grass2_m"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#p2 <- fx_cforest(lindsay_B,modelName)
rf_df <- readRDS(fileName)
d2 <- rf_df

#Productivity
modelName = "Grass2_m_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#p2 <- fx_cforest(lindsay_P,modelName)
rf_df <- readRDS(fileName)
d2_P <- rf_df

### Grass3 (IBC-grass)
source(paste0(scripts_dir, "/to_test/readIBC.R"))

#IBC_grass <- models$Grass3 %>%
IBC_grass <- models$IBC_grass.noNDD %>%
 ungroup() %>% # There shouldn't be groups anyways
 mutate(Biomass = scales::rescale(Biomass,
                                  to = c(0, 100)))%>%
 select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)
 IBC_grass_B <- IBC_grass %>% select(-Productivity)
 IBC_grass_P <- IBC_grass %>% select(-Biomass)
 setnames(IBC_grass_P, "Productivity", "Biomass")

#Biomass
model<-IBC_grass_B
C3 <- fx_single_cond(model,32,meta)
C4 <- fx_single_cond(model,32,iso)
model<-subset(model, select=c(-SLA, -meanSpacerLength))
C1 <- fx_single_cond(model,1,meta)
C2 <- fx_single_cond(model,1,iso)

modelName = "Grass3_m_noNDD"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
#rf_df <- readRDS(fileName)
p3 <- fx_plot(rf_df,modelName)
d3 <- rf_df

#Productivity
model<-IBC_grass_P
C3 <- fx_single_cond(model,32,meta)
C4 <- fx_single_cond(model,32,iso)
model<-subset(model, select=c(-SLA, -meanSpacerLength))
C1 <- fx_single_cond(model,1,meta)
C2 <- fx_single_cond(model,1,iso)

modelName = "Grass3_m_P_noNDD"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
#rf_df <- readRDS(fileName)
p3 <- fx_plot(rf_df,modelName)
d3_P <- rf_df

### Forest1 (PPA)
#source(paste0(scripts_dir, "/to_test/readPPA.R"))
#
#PPA <- models$Forest1 %>%
# ungroup() %>% # There shouldn't be groups anyways
# mutate(Biomass = scales::rescale(Biomass,
#                                  to = c(0, 100)))%>%
# select(-Model, -SeedRain)
# PPA_B <- PPA %>% select(-Productivity)
# PPA_P <- PPA %>% select(-Biomass)
# setnames(PPA_P, "Productivity", "Biomass")
#
##Biomass
modelName = "Forest1_m"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#p4 <- fx_cforest(PPA_B,modelName)
rf_df <- readRDS(fileName)
d4 <- rf_df

#Productivity
modelName = "Forest1_m_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#p4 <- fx_cforest(PPA_P,modelName)
rf_df <- readRDS(fileName)
d4_P <- rf_df

### Forest2 (TROLL)
#source(paste0(scripts_dir, "/to_test/readTROLL.R"))
#
#troll <- models$Forest2 %>%
#  ungroup() %>% # There shouldn't be groups anyways
#  mutate(Biomass = scales::rescale(Biomass,
#                                   to = c(0, 100)))%>%
#  select(-Model, -SeedRain)
#troll_B <- troll %>% select(-Productivity)
#troll_P <- troll %>% select(-Biomass)
#setnames(troll_P, "Productivity", "Biomass")

##Biomass
#model<-troll_B
#C3 <- fx_single_cond(model,32,meta)
#C4 <- fx_single_cond(model,32,iso)
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

##Productivity
#model<-troll_P
#C3 <- fx_single_cond(model,32,meta)
#C4 <- fx_single_cond(model,32,iso)
#model<-subset(model, select=c(-pmass,-lma,-dmax))
#C2 <- fx_single_cond(model,1,iso)
#model<-subset(troll_P, select=c(-hmax,-ah,-wsg))
#C1 <- fx_single_cond(model,1,meta)

modelName = "Forest2_m_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
rf_df <- readRDS(fileName)
#p5 <- fx_plot(rf_df,modelName)
d5_P <- rf_df

### Forest2 (TROLL) h_realmax
source(paste0(scripts_dir, "/to_test/readTROLL.R"))

troll_hrealmax <- models$Forest2 %>%
   ungroup() %>% # There shouldn't be groups anyways
   mutate(Biomass = scales::rescale(Biomass,
                                    to = c(0, 100)))%>%
   select(-Model, -SeedRain) %>%
   mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
   select(-hmax, -ah, -dmax)
troll_hrealmax_B <- troll_hrealmax %>% select(-Productivity)
troll_hrealmax_P <- troll_hrealmax %>% select(-Biomass)
setnames(troll_hrealmax_P, "Productivity", "Biomass")

##Biomass
#model<-troll_hrealmax_B
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

#Productivity
model<-troll_hrealmax_P
C3 <- fx_single_cond(model,32,meta)
C4 <- fx_single_cond(model,32,iso)
model<-subset(model, select=-wsg)
C1 <- fx_single_cond(model,1,meta)
model<-subset(troll_hrealmax_P, select=c(-nmass,-lma,-h_realmax))
C2 <- fx_single_cond(model,1,iso)

#See what the data frame looks like for only monoculture-isolation:
t2 <- model[which(model$Stage == 'Without seed inflow' & model$Ninitial == 1), ]
mean(t2$Biomass)
sd(t2$Biomass)
#Compare with respect to Mixture-isolation
t3 <- model[which(model$Stage == 'Without seed inflow' & model$Ninitial == 32), ]
mean(t3$Biomass)
sd(t3$Biomass)

modelName = "Forest2_hrealmax_m_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
rf_df <- readRDS(fileName)
#p5h <- fx_plot(rf_df,modelName)
d5h_P <- rf_df

### Dryland (Bjoern)
#source(paste0(scripts_dir, "/to_test/readBjoern.R"))
#
#bjoern <- models$bjoern %>% 
#    ungroup() %>% # There shouldn't be groups anyways
#    mutate(Biomass = scales::rescale(Biomass,
#                                     to = c(0, 100)))%>%
#    select(-Model, -SeedRain) 
# bjoern_B <- bjoern %>% select(-Productivity)
# bjoern_P <- bjoern %>% select(-Biomass)
# setnames(bjoern_P, "Productivity", "Biomass")
#
##Biomass
#model<-bjoern_B
#C1 <- fx_single_cond(model,1,meta)
#model<-subset(model, select=-pLeaf)
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

##Productivity
#model<-bjoern_P
#C1 <- fx_single_cond(model,1,meta)
#model<-subset(model, select=-pRoot)
#C2 <- fx_single_cond(model,1,iso)
#C4 <- fx_single_cond(model,32,iso)
#model<-subset(model, select=-pLeaf)
#C3 <- fx_single_cond(model,32,meta)

modelName = "Dryland_m_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
#rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
rf_df <- readRDS(fileName)
#p6 <- fx_plot(rf_df,modelName)
d6_P <- rf_df

# All models
all_d <- rbind(d1,d2,d3,d4,d5,d5h,d6) %>%
mutate(condition = recode(condition, "1 80" = "Mono.-Meta."
			  ,"1 180" = "Mono.-Iso."
			  ,"32 80" = "Mix.-Meta."
			  ,"32 180" = "Mix.-Iso.")
       ,mName = recode(mName, "Grass1_m" = "Grass 1"
		     ,"Grass2_m" = "Grass 2"
		     ,"Grass3_m" = "Grass 3"
		     ,"Forest1_m" = "Forest 1"
		     ,"Forest2_m" = "Forest 2"
		     ,"Forest2_hrealmax_m" = "Forest 2h"
		     ,"Dryland_m" = "Dryland")
)

ggplot(all_d, aes(x=reorder(varnames, CPI), y=sCPI, fill=type)) +
#ggplot(all_d, aes(x=reorder(varnames, var_categ), y=sCPI, fill=type)) +
#ggplot(all_d, aes(x=varnames, y=sCPI, fill=type)) +
    geom_bar(position='dodge',stat='identity') +
#            geom_text(aes(label=scientific(CPI, digits = 2),size=10)
#    		  ,position = position_dodge(width = 1)
#                      ,hjust=0,vjust=0.2,color="black",size=2
#                      ,show.legend = FALSE,angle = 90) +
    ylab("Conditional permutation importance, scaled to correlation") +
    xlab("Traits") +
    scale_fill_brewer(palette = "Set1",name = "Trait type") +
    scale_color_brewer(palette = "Set1") +
    ylim(NA,1) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5)) + 
#    ylim(NA,1.25) +
    facet_grid(condition ~ mName, scales = "free_x", 
     space = "free_x") + # Let the width of facets vary and force all bars to have the same width.
    theme_bw() +
    theme(text = element_text(size = 14),legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file=paste0(tmp_dir,"/randomForest/all_models_cforest_grid.pdf")
       , width=9.5, height=7, dpi=300
)
while (!is.null(dev.list()))  dev.off()

#Productivity
all_d_P <- rbind(d1_P,d2_P,d3_P,d4_P,d5_P,d5h_P,d6_P) %>%
mutate(condition = recode(condition, "1 80" = "Mono.-Meta."
			  ,"1 180" = "Mono.-Iso."
			  ,"32 80" = "Mix.-Meta."
			  ,"32 180" = "Mix.-Iso.")
       ,mName = recode(mName, "Grass1_m_P" = "Grass 1"
		     ,"Grass2_m_P" = "Grass 2"
		     ,"Grass3_m_P" = "Grass 3"
		     ,"Forest1_m_P" = "Forest 1"
		     ,"Forest2_m_P" = "Forest 2"
		     ,"Forest2_hrealmax_m_P" = "Forest 2h"
		     ,"Dryland_m_P" = "Dryland")
)

ggplot(all_d_P, aes(x=reorder(varnames, CPI), y=sCPI, fill=type)) +
    geom_bar(position='dodge',stat='identity') +
    ylab("Conditional permutation importance, scaled to correlation") +
    xlab("Traits") +
    scale_fill_brewer(palette = "Set1",name = "Trait type") +
    scale_color_brewer(palette = "Set1") +
    ylim(NA,1) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5)) + 
    facet_grid(condition ~ mName, scales = "free_x", 
     space = "free_x") + # Let the width of facets vary and force all bars to have the same width.
    theme_bw() +
    theme(text = element_text(size = 14),legend.position = "top",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file=paste0(tmp_dir,"/randomForest/all_models_cforest_grid_P.pdf")
       , width=9.5, height=7, dpi=300
)
while (!is.null(dev.list()))  dev.off()
