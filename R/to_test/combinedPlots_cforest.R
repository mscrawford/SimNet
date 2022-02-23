library(gridExtra)

set.seed(1987)
#base_dir          <- setwd("~/Documents/SimNet")
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

#source(paste0(scripts_dir, "/readModels.R"))
source(paste0(scripts_dir, "/to_test/fx_cforest_party.R"))

switch = 0
if(switch==0){model <- c()}else{}

fx_run_cforest <- function(model,modelName,fileName){
	if(switch==1){rf_df <- fx_cforest(model,modelName)}
	else{rf_df <- readRDS(fileName)}
        return(rf_df)
}	

fx_read_model <- function(raw_file,models_id){
	source(paste0(scripts_dir, "/to_test/",raw_file), local = TRUE)
	Ms <- as.data.frame(models)
	colnames(Ms) <- gsub(paste0(models_id,"."),"",colnames(Ms))
	select_1 <- 
	if(models_id == "Grass3"){
		M <- Ms %>%
			select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)}
	else{
		M <- Ms %>%
			select(-Model, -SeedRain)}
	Ml <- list()
	Ml[["B"]] <- M %>% select(-Productivity)
	Ml[["P"]] <- M %>% select(-Biomass)
	setnames(Ml[["P"]], "Productivity", "Biomass")
        return(Ml)
}	

fx_edit_final_df <- function(df){
	df <- df %>%
		mutate(mName = gsub('_P$', '',mName)) %>%
		mutate(condition = recode(condition, "1 80" = "Mono.-Meta."
					  ,"1 180" = "Mono.-Iso."
					  ,"32 80" = "Mix.-Meta."
					  ,"32 180" = "Mix.-Iso.")
	,typen = case_when(type == "Size"~0,
			   type == "Resource acq."~1)
	,type = recode(type, "Size"="Size related"
		       ,"Resource acq."="Resource related")
	,funcdom = case_when(mName == "Grass1"~0.8,#1,
			     mName == "Grass2"~0.6,#3, 
			     mName == "Grass2_v"~0.6,#3, 
			     mName == "Grass3"~0.5,#4, 
			     mName == "Grass3_PCA2"~0.5,#4, 
			     mName == "Grass3_PCA3"~0.5,#4, 
			     mName == "Forest1"~0.65,#2, 
			     mName == "Forest2_hrealmax"~0.45,#5, 
			     mName == "Forest2_hrealmax_PCA"~0.45,#5, 
			     mName == "Dryland"~0.05)#6)
	,funcdom_p = case_when(mName == "Grass1"~0.8,#1,
			     mName == "Grass2"~0.6,#2,
			     mName == "Grass2_v"~0.6,#2,
			     mName == "Grass3"~0.5,#3,
			     mName == "Grass3_PCA2"~0.5,#3,
			     mName == "Grass3_PCA3"~0.5,#3,
			     mName == "Forest1"~0.25,#5,
			     mName == "Forest2_hrealmax"~0.35,#4,
			     mName == "Forest2_hrealmax_PCA"~0.35,#4,
			     mName == "Dryland"~0.1)#6)
	,modeln = case_when(condition == "Mono.-Meta."~1,
			   condition == "Mono.-Iso."~2,
			   condition == "Mix.-Meta."~3,
			   condition == "Mix.-Iso."~4)
	,mName = recode(mName, "Grass1" = "Grass 1"
			,"Grass2" = "Grass 2"
			,"Grass2_v" = "Grass 2"
			,"Grass3" = "Grass 3"
			,"Grass3_PCA2" = "Grass 3"
			,"Grass3_PCA3" = "Grass 3"
			,"Forest1" = "Forest 1"
			#,"Forest2" = "Forest 2"
			,"Forest2_hrealmax" = "Forest 2"
			,"Forest2_hrealmax_PCA" = "Forest 2"
			,"Dryland" = "Dryland"))
	df <- df[complete.cases(df), ] #remove NA
	return(df)
}

fx_plot_all <- function(df,resvar,plot_name){	
	plot_fdc <- ggplot(df, aes(x=reorder(mName, if(resvar=="Biomass"){funcdom}else{funcdom_p}),
				   y=if(resvar=="Biomass"){funcdom}else{funcdom_p})) +
		geom_bar(position='dodge',stat='identity')
#	ggsave(file=paste0(tmp_dir,"/randomForest/",plot_name,"_top.pdf")
#	       , width=9.5, height=7, dpi=300
#	)
#	while (!is.null(dev.list()))  dev.off()

	p <- ggplot(df, aes(x=reorder(varnames,typen), y=sCPI, fill=type)) +
	    geom_bar(position='dodge',stat='identity') +
	#            geom_text(aes(label=scientific(CPI, digits = 2),size=10)
	#    		  ,position = position_dodge(width = 1)
	#                      ,hjust=0,vjust=0.2,color="black",size=2
	#                      ,show.legend = FALSE,angle = 90) +
	    ggtitle(paste0("Trait importance in determining ",resvar)) +
	    ylab("Conditional permutation importance, scaled to correlation") +
	    xlab("Traits") +
	    scale_fill_brewer(palette = "Set1",name = "Trait type") +
	    scale_color_brewer(palette = "Set1") +
	    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5)) + 
	    facet_grid(reorder(condition,modeln) ~ reorder(mName, if(resvar=="Biomass"){funcdom}else{funcdom_p}), scales = "free_x", 
	     space = "free_x") +  # Let the width of facets vary and force all bars to have the same width.
	    theme_bw() +
	    theme(text = element_text(size = 14),legend.position = "top",
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
	    plot.title = element_text(hjust = 0.5))

#	    grid.arrange(plot_fdc, empty, p, ncol=1, nrow=2, widths=c(4,1), heights=c(1,4))

	ggsave(file=paste0(tmp_dir,"/randomForest/",plot_name,".pdf")
	       , width=9.5, height=7, dpi=300
	)
	while (!is.null(dev.list()))  dev.off()
	return(p)
}
### Grass2 (Lindsay's model) - variable
if(switch==1){lindsay <- fx_read_model("readLindsay_variable.R","Grass2")}else{}

#Biomass
if(switch==1){model<-lindsay[["B"]] #%>%
	#select(-thetai)
}else{}
modelName = "Grass2_v"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_v <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-lindsay[["P"]] #%>%
#	mutate(Vi = thetai*0)
#	select(-thetai)
}else{}
#modelName = "Grass2_P_t"
modelName = "Grass2_v_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_v_P <- fx_run_cforest(model,modelName,fileName)

### Grass2 (Lindsay's model)
if(switch==1){lindsay <- fx_read_model("readLindsay.R","Grass2")}else{}

#Biomass
if(switch==1){model<-lindsay[["B"]] #%>%
	#select(-thetai)
}else{}
modelName = "Grass2"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-lindsay[["P"]] #%>%
#	mutate(Vi = thetai*0)
#	select(-thetai)
}else{}
#modelName = "Grass2_P_t"
modelName = "Grass2_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_P <- fx_run_cforest(model,modelName,fileName)

switch = 0

### Grass1 (Adam's model)
if(switch==1){adam <- fx_read_model("readAdam.R", "Grass1")}else{}

#Biomass
if(switch==1){model<-adam[["B"]] %>%
	mutate(pNi = 1/pNi) %>%
	select(-no3i,-pNi)}else{}
modelName = "Grass1"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d1 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-adam[["P"]] %>%
	mutate(pNi = 1/pNi) %>%
	select(-no3i,-pNi)}else{}
modelName = "Grass1_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d1_P <- fx_run_cforest(model,modelName,fileName)

### Grass3 (IBC-grass)
if(switch==1){IBC_grass <- fx_read_model("readIBC.R","Grass3")}else{}
#if(switch==1){IBC_grass <- fx_read_model("readIBC.R","IBC_grass.noNDD")}else{}

#Biomass
if(switch==1){model<-IBC_grass[["B"]]}else{}
modelName = "Grass3"
#modelName = "Grass3_noNDD"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-IBC_grass[["P"]]}else{}
modelName = "Grass3_P"
#modelName = "Grass3_P_noNDD"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3_P <- fx_run_cforest(model,modelName,fileName)

### Grass3 (IBC-grass) PCA - 3 components
if(switch==1){grass3 <- readRDS(paste0(tmp_dir,"/PCA/Grass3_PCAcoord.Rda")) %>%
	select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, PC1score, PC2score, PC3score)
Ml <- list()
	Ml[["B"]] <- grass3 %>% select(-Productivity)
	Ml[["P"]] <- grass3 %>% select(-Biomass)
	data.table::setnames(Ml[["P"]], "Productivity", "Biomass")}else{} 

#Biomass
if(switch==1){model<-Ml[["B"]]}else{}
modelName = "Grass3_PCA3"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3_3 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-Ml[["P"]]}else{}
modelName = "Grass3_PCA3_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d3_3_P <- fx_run_cforest(model,modelName,fileName)

### Forest1 (PPA)
if(switch==1){PPA <- fx_read_model("readPPA.R","Forest1")}else{}

#Biomass
if(switch==1){model<-PPA[["B"]]}else{}
modelName = "Forest1"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d4 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-PPA[["P"]]}else{}
modelName = "Forest1_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d4_P <- fx_run_cforest(model,modelName,fileName)

### Forest2 (TROLL)
if(switch==1){troll <- fx_read_model("readTROLL.R","Forest2")}else{}

#Biomass
if(switch==1){model<-troll[["B"]]}else{}
modelName = "Forest2"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-troll[["P"]]}else{}
modelName = "Forest2_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5_P <- fx_run_cforest(model,modelName,fileName)

### Forest2 (TROLL) h_realmax
if(switch==1){troll_hrealmax <- fx_read_model("readTROLL.R","Forest2")}else{}

#Biomass
if(switch==1){model<-troll_hrealmax[["B"]] %>%
   mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
   select(-hmax, -ah, -dmax)}else{}

modelName = "Forest2_hrealmax"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-troll_hrealmax[["P"]] %>%
   mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
   select(-hmax, -ah, -dmax)}else{} 

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
if(switch==1){forest2 <- readRDS(paste0(tmp_dir,"/PCA/Forest2_PCAcoord.Rda"))
Ml <- list()
	Ml[["B"]] <- forest2 %>% select(c(-lma, -nmass, -pmass, -wsg, -h_realmax,-Productivity))
	Ml[["P"]] <- forest2 %>% select(c(-lma, -nmass, -pmass, -wsg, -h_realmax,-Biomass))
	data.table::setnames(Ml[["P"]], "Productivity", "Biomass")}else{} 

#Biomass
if(switch==1){model<-Ml[["B"]]}else{}
modelName = "Forest2_hrealmax_PCA"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_2 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-Ml[["P"]]}else{}
modelName = "Forest2_hrealmax_PCA_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_2_P <- fx_run_cforest(model,modelName,fileName)

### Dryland (Bjoern)
if(switch==1){bjoern <- fx_read_model("readBjoern.R","bjoern")}else{}

#Biomass
if(switch==1){model<-bjoern[["B"]] %>%
select(-pRoot)}else{}
modelName = "Dryland"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d6 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-bjoern[["P"]] %>%
select(-pRoot)}else{}
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
all_d <- rbind(d1,d2,d3_3,d4,d5h_2,d6)
all_d <- fx_edit_final_df(all_d)
fx_plot_all(all_d,resvar,pN5)
all_d <- all_d[all_d$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d,resvar,pN6)


## PCA 3 components for Grass 3, Grass 2 variable
all_d <- rbind(d1,d2_v,d3_3,d4,d5h_2,d6)
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
all_d_P <- rbind(d1_P,d2_P,d3_3_P,d4_P,d5h_2_P,d6_P)
all_d_P <- fx_edit_final_df(all_d_P)
fx_plot_all(all_d_P,resvar,paste0(pN5,"_P"))
all_d_P <- all_d_P[all_d_P$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d_P,resvar,paste0(pN6,"_P"))

## PCA 3 components for Grass 3, Grass 2 variable
all_d_P <- rbind(d1_P,d2_v_P,d3_3_P,d4_P,d5h_2_P,d6_P)
all_d_P <- fx_edit_final_df(all_d_P)
fx_plot_all(all_d_P,resvar,paste0(pN5,"_v_P"))
all_d_P <- all_d_P[all_d_P$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d_P,resvar,paste0(pN6,"_v_P"))
