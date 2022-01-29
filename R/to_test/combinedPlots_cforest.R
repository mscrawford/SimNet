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
	,funcdom = case_when(mName == "Grass1"~1,
			     mName == "Forest1"~2, 
			     mName == "Grass2"~3, 
			     mName == "Grass3"~4, 
			     mName == "Dryland"~5, 
			     mName == "Forest2_hrealmax"~6)
	,modeln = case_when(condition == "Mono.-Meta."~1,
			   condition == "Mono.-Iso."~2,
			   condition == "Mix.-Meta."~3,
			   condition == "Mix.-Iso."~4)
#	,modeln = gsub('Mono.*',1,condition)
#	,modeln = gsub('Mix.*',2,modeln)
	,mName = recode(mName, "Grass1" = "Grass 1"
			,"Grass2" = "Grass 2"
			,"Grass3" = "Grass 3"
			,"Forest1" = "Forest 1"
			#,"Forest2" = "Forest 2"
			,"Forest2_hrealmax" = "Forest 2"
			,"Dryland" = "Dryland"))
	print(df)
	df <- df[complete.cases(df), ] #remove NA
	print(df)
	return(df)
}

fx_plot_all <- function(df,plot_title,plot_name){	
	p <- ggplot(df, aes(x=reorder(varnames,typen), y=sCPI, fill=type)) +
	    geom_bar(position='dodge',stat='identity') +
	#            geom_text(aes(label=scientific(CPI, digits = 2),size=10)
	#    		  ,position = position_dodge(width = 1)
	#                      ,hjust=0,vjust=0.2,color="black",size=2
	#                      ,show.legend = FALSE,angle = 90) +
	    ggtitle(plot_title) +
	    ylab("Conditional permutation importance, scaled to correlation") +
	    xlab("Traits") +
	    scale_fill_brewer(palette = "Set1",name = "Trait type") +
	    scale_color_brewer(palette = "Set1") +
	    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5)) + 
	    facet_grid(reorder(condition,modeln) ~ reorder(mName, funcdom), scales = "free_x", 
	     space = "free_x") + # Let the width of facets vary and force all bars to have the same width.
	    theme_bw() +
	    theme(text = element_text(size = 14),legend.position = "top",
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
	    plot.title = element_text(hjust = 0.5))
	ggsave(file=paste0(tmp_dir,"/randomForest/",plot_name,".pdf")
	       , width=9.5, height=7, dpi=300
	)
	while (!is.null(dev.list()))  dev.off()
	return(p)
}

### Grass1 (Adam's model)
if(switch==1){adam <- fx_read_model("readAdam.R", "Grass1")}else{}

#Biomass
if(switch==1){model<-adam[["B"]]}else{}
modelName = "Grass1"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d1 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-adam[["P"]]}else{}
modelName = "Grass1_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d1_P <- fx_run_cforest(model,modelName,fileName)

### Grass2 (Lindsay's model)
if(switch==1){lindsay <- fx_read_model("readLindsay.R","Grass2")}else{}

#Biomass
if(switch==1){model<-lindsay[["B"]]}else{}
modelName = "Grass2"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-lindsay[["P"]]}else{}
modelName = "Grass2_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d2_P <- fx_run_cforest(model,modelName,fileName)

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

##See what the data frame looks like for only monoculture-isolation:
#t2 <- model[which(model$Stage == 'Without seed inflow' & model$Ninitial == 1), ]
#mean(t2$Biomass)
#sd(t2$Biomass)
##Compare with respect to Mixture-isolation
#t3 <- model[which(model$Stage == 'Without seed inflow' & model$Ninitial == 32), ]
#mean(t3$Biomass)
#sd(t3$Biomass)
modelName = "Forest2_hrealmax_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d5h_P <- fx_run_cforest(model,modelName,fileName)

### Dryland (Bjoern)
if(switch==1){bjoern <- fx_read_model("readBjoern.R","bjoern")}else{}

#Biomass
if(switch==1){model<-bjoern[["B"]]}else{}
modelName = "Dryland"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d6 <- fx_run_cforest(model,modelName,fileName)

#Productivity
if(switch==1){model<-bjoern[["P"]]}else{}
modelName = "Dryland_P"
fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
d6_P <- fx_run_cforest(model,modelName,fileName)

# All models
#Biomass
all_d <- rbind(d1,d2,d3,d4,d5h,d6)
all_d <- fx_edit_final_df(all_d)
fx_plot_all(all_d,"Trait importance in determining Biomass","all_models_cforest_grid_c")
all_d <- all_d[all_d$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d,"Trait importance in determining Biomass","all_models_cforest_grid_c_meta")

#Productivity
all_d_P <- rbind(d1_P,d2_P,d3_P,d4_P,d5h_P,d6_P)
all_d_P <- fx_edit_final_df(all_d_P)
fx_plot_all(all_d_P,"Trait importance in determining Productivity","all_models_cforest_grid_c_P")
all_d_P <- all_d_P[all_d_P$condition %in% c('Mono.-Meta.','Mix.-Meta.'), ]
fx_plot_all(all_d_P,"Trait importance in determining Productivity","all_models_cforest_grid_c_P_meta")
