library(party)
library(ggplot2)
library(dplyr)
library(scales)
set.seed(1987)

meta <- seq(80, 100)
iso <- seq(180, 200)

fx_single_cond <- function(modelName,model,NoSpp,stage){
# Function to perform random forest analysis in a single condition (e.g. Mixture-Metacommuity)
	is_productivity = grepl("_P$",modelName)
	print("Single cond")
	sizeT <- c("abmi","Vi","MaxMass","PC2score","dmax","ah","hmax","h_realmax","maxSize")
        model <- model %>%
            filter(Ninitial == NoSpp,
                   Year %in% stage) %>%
            mutate(id = row_number()) %>%
            mutate_if(is.character, as.factor) %>%
	    mutate(Biomass = scales::rescale(Biomass, to = c(0, 100)))%>%
	    mutate(Productivity = scales::rescale(Productivity, to = c(0, 100)))%>%
            select(-SpeciesID, -Ninitial, -Stage, -Rep, -Year)
	model <- subset (model, select = -(if(is_productivity){Biomass}else{Productivity})) 
        
        train <- model %>% sample_frac(.70)
        test <- anti_join(model, train, by = 'id')
        
        train <- train %>% select(-id)
        test <- test %>% select(-id)
        
        rf <- if(is_productivity){cforest(Productivity ~ .,
                      data = as.data.frame(train),
		      control = cforest_unbiased(mtry = 2, ntree = 501))}
		else{cforest(Biomass ~ .,
                      data = as.data.frame(train),
		      control = cforest_unbiased(mtry = 2, ntree = 501))}
		      
	print(rf)
        pred <- data.frame(pred = predict(rf, newdata=test,OOB=TRUE))
	if(is_productivity){test_res = test$Productivity}else{test_res = test$Biomass}
	corr <- round(cor(pred, test_res)[[1]], 2) #, method = "kendall" #pearson (default)
        title = paste("Correlation: ", corr,sep = "")

        #Conditional permutation importance:
        CPI <- varimp(rf, conditional = TRUE)
        #Prepare data frame for plotting
        rf_df <- as.data.frame(CPI)
        rf_df$varnames <- rownames(rf_df)
        rownames(rf_df) <- NULL
	if(grepl("^Grass3_PCA",modelName)){
	rf_df$type <- ifelse(rf_df$varnames == "PC2score","Size/Resource related","Resource related")}
	else{rf_df$type <- ifelse(rf_df$varnames %in% sizeT,"Size related","Resource related")}
        rf_df <- arrange(rf_df,desc(CPI))
        rf_df$var_categ <- c(1: dim(rf_df)[1])
	cond <- paste(NoSpp,stage[1])
	rf_df$condition <- rep(cond,times=dim(rf_df)[1])
	rf_df <- rf_df %>% 
		mutate(CPI = replace(CPI, which(CPI<0), NA) ) %>%
		mutate(sCPI = (CPI*corr/sum(CPI,na.rm=TRUE)))
        print(rf_df)
	print(paste0("corr=",corr,"; sum of scaled importance=",sum(rf_df$CPI,na.rm=TRUE)))
        return(rf_df)
}

fx_cforest <- function(model,modelName){
# Function to perform random forest for all four conditions (for a single model)
	print(paste(modelName,"_fx_cforest"))
	C1 <- fx_single_cond(modelName,model,1,meta)
	C2 <- fx_single_cond(modelName,model,1,iso)
	C3 <- fx_single_cond(modelName,model,32,meta)
	C4 <- fx_single_cond(modelName,model,32,iso)
        rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
	plot <- fx_plot(rf_df,modelName)
        return(rf_df)
}

fx_cforest_df <- function(C1,C2,C3,C4,modelName){
# Function to make a data frame out of the random forest results for all four conditions (for a single model)
	print(paste(modelName,"_4conditions_df"))
	rf_df <- rbind(C1,C2,C3,C4)
	rf_df$mName <- rep(modelName,times=dim(rf_df)[1])
	fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
	saveRDS(rf_df,file=fileName)
	return(rf_df)}

fx_plot <- function(rf_df,modelName){
# Function to plot the random forest results for all four conditions (for a single model)
	print(paste(modelName,"_plot"))
	rf_df <- rf_df %>%
	mutate(condition = recode(condition, "1 80" = "Monoculture-Metacommunity"
				  ,"1 180" = "Monoculture-Isolation"
				  ,"32 80" = "Mixture-Metacommunity"
				  ,"32 180" = "Mixture-Isolation"))        
	
	llim <- (min(rf_df$sCPI))
        llim <- llim-(llim*.05)
        hlim <- (max(rf_df$sCPI))
        hlim <- hlim+(hlim*.1)
	
	resvar = if(grepl("_P$",modelName)){"Productivity"}else{"Biomass"} 
	plot <- ggplot(rf_df, aes(x=reorder(varnames, var_categ), y=sCPI, fill=type)) +
	    geom_bar(position='dodge',stat='identity') +
	            geom_text(aes(label=scientific(CPI, digits = 2))
	    		  ,position = position_dodge(width = 1)
	                      ,hjust=0.5,vjust=-0.5,color="black"#,size=7.5
	                      ,show.legend = FALSE) +#,angle = 45
	    ggtitle(paste0("Trait importance in determining ",resvar)) +
	    ylab("Conditional permutation importance, scaled") +
	    xlab("Traits") +
	    scale_fill_manual(name = "Trait type",
			       values = c("Resource related" = "red",
					  "Size related" = "blue",
					  "Size/Resource related" = "purple")) +
	    ylim(NA,hlim) +
	    facet_grid(~condition, 
             space = "free_x") + # Let the width of facets vary and force all bars to have the same width.
	    theme_bw() +
	    theme(text = element_text(size = 20),legend.position = "top",
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
	    plot.title = element_text(hjust = 0.5))
	ggsave(file=paste0(tmp_dir,"/randomForest/",modelName,".pdf")
	       , width=15, height=7, dpi=300)
	while (!is.null(dev.list()))  dev.off()
        return(plot)
}

fx_run_cforest <- function(model,modelName,fileName){
	if(SAVE_CACHE){rf_df <- fx_cforest(model,modelName)}
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
        return(M)
}	

fx_edit_final_df <- function(df){
# Function to edit the data frame of the random forest results for all four conditions (for all models)
	df <- df %>%
		mutate(mName = gsub('_P$', '',mName)) %>%
		mutate(condition = recode(condition, "1 80" = "Mono.-Meta."
					  ,"1 180" = "Mono.-Iso."
					  ,"32 80" = "Mix.-Meta."
					  ,"32 180" = "Mix.-Iso.")
	,typen = case_when(type == "Size related"~0,
			   type == "Resource related"~1,
			   type == "Size/Resource related"~0)
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
# Function to plot the random forest results for all four conditions (for all models)
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
	    scale_fill_manual(name = "Trait type",
			       values = c("Resource related" = "red",
					  "Size related" = "blue",
					  "Size/Resource related" = "purple")) +
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
