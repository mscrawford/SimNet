library(party)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(ggpubr)
library(ggrepel)

set.seed(1987)

meta <- seq(80, 100)
iso <- seq(180, 200)

#Rcolorbrewer colors:
red <- "#e41a1cff"
blue <- "#377eb8ff"
purple <- "#984ea3ff"

fx_cforest_single_condition <- function(modelName,model,NoSpp,stage){
# Function to perform random forest analysis in a single condition (e.g. Mixture-Metacommuity)
	# Prepare data before cforest analysis
	is_productivity = grepl("_P$",modelName)
	#print("Single cond")
	sizeT <- c("abmi","Vi","MaxMass","PC2score","dmax","ah","hmax","h_realmax","maxSize")
        model <- model %>%
            filter(Ninitial == NoSpp,
                   Year %in% stage) %>%
            select(-SpeciesID, -Ninitial, -Stage, -Rep, -Year)
	model <- subset (model, select = -(if(is_productivity){Biomass}else{Productivity})) 
        
	# Cforest analysis
	set.seed(1987)
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
		      
	#print(rf)
        pred <- data.frame(pred = predict(rf, newdata=test,OOB=TRUE))
	if(is_productivity){test_res = test$Productivity}else{test_res = test$Biomass}
	corr <- cor(pred, test_res)[[1]] #, method = "kendall" #pearson (default)
	var <- round((corr^2), 2) #, method = "kendall" #pearson (default)
        title = paste("Correlation: ", corr,sep = "")

        ##Conditional permutation importance:
	# We used conditional = TRUE if predictor variables are correlated
        CPI <- varimp(rf, conditional = TRUE)
        #CPI <- varimp(rf, conditional = FALSE)

        #Prepare data frame for plotting
        rf_df <- as.data.frame(CPI)
        rf_df$varnames <- rownames(rf_df)
        rownames(rf_df) <- NULL
	if(grepl("^Grass3_PCA",modelName)){
	rf_df$type <- ifelse(rf_df$varnames == "PC2score","Mixed","Resource related")}
	else{rf_df$type <- ifelse(rf_df$varnames %in% sizeT,"Size related","Resource related")}
        rf_df <- arrange(rf_df,desc(CPI))
        rf_df$var_categ <- c(1: dim(rf_df)[1])
	cond <- paste(NoSpp,stage[1])
	rf_df$condition <- rep(cond,times=dim(rf_df)[1])
	rf_df <- rf_df %>% 
		mutate(CPI = replace(CPI, which(CPI<0), NA) ) %>%
		mutate(sCPI = (CPI*var/sum(CPI,na.rm=TRUE)))
		#mutate(sCPI = (CPI*corr/sum(CPI,na.rm=TRUE)))
        #print(rf_df)
	#print(paste0("var=",var,"; sum of scaled importance=",sum(rf_df$CPI,na.rm=TRUE)))
	#print(paste0("corr=",corr,"; sum of scaled importance=",sum(rf_df$CPI,na.rm=TRUE)))
        return(rf_df)
}

fx_cforest <- function(model,modelName,fileName){
# Function to perform random forest for all four conditions (for a single model)
	#print(paste(modelName,"_fx_cforest"))
	if(SAVE_CACHE){
	C1 <- fx_cforest_single_condition(modelName,model,1,meta)
	C2 <- fx_cforest_single_condition(modelName,model,1,iso)
	C3 <- fx_cforest_single_condition(modelName,model,32,meta)
	C4 <- fx_cforest_single_condition(modelName,model,32,iso)
        rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
	}
	else{rf_df <- readRDS(fileName)}
	plot <- fx_plot(rf_df,modelName)
        return(rf_df)
}

fx_cforest_df <- function(C1,C2,C3,C4,modelName){
# Function to make a data frame out of the random forest results for all four conditions (for a single model)
	#print(paste(modelName,"_4conditions_df"))
	rf_df <- rbind(C1,C2,C3,C4)
	rf_df$mName <- rep(modelName,times=dim(rf_df)[1])
	rf_df <- fx_edit_final_df(rf_df)
	#print(rf_df)
	fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
	#fileName = paste0(tmp_dir,"/randomForest/conditionalFALSE/",modelName,".Rda")
	saveRDS(rf_df,file=fileName)
	return(rf_df)}

fx_plot <- function(rf_df,modelName){
# Function to plot the random forest results for all four conditions (for a single model)
	#print(paste(modelName,"_plot"))
	vals1 = c("Resource related" = red, "Size related" = blue)
	vals2 = c("Resource related" = red,"Mixed" = purple)
	rf_df <- rf_df %>%
		mutate(condition = recode(condition, "Mono.-Meta." = "Monoculture-Metacommunity"
				  ,"Mono.-Iso." = "Monoculture-Isolation"
				  ,"Mix.-Meta." = "Mixture-Metacommunity"
				  ,"Mix.-Iso." = "Mixture-Isolation"))        
	
	llim <- (min(rf_df$sCPI))
        llim <- llim-(llim*.05)
        hlim <- (max(rf_df$sCPI))
        hlim <- hlim+(hlim*.1)
	
	resvar = if(grepl("_P$",modelName)){"Productivity"}else{"Biomass"} 
	plot <- ggplot(rf_df, aes(x=reorder(varnames, typen), y=sCPI, fill=type)) +
	    geom_bar(position='dodge',stat='identity') +
	            geom_text(aes(label=scientific(sCPI, digits = 2))
	    		  ,position = position_dodge(width = 1)
	                      ,hjust=0.5,vjust=-0.5,color="black"#,size=7.5
	                      ,show.legend = FALSE) +#,angle = 45
	    #ggtitle(paste0("Trait importance in determining ",resvar)) +
	    ylab("Proportion of variance explained") +
	    #ylab("Conditional permutation importance, scaled") +
	    xlab("Traits") +
	    scale_fill_manual(name = "Trait type",
			      values = if("Mixed" %in% rf_df$type){vals2}else{vals1}) +
	    ylim(NA,hlim) +
	    facet_grid(~condition, 
             space = "free_x") + # Let the width of facets vary and force all bars to have the same width.
	    theme_bw() +
	    theme(text = element_text(size = 20),legend.position = "top",
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
	    plot.title = element_text(hjust = 0.5))
	ggsave(file=paste0(tmp_dir,"/randomForest/",modelName,".pdf")
	#ggsave(file=paste0(tmp_dir,"/randomForest/conditionalFALSE/",modelName,".pdf")
	       , width=15, height=7, dpi=300)
	while (!is.null(dev.list()))  dev.off()
        return(plot)
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
	M <- M %>%
		mutate(id = row_number()) %>%
		mutate_if(is.character, as.factor) %>%
		mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
		mutate(Productivity = scales::rescale(Productivity, to = c(0, 100)))
        return(M)
}	

fx_edit_final_df <- function(df){
	#print("Edited df")
# Function to edit the data frame of the random forest results for all four conditions (for all models)
	df <- df %>%
		mutate(mName = gsub('_P$', '',mName)) %>%
		mutate(condition = recode(condition, "1 80" = "Mono.-Meta."
					  ,"1 180" = "Mono.-Iso."
					  ,"32 80" = "Mix.-Meta."
					  ,"32 180" = "Mix.-Iso.")
	,typen = case_when(type == "Size related"~0,
			   type == "Resource related"~1,
			   type == "Mixed"~0)
	,funcdom = case_when(mName == "Grass1"~0.8,#1,
			     grepl("^Grass2",mName)~0.6,#3, 
			     grepl("^Grass3",mName)~0.5,#4, 
			     mName == "Forest1"~0.65,#2, 
			     grepl("^Forest2",mName)~0.45,#5, 
			     mName == "Dryland"~0.05)#6)
	,funcdom_p = case_when(mName == "Grass1"~0.8,#1,
			     grepl("^Grass2",mName)~0.6,#2, 
			     grepl("^Grass3",mName)~0.5,#3, 
			     mName == "Forest1"~0.25,#5,
			     grepl("^Forest2",mName)~0.35,#4, 
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
			,"Forest2" = "Forest 2"
			,"Forest2_hrealmax" = "Forest 2"
			,"Forest2_hrealmax_PCA" = "Forest 2"
			,"Forest2_PCA" = "Forest 2"
			,"Dryland" = "Dryland"))
	df <- df[complete.cases(df), ] #remove NA
	#print(df)
	return(df)
}

fx_plot_all <- function(df,resvar,plot_name){	
# Function to plot the random forest results for all four conditions (for all models)
	# include function-dominance correlation in model name
	df$mNameFDC <-  paste0(df$mName,'\n (',if(resvar=="Biomass"){df$funcdom}else{df$funcdom_p},')')
	p <- ggplot(df, aes(x=reorder(varnames,typen), y=sCPI, fill=type)) +
	    geom_bar(position='dodge',stat='identity') +
	#            geom_text(aes(label=scientific(CPI, digits = 2),size=10)
	#    		  ,position = position_dodge(width = 1)
	#                      ,hjust=0,vjust=0.2,color="black",size=2
	#                      ,show.legend = FALSE,angle = 90) +
	    #ggtitle(paste0("Trait importance in determining ",resvar)) +
	    ylab("Proportion of variance explained") +
	    #ylab("Conditional permutation importance, \n scaled to r^2") +
	    #ylab("Conditional permutation importance, \n scaled to correlation") +
	    xlab("Traits") +
	    scale_fill_manual(name = "Trait type",
			       values = c("Resource related" = red,
					  "Size related" = blue,
					  "Mixed" = purple)) +
	    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5)) + 
	    facet_grid(reorder(condition,modeln) ~ reorder(mNameFDC, -if(resvar=="Biomass"){funcdom}else{funcdom_p}), scales = "free_x") +
	    theme_bw() +
	    theme(text = element_text(size = 26),legend.position = "right",
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
	    plot.title = element_text(hjust = 0.5))

	ggsave(file=paste0(tmp_dir,"/randomForest/",plot_name,".png")
	       , width=13, height=8, dpi=300
	)
	while (!is.null(dev.list()))  dev.off()
	return(p)
}

fx_plot_all_fdc_plot <- function(df,resvar,plot_name){	
# Function to plot the random forest results for all four conditions (for all models)
	plot_fdc <- ggplot(df, aes(x=reorder(mName, -if(resvar=="Biomass"){funcdom}else{funcdom_p}),
				   y=if(resvar=="Biomass"){funcdom}else{funcdom_p},
				   fill=-if(resvar=="Biomass"){funcdom}else{funcdom_p})) +
		geom_bar(position='dodge',stat='identity') +
		ylab("Func.-\nDom.\nCorr.") +
		xlab("Models") +
		theme_classic() +
		theme(text = element_text(size = 20),legend.position = "none",
		      axis.title.y = element_text(angle = 0))

	p <- ggplot(df, aes(x=reorder(varnames,typen), y=sCPI, fill=type)) +
	    geom_bar(position='dodge',stat='identity') +
	#            geom_text(aes(label=scientific(CPI, digits = 2),size=10)
	#    		  ,position = position_dodge(width = 1)
	#                      ,hjust=0,vjust=0.2,color="black",size=2
	#                      ,show.legend = FALSE,angle = 90) +
	    #ggtitle(paste0("Trait importance in determining ",resvar)) +
	    ylab("Proportion of variance explained") +
	    #ylab("Conditional permutation importance, \n scaled to r^2") +
	    #ylab("Conditional permutation importance, \n scaled to correlation") +
	    xlab("Traits") +
	    scale_fill_manual(name = "Trait type",
			       values = c("Resource related" = red,
					  "Size related" = blue,
					  "Mixed" = purple)) +
	    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.5)) + 
	    facet_grid(reorder(condition,modeln) ~ reorder(mName, -if(resvar=="Biomass"){funcdom}else{funcdom_p}), scales = "free_x") +#, space = "free_x") +  # Let the width of facets vary and force all bars to have the same width.
	    theme_bw() +
	    theme(text = element_text(size = 28),legend.position = "top",
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
	    plot.title = element_text(hjust = 0.5))

	    fileName = paste0(tmp_dir,"/randomForest/",plot_name,".png")
	    #fileName = paste0(tmp_dir,"/randomForest/",plot_name,".pdf")
	    #fileName = paste0(tmp_dir,"/randomForest/conditionalFALSE/",plot_name,".pdf")
	    png(fileName, width=864, height=768, units = "px")
	    #pdf(fileName, width=9, height=8)
	    grid.arrange(plot_fdc, p, nrow=2, heights=c(1,5))#, widths=c(4,1))
	    dev.off()
#	ggsave(file=paste0(tmp_dir,"/randomForest/",plot_name,".pdf")
#	       , width=9.5, height=7, dpi=300
#	)
#	while (!is.null(dev.list()))  dev.off()
	return(p)
}

fx_diff <- function(df){	
# Calculate difference between biomass variance explained by size (and resource related) traits in monoculture vs. mixture. 
	mono <- df[df$condition %in% c('Mono.-Meta.'), ] %>%
		group_by(type, mName, funcdom, funcdom_p) %>%
		summarize(sCPI=sum(sCPI)) 
	mix <- df[df$condition %in% c('Mix.-Meta.'), ] %>%
		group_by(type, mName, funcdom, funcdom_p) %>%
		summarize(sCPI=sum(sCPI)) 
	new_df = mono
	new_df['sCPI'] = mono['sCPI'] - mix['sCPI'] 
	return(new_df)
}

fx_diff_v2 <- function(df){	
  # Calculate difference between biomass variance explained by size (and resource related) traits in monoculture vs. mixture. 
  mono <- df[df$condition %in% c('Mono.-Meta.'), ] %>%
    group_by(type, mName, funcdom, funcdom_p) %>%
    summarize(sCPI=sum(sCPI)) 
  mix <- df[df$condition %in% c('Mix.-Meta.'), ] %>%
    group_by(type, mName, funcdom, funcdom_p) %>%
    summarize(sCPI=sum(sCPI)) 
  new_df = mono
  new_df['sCPI'] = sqrt(mono['sCPI'] * mix['sCPI']) 
  return(new_df)
}

fx_plot_diff_mono_mix <- function(plotName,df){
# Plot product of mono and mixture % of biomass or productivity variance explained per trait type
	print('Difffffff')
	print(head(df))
#	df <- df %>%
#		mutate(sCPI = abs(sCPI))
	print(head(df))
	is_productivity = grepl("_P$",plotName)
	response <- if(is_productivity){"Productivity"}else{"Biomass"}
	df_stat_smooth <- df
	print(df_stat_smooth)

	xlab <- 'Function-dominance correlation'
	ylab <- expression(sqrt("Monoculture VE * Mixture VE"))
  p1 <- ggplot(df,
               aes(x = if(is_productivity){funcdom_p}else{funcdom},
                          y = sCPI, label = mName)) + #size=3 
                          #color = type)) +
	    geom_text_repel(color="black", 
	      guide="none") +
    #guides(size="none", fill="none") + 
    labs(shape = "Model", x = xlab, y = ylab) +
    #stat_regline_equation(label.x = c(0.25,0.55), label.y = c(1,1),aes(label =  ..adj.rr.label..)) +
    stat_regline_equation(data = df_stat_smooth, label.x = 0.1, label.y = 0.9, show.legend=FALSE,
			  aes(label =  ..adj.rr.label..), color=blue, fill=blue) +
    geom_point(color="black", shape=21, fill=blue) +
    scale_color_manual(values=blue) +
    scale_fill_manual(values=blue, guide="none")+
    #geom_hline(yintercept=0, linetype='dotted') +
    theme_bw() +
    theme_classic() +
    theme(text = element_text(size = 14), legend.position = "top", legend.direction = "horizontal") 
  
  filename <- paste0(plotName,".png")
  #filename <- paste0(plotName,".pdf")
  path <- paste0(tmp_dir,"/randomForest/")
  ggsave(filename = filename, path = path, plot = p1
         ,height = 13, width = 15, units = "cm")
  return(p1)
}

fx_plot_diff_mono_mix_smooth <- function(plotName,df){
# Add regression line to fx_plot_diff_mono_mix
#	df <- df %>%
#		mutate(sCPI = abs(sCPI))
	#exclude data from Dryland-Resource related for the regression
	df_stat_smooth <- subset(df,type!="Resource related" | mName!="Dryland")
#	df_stat_smooth <- df #%>%
#		filter(type == 'Size related')
	p1 = fx_plot_diff_mono_mix(plotName,df) + 
		stat_smooth(data = df_stat_smooth, 
		method="lm", se= FALSE, linetype='dashed',  color=blue,
	       	method.args = list(family = "symmetric"))
#    stat_smooth(method="lm", linetype='dashed', alpha = 0.2, aes(fill=type, color=type)) +
#    scale_fill_manual(name = "Trait type",
#		      values = c("Resource related" = red,
#				 "Size related" = blue)) +
    filename <- paste0(plotName,".png")
    path <- paste0(tmp_dir,"/randomForest/")
    ggsave(filename = filename, path = path, plot = p1
           ,height = 13, width = 15, units = "cm")
    return(p1)
}
