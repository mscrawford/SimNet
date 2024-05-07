library(ggplot2)
library(tidyr)
library(party)
#library(scales)
#library(dplyr)

set.seed(1987)

fx_plot_trait_Vs_biomass <- function(df, condition, plot_name){
df <- na.omit(df)
df <- df %>%
    filter(if(condition=="Monoculture"){NumSp == 1}else{NumSp > 1}) %>%
    mutate(Year = as.character(Year), log_bm = log(species_biomass_m2)) %>%
    #gather(-Species, -Year, -id, -NumSp, -species_biomass_m2,-log_bm, key= "var", value = "value") %>%
    gather(-Species, -Plot, -Year, -id, -NumSp, -species_biomass_m2,-log_bm, key= "var", value = "value") %>%
    ggplot(aes(x=as.numeric(value), y=log_bm, color=Year, shape=Year)) +
    geom_point() +
    labs(y = "log biomass", x = "Trait value") +
    geom_smooth(method="lm", fill=NA) +
    #stat_smooth(fullrange = TRUE, color="red", method="loess", se=FALSE) +
    facet_wrap(~var, scales ="free") +
    scale_x_continuous(n.breaks = 4) +
    #scale_x_continuous(breaks = equal_breaks(n = 4, s=0.05)) +
    theme_bw() +
    theme(text = element_text(size = 25), strip.text = element_text(size=25),
          legend.text = element_text(size=20), legend.title = element_text(size=20))
    ggsave(file=paste0(tmp_dir, plot_name)
           , width=18*2, height=14*2, dpi=300)
    while (!is.null(dev.list()))  dev.off()
    return()
}

fx_cforest_data_sets <- function(model,condition){
#Adapted from fx_cforest_single_condition <- function(modelName,model,condition,stage){
sizeT <- c("height_.m.", "shoot.height", "RootingDepth_Target", "Height", "HEIGHT","MAXHEIGHT_AVG")
print(typeof(model))
	model <- model %>%
            filter(if(condition=="Monoculture"){NumSp == 1}else{NumSp > 1}) %>%
            #select(-one_of("Species", "Year", "NumSp"))
            select(-one_of("Species", "Plot", "Year", "NumSp"))
print("Model cforest")
print(head(model))
        
	# Cforest analysis
	set.seed(1987)
#        train <- model %>% sample_frac(.70)
#        test <- anti_join(model, train, by = 'id')
        
#        train <- train %>% select(-id)
#        test <- test %>% select(-id)
        model <- model %>% select(-id)
        
        rf <- cforest(species_biomass_m2 ~ .,
                      data = as.data.frame(model),
                      #data = as.data.frame(train),
		      #control = cforest_unbiased(mtry = 1, ntree = 2000))
		      control = cforest_unbiased(mtry = 2, ntree = 501))
		      #control = cforest_unbiased(mtry = 2, ntree = 1001))
		     
	print("%%%%%%%%%%%%%%%%%%%%% summary(rf) %%%%%%%%%%%%%%%%%%%%%%%%%")
	summary(rf)
        pred <- data.frame(pred = predict(rf, OOB=TRUE))
	test_res = model$species_biomass_m2
	corr <- cor(pred, test_res)[[1]] #, method = "kendall" #pearson (default)
	var <- round((corr^2), 2) #, method = "kendall" #pearson (default)
        title = paste("Correlation: ", corr,sep = "")

        ##Conditional permutation importance:
	# We used conditional = TRUE if predictor variables are correlated
        CPI <- varimp(rf, conditional = TRUE)
        #CPI <- varimp(rf, conditional = FALSE)
	print(head(CPI))

        #Prepare data frame for plotting
        rf_df <- as.data.frame(CPI)
        rf_df$varnames <- rownames(rf_df)
        rownames(rf_df) <- NULL
	rf_df$type <- ifelse(rf_df$varnames %in% sizeT,"Size related","Resource related")
        rf_df <- arrange(rf_df,desc(CPI))
        rf_df$var_categ <- c(1: dim(rf_df)[1])
	cond <- NoSpp
	rf_df$condition <- rep(cond,times=dim(rf_df)[1])
	rf_df <- rf_df %>% 
		mutate(CPI = replace(CPI, which(CPI<0), NA) ) %>%
		#mutate(CPI = replace(CPI, which(CPI<0), 0) ) %>%
		#mutate(sCPI = (CPI/sum(CPI,na.rm=TRUE)))
	  mutate(Acc = rep(var,times=dim(rf_df)[1]) ) %>%
		mutate(sCPI = (CPI*var/sum(CPI,na.rm=TRUE)))
		#mutate(sCPI = (CPI*corr/sum(CPI,na.rm=TRUE)))
        print(rf_df)
	#print(paste0("var=",var,"; sum of scaled importance=",sum(rf_df$CPI,na.rm=TRUE)))
	#print(paste0("corr=",corr,"; sum of scaled importance=",sum(rf_df$CPI,na.rm=TRUE)))
        return(rf_df)
}
