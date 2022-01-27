library(party)
library(ggplot2)
library(dplyr)
library(scales)
library(RColorBrewer)
set.seed(1987)

meta <- seq(80, 100)
iso <- seq(180, 200)

fx_single_cond <- function(model,NoSpp,stage){
	print("Single cond")
	sizeT <- c("abmi","Vi","MaxMass","PC2score","dmax","ah","hmax","h_realmax","maxSize")
        model <- model %>%
            filter(Ninitial == NoSpp,
                   Year %in% stage) %>%
            mutate(id = row_number()) %>%
            mutate_if(is.character, as.factor) %>%
	    mutate(Biomass = scales::rescale(Biomass, to = c(0, 100)))%>%
            select(-SpeciesID, -Ninitial, -Stage, -Rep, -Year) 
        
        set.seed(1987)
        train <- model %>% sample_frac(.70)
        test <- anti_join(model, train, by = 'id')
        
        train <- train %>% select(-id)
        test <- test %>% select(-id)
        
        rf <- cforest(Biomass ~ .
                      ,data = as.data.frame(train)
		      ,control = cforest_unbiased(mtry = 2, ntree = 501)
		      )
        pred <- data.frame(pred = predict(rf, newdata=test,OOB=TRUE))
        #print(pred)
        #print(test$Biomass)
	corr <- round(cor(pred, test$Biomass)[[1]], 2) #, method = "kendall" #pearson (default)
        title = paste("Correlation: ", corr,sep = "")
        print(rf)
        # print(pred)
        set.seed(1987)
        #Conditional permutation importance:
        set.seed(1987)
        CPI <- varimp(rf, conditional = TRUE)
        #Prepare data frame for plotting
        rf_df <- as.data.frame(CPI)
        rf_df$varnames <- rownames(rf_df)
        rownames(rf_df) <- NULL
	rf_df$type <- ifelse(rf_df$varnames %in% sizeT,"Size","Resource acq.")
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

fx_cforest_df <- function(C1,C2,C3,C4,modelName){
	print(paste(modelName,"_4conditions_df"))
	rf_df <- rbind(C1,C2,C3,C4)
	rf_df$mName <- rep(modelName,times=dim(rf_df)[1])
	fileName = paste0(tmp_dir,"/randomForest/",modelName,".Rda")
	saveRDS(rf_df,file=fileName)
	return(rf_df)}

fx_plot <- function(rf_df,modelName){
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

	plot <- ggplot(rf_df, aes(x=reorder(varnames, var_categ), y=sCPI, fill=type)) +
	#plot <- ggplot(rf_df, aes(x=condition, y=sCPI, fill=reorder(varnames, var_categ))) +
	    #ggplot(rf_df, aes(fill=condition, y=sCPI, x=reorder(varnames, sCPI))) +
	    geom_bar(position='dodge',stat='identity') +
	            geom_text(aes(label=scientific(CPI, digits = 2))
	    		  ,position = position_dodge(width = 1)
	                      ,hjust=0,vjust=0.2,color="black"#,size=7.5
	                      ,show.legend = FALSE,angle = 45) +
	    ylab("Conditional permutation importance, scaled") +
	    xlab("Traits") +
	    scale_fill_brewer(palette = "Set1",name = "Trait type") +
	    scale_color_brewer(palette = "Set1") +
	    ylim(NA,hlim) +
	    facet_grid(~condition, 
             space = "free_x") + # Let the width of facets vary and force all bars to have the same width.
	    theme_bw() +
	    theme(text = element_text(size = 20),legend.position = "top",
	    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
	ggsave(file=paste0(tmp_dir,"/randomForest/",modelName,".pdf")
	       , width=15, height=7, dpi=300)
	while (!is.null(dev.list()))  dev.off()
        return(plot)
}

fx_cforest <- function(model,modelName){
	print(paste(modelName,"_fx_cforest"))
	C1 <- fx_single_cond(model,1,meta)
	C2 <- fx_single_cond(model,1,iso)
	C3 <- fx_single_cond(model,32,meta)
	C4 <- fx_single_cond(model,32,iso)
        rf_df <- fx_cforest_df(C1,C2,C3,C4,modelName)
	plot <- fx_plot(rf_df,modelName)
        return(plot)
}
