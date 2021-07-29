library(party)
library(ggplot2)
library(dplyr)
set.seed(1987)

meta = seq(80, 100)
iso = seq(180, 200)

fx_cforest_party <- function(modelName,model,NoSpp,stage){
	print(modelName)
        model <- model %>%
            filter(Ninitial == NoSpp,
                   Year %in% stage) %>%
            mutate(id = row_number()) %>%
            mutate_if(is.character, as.factor) %>%
            select(-Productivity,-SpeciesID, -Ninitial, -Stage, -Rep, -Year) 
        
        train <- model %>% sample_frac(.50)
        test <- anti_join(model, train, by = 'id')
        
        train <- train %>% select(-id)
        test <- test %>% select(-id)
        
        rf <- cforest(Biomass ~ .
                      ,data = as.data.frame(train)
		      ,control = cforest_unbiased(mtry = 2, ntree = 501)
		      )
        pred <- data.frame(pred = predict(rf, test,OOB=TRUE))
        title = paste("Correlation: ", round(cor(pred, test$Biomass)[[1]], 2),
        #               "  |  mean squared error: ", round(mean(rf$mse), 2),
        #               "  |  R-squared: ", round(mean(rf$rsq), 2),
        sep = "")
        set.seed(1987)
	#Permutation importance:
        PI <- varimp(rf)
        #Conditional permutation importance:
        CPI <- varimp(rf, conditional = TRUE)
        #Prepare data frame for plotting
        df1 <- as.data.frame(PI)
        df1$varnames <- rownames(df1)
        rownames(df1) <- NULL
        df1$var_categ <- c(1: dim(df1)[1])
        #Plot
        ggplot(df1, aes(x=reorder(varnames, PI)
                       ,y=PI
                       ,color=as.factor(var_categ))) +
                geom_point(show.legend = FALSE) +
          geom_text(aes(label=sprintf("%.2f", PI))
                    ,hjust=0,vjust=0,color="black",size=7
                    ,angle = -90,show.legend = FALSE) +
          geom_segment(aes(x=varnames,xend=varnames,y=0,yend=PI)
                       ,show.legend = FALSE) +
          scale_color_discrete(name="Variable Group") +
          ylab("Origional permutation importance") +
          xlab("Trait") +
          theme_bw() +
          theme(panel.border = element_blank()
                ,axis.line = element_line(colour = "black")
                ,text = element_text(size = 22)
                ,axis.text.x = element_text(size = 20)
                ,axis.text.y = element_text(size = 20)) +
                coord_flip() +
                labs(title = title
                     #,subtitle = ""
                     #,caption = ""
                     )
        ggsave(paste0(tmp_dir,"/randomForest/",modelName,"_orig.pdf"))
        while (!is.null(dev.list()))  dev.off()
        #Prepare data frame for plotting
	df <- as.data.frame(CPI)
	df$varnames <- rownames(df)
	rownames(df) <- NULL
	df$var_categ <- c(1: dim(df)[1])
	#Plot
	ggplot(df, aes(x=reorder(varnames, CPI)
	               ,y=CPI
	               ,color=as.factor(var_categ))) +
	        geom_point(show.legend = FALSE) +
	        geom_text(aes(label=sprintf("%.2f", CPI))
	                  ,hjust=0,vjust=0,color="black",size=7
	                  ,angle = -90,show.legend = FALSE) +
	        geom_segment(aes(x=varnames,xend=varnames,y=0,yend=CPI)
	                     ,show.legend = FALSE) +
	        scale_color_discrete(name="Variable Group") +
	        ylab("Conditional permutation importance") +
	        xlab("Trait") +
	        theme_bw() +
	        theme(panel.border = element_blank()
	              ,axis.line = element_line(colour = "black")
	              ,text = element_text(size = 22)
	              ,axis.text.x = element_text(size = 20)
	              ,axis.text.y = element_text(size = 20)) +
	        coord_flip() +
	        labs(title = title)
	ggsave(paste0(tmp_dir,"/randomForest/",modelName,".pdf"))
        while (!is.null(dev.list()))  dev.off()
        return(rf)
}
