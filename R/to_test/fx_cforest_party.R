library(party)
set.seed(42)

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
        
        train <- model %>% sample_frac(.70)
        test <- anti_join(model, train, by = 'id')
        
        train <- train %>% select(-id)
        test <- test %>% select(-id)
        
        rf <- cforest(Biomass ~ .
                      ,data = as.data.frame(train)
		      ,control = cforest_unbiased(mtry = 2, ntree = 501)
		      )
        #predict(rf,OOB = TRUE)
#        pred <- data.frame(pred = predict(rf, test),OOB=TRUE)
 #       title = paste("correlation: ", round(cor(pred, test$Biomass)[[1]], 2),
 #                     "  |  mean squared error: ", round(mean(rf$mse), 2),
 #                     "  |  R-squared: ", round(mean(rf$rsq), 2),
 #                     sep = "")
        pdf(paste0(tmp_dir,"/randomForest/",modelName,".pdf"))
#        varImpPlot(rf,main = title)#, cex=1.2)
	varimp(rf, conditional = TRUE)
        #dev.set(dev.next())
        while (!is.null(dev.list()))  dev.off()
	print(varimp(rf))
	print(varimp(rf, conditional = TRUE))
        return(rf)
}
