library(randomForest)

meta = seq(80, 100)
iso = seq(180, 200)

randomForest_4conditions <- function(modelName,model,NoSpp,stage){
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
        
        rf <- randomForest(data = as.data.frame(train),
                                        Biomass ~ .,
                                        importance = TRUE)
        
        pred <- data.frame(pred = predict(rf, test))
        title = paste("correlation: ", round(cor(pred, test$Biomass)[[1]], 2),
                      "  |  mean squared error: ", round(mean(rf$mse), 2),
                      "  |  R-squared: ", round(mean(rf$rsq), 2),
                      sep = "")
        pdf(paste0(tmp_dir,"/randomForest/",modelName,".pdf"))
        varImpPlot(rf,main = title)#, cex=1.2)
        #dev.set(dev.next())
        while (!is.null(dev.list()))  dev.off()
        print(rf$importance)
        return(rf)
}
