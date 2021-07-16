library(randomForest)

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path)); setwd("../../")
#base_dir          <- getwd()
base_dir          <- setwd("../../")
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")

#source(paste0(scripts_dir, "/readModels.R"))

meta = seq(80, 100)
iso = seq(180, 200)

randomForest_4conditions <- function(modelName,model,NoSpp,stage){
	print(modelName)
        model <- model %>%
            filter(Ninitial == NoSpp,
                   Year %in% stage) %>%
            mutate(id = row_number()) %>%
            mutate_if(is.character, as.factor) %>%
            select(-Productivity,-SpeciesID, -Ninitial, -Stage, -Rep) 
        
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
        varImpPlot(rf, main = title)
        #dev.set(dev.next())
        while (!is.null(dev.list()))  dev.off()
        return(rf)
}

### Dryland (Bjoern)
source(paste0(scripts_dir, "/to_test/readBjoern.R"))

bjoern <- models$bjoern %>%#Dryland %>% 
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))

bjoern <- bjoern %>%
    select(-Model, -SeedRain) %>%
    mutate(SpeciesID = as.factor(SpeciesID))

bjoern_traits <- model_traits[[6]] %>%
    mutate(SpeciesID = as.factor(SpeciesID))

bjoern <- inner_join(bjoern, bjoern_traits, by = c("SpeciesID"))

DC1 <- randomForest_4conditions("DC1",bjoern,1,meta)
DC2 <- randomForest_4conditions("DC2",bjoern,1,iso)
DC3 <- randomForest_4conditions("DC3",bjoern,32,meta)
DC4 <- randomForest_4conditions("DC4",bjoern,32,iso)

