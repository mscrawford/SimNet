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

### Grass2 (Lindsay's model)
source(paste0(scripts_dir, "/to_test/readLindsay.R"))

lindsay <- models$Grass2

lindsay <- lindsay %>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))
lindsay <- lindsay %>%
    select(-Model, -SeedRain) %>%
    mutate(SpeciesID = as.factor(SpeciesID))

lindsay_traits <- lindsay_traits %>%
    mutate(SpeciesID = as.factor(SpeciesID))

lindsay <- inner_join(lindsay, lindsay_traits, by = c("SpeciesID"))

G2C1 <- randomForest_4conditions("G2C1",lindsay,1,meta)
G2C2 <- randomForest_4conditions("G2C2",lindsay,1,iso)
#G2C3 <- randomForest_4conditions("G2C3",lindsay,32,meta)
#G2C4 <- randomForest_4conditions("G2C4",lindsay,32,iso)
