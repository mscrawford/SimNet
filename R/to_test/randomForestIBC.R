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

### Grass3 (IBC-grass)
source(paste0(scripts_dir, "/to_test/readIBC.R"))

IBC_grass <- models$Grass3 
IBC_grass_traits <- IBC_grass_traits %>% mutate(SpeciesID = as.character(SpeciesID))

IBC_grass <- inner_join(IBC_grass, IBC_grass_traits) %>%
    select(Rep, Ninitial, SpeciesID, Year, Stage, Productivity, Biomass, LMR, MaxMass, Gmax, SLA, meanSpacerLength)

#IBC_grass <- models$Grass3 %>%
#    select(-Model, -SeedRain) %>%
#    mutate(SpeciesID = as.factor(SpeciesID))
#
#IBC_grass_traits <- IBC_grass_traits %>%
#    mutate(SpeciesID = as.factor(SpeciesID))
#
#IBC_grass <- inner_join(IBC_grass, IBC_grass_traits, by = c("SpeciesID"))

G3C1 <- randomForest_4conditions("G3C1",IBC_grass,1,meta)
G3C2 <- randomForest_4conditions("G3C2",IBC_grass,1,iso)
G3C3 <- randomForest_4conditions("G3C3",IBC_grass,32,meta)
G3C4 <- randomForest_4conditions("G3C4",IBC_grass,32,iso)
