#library(data.table)
#library(tidyverse)
#library(ggthemes)
#library(gganimate)
#library(viridis) 
#library(cowplot)
#library(scales)
#library(plotly)
library(randomForest)
#library(pdp)

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

### Grass1 (Adam's model)
source(paste0(scripts_dir, "/to_test/readAdam.R"))

adam <- models$Grass1

adam <- adam %>%
  ungroup() %>% # There shouldn't be groups anyways
  mutate(Biomass = scales::rescale(Biomass,
                                   to = c(0, 100)))

adam <- adam %>%
    select(-Model, -SeedRain) %>%
    mutate(SpeciesID = as.factor(SpeciesID))

adam_traits <- adam_traits %>%
    mutate(SpeciesID = as.factor(SpeciesID))

adam <- inner_join(adam, adam_traits, by = c("SpeciesID"))

G1C1 <- randomForest_4conditions("G1C1",adam,1,meta)
#G1C2 <- randomForest_4conditions("G1C2",adam,1,iso)
#G1C3 <- randomForest_4conditions("G1C3",adam,32,meta)
#G1C4 <- randomForest_4conditions("G1C4",adam,32,iso)

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
#G2C2 <- randomForest_4conditions("G2C2",lindsay,1,iso)
#G2C3 <- randomForest_4conditions("G2C3",lindsay,32,meta)
#G2C4 <- randomForest_4conditions("G2C4",lindsay,32,iso)

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
#G3C2 <- randomForest_4conditions("G3C2",IBC_grass,1,iso)
#G3C3 <- randomForest_4conditions("G3C3",IBC_grass,32,meta)
#G3C4 <- randomForest_4conditions("G3C4",IBC_grass,32,iso)

### Forest1 (PPA)
source(paste0(scripts_dir, "/to_test/readPPA.R"))

PPA <- models$Forest1 %>%
    select(-Model, -SeedRain)

F1C1 <- randomForest_4conditions("F1C1",PPA,1,meta)
#F1C2 <- randomForest_4conditions("F1C2",PPA,1,iso)
#F1C3 <- randomForest_4conditions("F1C3",PPA,32,meta)
#F1C4 <- randomForest_4conditions("F1C4",PPA,32,iso)

### Forest2 (TROLL)
source(paste0(scripts_dir, "/to_test/readTROLL.R"))
troll <- models$Forest2 

troll <- troll %>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))

troll <- troll %>%
    select(-Model, -SeedRain) %>%
    mutate(SpeciesID = as.factor(SpeciesID))

troll_traits <- troll_traits %>%
    mutate(SpeciesID = as.factor(SpeciesID))

troll <- inner_join(troll, troll_traits, by = c("SpeciesID"))

troll <- troll %>%
    mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
    select(-hmax, -ah)

F2C1 <- randomForest_4conditions("F2C1",troll,1,meta)
#F2C2 <- randomForest_4conditions("F2C2",troll,1,iso)
#F2C3 <- randomForest_4conditions("F2C3",troll,32,meta)
#F2C4 <- randomForest_4conditions("F2C4",troll,32,iso)

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
#DC2 <- randomForest_4conditions("DC2",bjoern,1,iso)
#DC3 <- randomForest_4conditions("DC3",bjoern,32,meta)
#DC4 <- randomForest_4conditions("DC4",bjoern,32,iso)

