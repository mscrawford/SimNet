library(data.table)
library(tidyverse)
library(ggthemes)
library(gganimate)
library(viridis) 
library(cowplot)
library(scales)
library(plotly)
library(randomForest)
library(pdp)

## Read models
setwd("../")
base_dir          <- getwd()
scripts_dir       <- paste0(base_dir, "/R")
tmp_dir           <- paste0(base_dir, "/tmp")
raw_data_dir      <- paste0(base_dir, "/data/raw")
meta = seq(80, 100)
iso = seq(180, 200)

source(paste0(scripts_dir, "/readTROLL.R"))

troll <- troll %>%
    select(-Productivity,-Model, -SeedRain) %>%
    mutate(SpeciesID = as.factor(SpeciesID))

troll_traits <- troll_traits %>%
    mutate(SpeciesID = as.factor(SpeciesID))

troll <- inner_join(troll, troll_traits, by = c("SpeciesID"))

troll <- troll %>%
    mutate(h_realmax = hmax * dmax / (dmax + ah)) %>%
    select(-hmax, -ah)

troll.32 <- troll %>%
    filter(Ninitial == 32,
           Year %in% meta) %>%
    mutate(id = row_number()) %>%
    mutate_if(is.character, as.factor) %>%
    select(-SpeciesID, -Ninitial, -Stage, -Rep) 

train <- troll.32 %>% sample_frac(.70)
test <- anti_join(troll.32, train, by = 'id')

train <- train %>% select(-id)
test <- test %>% select(-id)

rf.32 <- randomForest(data = as.data.frame(train),
                                Biomass ~ .,
                                importance = TRUE)

pred <- data.frame(pred = predict(rf.32, test))
title = paste("correlation: ", round(cor(pred, test$Biomass)[[1]], 2),
              "  |  mean squared error: ", round(mean(rf.32$mse), 2),
              "  |  R-squared: ", round(mean(rf.32$rsq), 2),
              sep = "")
pdf(paste0(tmp_dir,"/randomForest/F2C3.pdf"))
varImpPlot(rf.32, main = title)
while (!is.null(dev.list()))  dev.off()

troll.32 <- troll %>%
    filter(Ninitial == 32,
           Year %in% iso) %>%
    mutate(id = row_number()) %>%
    mutate_if(is.character, as.factor) %>%
    select(-SpeciesID, -Ninitial, -Stage, -Rep) 

train <- troll.32 %>% sample_frac(.70)
test <- anti_join(troll.32, train, by = 'id')

train <- train %>% select(-id)
test <- test %>% select(-id)

rf.32 <- randomForest(data = as.data.frame(train),
                                Biomass ~ .,
                                importance = TRUE)

pred <- data.frame(pred = predict(rf.32, test))
title = paste("correlation: ", round(cor(pred, test$Biomass)[[1]], 2),
              "  |  mean squared error: ", round(mean(rf.32$mse), 2),
              "  |  R-squared: ", round(mean(rf.32$rsq), 2),
              sep = "")
pdf(paste0(tmp_dir,"/randomForest/F2C4.pdf"))
varImpPlot(rf.32, main = title)
while (!is.null(dev.list()))  dev.off()
