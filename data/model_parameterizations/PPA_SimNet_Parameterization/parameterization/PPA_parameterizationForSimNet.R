library(tidyverse)

set.seed(1)

# Parameterization for SimNet of the PPA model

# This script reads in the SimNet species' vita rates and runs a kmeans clustering on them,
# to sample 64 species from the mixture while maintaining maximal functional diversity.


# Format species dataset ----------------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Format species dataset ----------------------------------------------------------------------

PPA_vitalRates = read.csv("../info/SimNet_PPA_vitalRates.txt", sep = "\t", header = TRUE)

PPA_traits <- read.csv("../info/Data_S1.txt", sep = "\t") %>%
    select(sp, PC1score, PC2score)

data <- inner_join(PPA_vitalRates, PPA_traits)

data <- data %>%
    rename(SpeciesID = sp) %>%
    mutate(SpeciesID = as.factor(SpeciesID),
           SpeciesID = as.numeric(SpeciesID)) %>%
    tbl_df()


# Select 64 species through k-means clustering -----------------------------------------------

clusters <- kmeans(data %>%
                       select(G1, G2, G3, G4,
                              mu1, mu2, mu3, mu4,
                              `F`,
                              wd),
                   64)

data <- cbind(data, clusters$cluster)

data <- data %>%
    group_by(`clusters$cluster`) %>%
    sample_n(1)

data <- data %>%
    ungroup() %>%
    select(-`clusters$cluster`)

data <- data %>%
    arrange(SpeciesID)


# Write PPA 64 species community to disk ------------------------------------------------------

vital_rates <- data %>%
    select(-PC1score, -PC2score)

write_csv(x = vital_rates, path = "../input/PPA_vitalRates.csv")

traits <- data %>%
    select(SpeciesID, PC1score, PC2score)

write_csv(x = traits, path = "../input/PPA_traits.csv")
