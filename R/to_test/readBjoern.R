library(data.table)
library(tidyverse)
library(scales)

setwd(raw_data_dir)

# -------------------------------------------------------------------------
# Succulent model

#succ <- readRDS("bjoern_Table1_averaged_smooth_NAreplaced0.rds") %>%
#    select(-Smooth) %>%
#    ungroup() %>%
#    mutate(Model = "Dryland")
#
#succ <- succ %>%
#    mutate(Productivity = replace(Productivity, Productivity < 0, 0))
#
#succ_traits <- readRDS("bjoern_Table2.rds") %>%
#    select(SpeciesID, maxSize, pLeaf, pRoot, pStorage)
#
#---
bjoern <- readRDS("bjoern_Table1_averaged_smooth_NAreplaced0.rds") %>%
    select(-Productivity, -Smooth)# %>%
#    mutate(Model = "Dryland")


bjoern_traits <- readRDS("bjoern_Table2.rds") %>%
    select(SpeciesID, maxSize, pLeaf, pRoot, pStorage)

bjoern <- bjoern %>%
    filter(!is.na(Biomass)) %>%
    filter(SeedRain %in% c(100))

bjoern <- bjoern %>%
    mutate(Stage = recode(Stage,
                          assembly = "metacommunity",
                          disassembly = "isolation"))

bjoern <- bjoern %>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))
# -------------------------------------------------------------------------
# General model formatting

model_runs <- list(bjoern)

model_runs <- map(.x = model_runs,
                  .f = ~ {
                      .x <- .x %>%
                          mutate(Model = as.character(Model),
                                 Ninitial = as.factor(Ninitial),
                                 Rep = as.numeric(Rep),
                                 SeedRain = as.factor(SeedRain),
                                 SpeciesID = as.character(SpeciesID),
                                 Stage = as.factor(Stage),
                                 Year = as.numeric(Year),
                                 Biomass = as.numeric(Biomass))#,
                                 #Productivity = as.numeric(Productivity))

                      .x <- .x %>%
                          mutate(Stage = recode(Stage,
                                                assembly = "With seed inflow",
                                                disassembly = "Without seed inflow"))

                      .x <- .x %>% filter(SeedRain %in% c(100))
                      .x$SeedRain <- droplevels(.x$SeedRain)

                      .x <- .x %>% filter(Year %in% c(100, 200))

                      as_tibble(.x)
                  })

model_traits <- list(bjoern_traits)

model_traits <- map(.x = model_traits,
                    .f = ~ .x %>% mutate(SpeciesID = as.character(SpeciesID)))


# -------------------------------------------------------------------------
# Formulation of model datasets

models <- map2(.x = model_runs,
               .y = model_traits,
               .f = ~ inner_join(.x, .y, by = c("SpeciesID")))

names(models) = map(.x = models,
                    .f = ~ unique(.x$Model))

setwd(base_dir)
