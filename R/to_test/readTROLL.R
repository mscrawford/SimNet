library(data.table)
library(tidyverse)
library(scales)

setwd(raw_data_dir)

# -------------------------------------------------------------------------
# TROLL

troll <- readRDS("TROLL_Table1.rds") %>%
    select(-LineNumber) %>%
    mutate(SpeciesID = as.factor(SpeciesID),
           SpeciesID = as.numeric(SpeciesID)) %>%
    mutate(Model = "Forest2")

troll_traits <- fread("TROLL_Table2.txt") %>%
    rename(SpeciesID = species_binomial) %>%
    mutate(SpeciesID = as.factor(SpeciesID),
           SpeciesID = as.numeric(SpeciesID))

# General model formatting

model_runs <- list(troll
                   )

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
                                 Biomass = as.numeric(Biomass),
                                 Productivity = as.numeric(Productivity))

                      .x <- .x %>%
                          mutate(Stage = recode(Stage,
                                                assembly = "With seed inflow",
                                                disassembly = "Without seed inflow"))

                      .x <- .x %>% filter(SeedRain %in% c(100))
                      .x$SeedRain <- droplevels(.x$SeedRain)

                      .x <- .x %>% filter(Year %in% c(100, 200))

                      as_tibble(.x)
                  })

model_traits <- list(troll_traits
                     )

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
