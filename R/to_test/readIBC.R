library(data.table)
library(tidyverse)
library(scales)

setwd(raw_data_dir)


# ---------------------------------------------------------------------------------------------
# IBC-grass

load("IBC-grass_Table1.rda")

IBC_grass <- d %>%
    as_tibble() %>%
    select(-SimID, -Stabilization) %>%
    mutate(Productivity = Biomass)

IBC_grass.noNDD <- IBC_grass %>% filter(Model == "IBC_grass.noNDD")
IBC_grass.NDD <- IBC_grass %>% filter(Model == "IBC_grass.NDD") %>%
    mutate(Model = "Grass3")

IBC_grass_traits <- read.csv("IBC-grass_Table2.csv")

remove_cols <- IBC_grass_traits %>%
    map_dfr(var) %>%
    gather() %>%
    filter(value == 0) %>%
    spread(key, value) %>%
    names()

IBC_grass_traits <- IBC_grass_traits[, setdiff(names(IBC_grass_traits), remove_cols)]

rm(d, IBC_grass, remove_cols)

# -------------------------------------------------------------------------
# General model formatting

model_runs <- list(IBC_grass.NDD
		   ,IBC_grass.noNDD)

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

model_traits <- list(IBC_grass_traits)

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
