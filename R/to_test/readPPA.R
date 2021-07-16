library(data.table)
library(tidyverse)
library(scales)

setwd(raw_data_dir)

# -------------------------------------------------------------------------
# PPA

PPA <- readRDS("PPA_Table1.rds") %>%
    filter(Year %in% seq(0, 2000, 10)) %>%
    mutate(Year = Year / 10) %>%
    select(-`F`, -N, -BasalArea)

PPA_traits <- fread("PPA_Table2.csv")

# PPA does not print species without biomass, so I include them by hand
PPA_initialCommunities <- fread("PPA_initialCommunities.csv")

PPA_initialCommunities <- PPA_initialCommunities %>%
    expand_grid(SeedRain = unique(PPA$SeedRain), # I only ran scenarios with these levels of seed rain
                Year = seq(1, 200)) %>%
    mutate(Model = "PPA",
           Stage = ifelse(Year > 100, "disassembly", "assembly"),
           Biomass = 0,
           Productivity = 0)

extirpated <- anti_join(PPA_initialCommunities, PPA,
                        by = c("Model", "Ninitial", "Rep", "SeedRain", "SpeciesID", "Stage", "Year"))

PPA <- bind_rows(PPA, extirpated) %>%
    mutate(Model = "Forest1")

rm(PPA_initialCommunities, extirpated)


# -------------------------------------------------------------------------
# General model formatting

model_runs <- list(PPA)

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

model_traits <- list(PPA_traits)

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
