library(data.table)
library(tidyverse)
library(scales)

setwd(raw_data_dir)

# -------------------------------------------------------------------------
# Lindsay's model

load("lindsaymod_trans_exp_out_HPC_variable.rda")

lindsay <- datout %>%
    mutate(Biomass = Productivity,
           Model = "Grass2")

# Adam ran the 64 species simulation too many times
lindsay_64 <- lindsay %>%
    filter(Ninitial == 64) %>%
    group_by(Model, SeedRain, SpeciesID, Stage, Year) %>%
    filter(Rep == first(Rep))

lindsay_not64 <- lindsay %>%
    filter(Ninitial != 64)

lindsay <- bind_rows(lindsay_64, lindsay_not64) %>%
    ungroup()

lindsay_traits <- read.csv("lindsaymod_trans_exp_speciesdata_variable.csv")

rm(datout, lindsay_64, lindsay_not64)

# General model formatting

model_runs <- list(lindsay
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

model_traits <- list(lindsay_traits
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
