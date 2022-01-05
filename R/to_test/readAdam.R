library(data.table)
library(tidyverse)
library(scales)

setwd(raw_data_dir)


# -------------------------------------------------------------------------
# Adam's model

load("adammod_trans_exp_out_HPC.rda")

adam <- datout %>%
    mutate(Biomass = Productivity,
           Model = "Grass1")

# Adam ran the 64 species simulation too many times
adam_64 <- adam %>%
    filter(Ninitial == 64) %>%
    group_by(Model, SeedRain, SpeciesID, Stage, Year) %>%
    filter(Rep == first(Rep))

adam_not64 <- adam %>%
    filter(Ninitial != 64)

adam <- bind_rows(adam_64, adam_not64) %>%
    ungroup()

# Remove all species with inappreciable biomass
adam <- adam %>%
    mutate(Biomass      = ifelse(Biomass < 1, 0, Biomass),
           Productivity = ifelse(Productivity < 1, 0, Productivity))

adam_traits <- read.csv("adammod_trans_exp_speciesdata.csv")

rm(datout, adam_64, adam_not64)

# General model formatting

model_runs <- list(adam
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

model_traits <- list(adam_traits
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
