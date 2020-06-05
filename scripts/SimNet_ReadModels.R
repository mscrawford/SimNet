library(data.table)
library(tidyverse)
library(scales)

setwd("../data/raw")

# Some analyses (those that include time to extinction in particular) want to keep each other
if (!exists("filterYears"))
{
    filterYears <- TRUE
}

# Read in model results -----------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# Adam

load("adammod_trans_exp_out_HPC.rda")

adam <- datout %>%
    mutate(Biomass = Productivity)

rm(datout)

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
    mutate(Biomass = ifelse(Biomass < 1,
                            0,
                            Biomass),
           Productivity = ifelse(Productivity < 1,
                            0,
                            Productivity))

adam_traits <- read.csv("adammod_trans_exp_speciesdata.csv")


# ---------------------------------------------------------------------------------------------
# Lindsay

load("lindsaymod_trans_exp_out_HPC.rda")

lindsay <- datout %>%
    mutate(Biomass = Productivity)

rm(datout)

# Adam ran the 64 species simulation too many times
lindsay_64 <- lindsay %>%
    filter(Ninitial == 64) %>%
    group_by(Model, SeedRain, SpeciesID, Stage, Year) %>%
    filter(Rep == first(Rep))

lindsay_not64 <- lindsay %>%
    filter(Ninitial != 64)

lindsay <- bind_rows(lindsay_64, lindsay_not64) %>%
    ungroup()

lindsay_traits <- read.csv("lindsaymod_trans_exp_speciesdata.csv")


# ---------------------------------------------------------------------------------------------
# IBC-grass

load("IBC-grass_Table1.rda")

IBC_grass <- d %>%
    tbl_df() %>%
    select(-SimID, -Stabilization) %>%
    mutate(Productivity = Biomass)

rm(d)

IBC_grass.noNDD <- IBC_grass %>% filter(Model == "IBC_grass.noNDD")
IBC_grass.NDD <- IBC_grass %>% filter(Model == "IBC_grass.NDD")

rm(IBC_grass)

IBC_grass_traits <- read.csv("IBC-grass_Table2.csv")

remove_cols <- IBC_grass_traits %>%
    map_dfr(var) %>%
    gather() %>%
    filter(value == 0) %>%
    spread(key, value) %>%
    names()

IBC_grass_traits <- IBC_grass_traits[, setdiff(names(IBC_grass_traits), remove_cols)]


# ---------------------------------------------------------------------------------------------
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
           Stage = ifelse(Year > 100,
                          "disassembly",
                          "assembly"),
           Biomass = 0,
           Productivity = 0)

extirpated <- anti_join(PPA_initialCommunities, PPA,
                        by = c("Model", "Ninitial", "Rep", "SeedRain", "SpeciesID", "Stage", "Year"))

PPA <- bind_rows(PPA, extirpated)


# ---------------------------------------------------------------------------------------------
# TROLL

troll <- readRDS("TROLL_Table1.rds") %>%
    select(-LineNumber) %>%
    mutate(SpeciesID = as.factor(SpeciesID),
           SpeciesID = as.numeric(SpeciesID))

troll_traits <- fread("TROLL_Table2.txt") %>%
    rename(SpeciesID = species_binomial) %>%
    mutate(SpeciesID = as.factor(SpeciesID),
           SpeciesID = as.numeric(SpeciesID))


# ---------------------------------------------------------------------------------------------
# SUCC

succ <- readRDS("bjoern_Table1_averagedSmooth.rds") %>%
    select(-Smooth) %>%
    ungroup()

succ <- succ %>%
    mutate(Productivity = replace(Productivity, Productivity < 0, 0))

succ_traits <- readRDS("bjoern_Table2.rds") %>%
    select(SpeciesID, maxSize, pLeaf, pRoot, pStorage)


# General model formatting --------------------------------------------------------------------

model_runs <- list(adam,
                   lindsay,
                   IBC_grass.NDD,
                   # IBC_grass.noNDD,
                   PPA,
                   troll,
                   succ)

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
                                                assembly = "metacommunity",
                                                disassembly = "isolation"))

                      .x <- .x %>%
                          filter(SeedRain %in% c(100))

                      .x$SeedRain <- droplevels(.x$SeedRain)

                      if (filterYears)
                      {
                          .x <- .x %>%
                              filter(Year %in% c(100, 200))
                      }

                      .x <- .x %>%
                          ungroup() %>% # There shouldn't be groups anyways
                          mutate(Biomass = scales::rescale(Biomass, to = c(0, 100))) %>%
                          mutate(Productivity = scales::rescale(Productivity, to = c(0, 100)))
                  })

model_traits <- list(adam_traits,
                     lindsay_traits,
                     IBC_grass_traits,
                     # IBC_grass_traits,
                     PPA_traits,
                     troll_traits,
                     succ_traits)

model_traits <- map(.x = model_traits,
                    .f = ~
                        {
                            .x <- .x %>%
                                mutate(SpeciesID = as.character(SpeciesID))
                        })


# Formulation of model datasets ---------------------------------------------------------------

models <- map2(.x = model_runs,
               .y = model_traits,
               .f = ~
                   {
                       inner_join(.x, .y, by = c("SpeciesID"))
                   })

names(models) = map(.x = models,
                    .f = ~
                        {
                            unique(.x$Model)
                        })

# Cleanup -------------------------------------------------------------------------------------

setwd("../../scripts")
rm(filterYears)
