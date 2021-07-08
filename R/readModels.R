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


# -------------------------------------------------------------------------
# Lindsay's model

load("lindsaymod_trans_exp_out_HPC.rda")

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

lindsay_traits <- read.csv("lindsaymod_trans_exp_speciesdata.csv")

rm(datout, lindsay_64, lindsay_not64)


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

model_runs <- list(adam,
                   lindsay,
                   IBC_grass.NDD,
                   # IBC_grass.noNDD,
                   PPA,
                   troll,
                   bjoern)

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

model_traits <- list(adam_traits,
                     lindsay_traits,
                     IBC_grass_traits,
                     PPA_traits,
                     troll_traits,
                     bjoern_traits)

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
