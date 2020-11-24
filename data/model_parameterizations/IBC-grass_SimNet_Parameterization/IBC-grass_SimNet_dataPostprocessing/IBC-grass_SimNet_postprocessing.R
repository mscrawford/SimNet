library(data.table)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(gganimate)
library(viridis)

read_data <- function(data_dir, file_type)
{
    main_dir = getwd()
    setwd(data_dir)

    files <- list.files(full.names = T)
    files <- files[which(grepl(file_type, files))]

    d <- bind_rows(map(.x = files,
                       .f = read_csv, col_names = TRUE, na = "NA", progress = TRUE))

    setwd(main_dir)
    return(d %>% as_tibble())
}

combine_data <- function(df_list, key)
{
    purrr::reduce(df_list,
                  left_join, by = key)
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# Derive SimNet tables from raw IBC-grass output ----------------------------------------------

data.dir = "./raw_data/"
trait.dir = "."

parameters <- read_data(data.dir, "parameter") %>%
    select(SimID, ComNr,
           nPFTs,
           Stabilization,
           SeedInput)

populations <- read_data(data.dir, "population") %>%
    select(SimID,
           PFT,
           Year,
           Shootmass,
           Rootmass,
           Repro)

traits <- read_data(trait.dir, "SimNetCommunity") %>%
    rename(SpeciesID = ID) %>%
    select(-Species)

d <- combine_data(list(parameters, populations), c("SimID"))

d <- d %>%
    mutate(Model = ifelse(Stabilization == 0,
                          "IBC_grass.noNDD",
                          "IBC_grass.NDD"),
           Stage = ifelse(Year <= 100,
                          "assembly",
                          "disassembly")) %>%
    rename(Ninitial = nPFTs,
           Rep = ComNr,
           SeedRain = SeedInput,
           SpeciesID = PFT) %>%
    mutate(Biomass = Shootmass + Rootmass + Repro) %>%
    select(-Shootmass, -Rootmass, -Repro)

save(d, file = "./IBC-grass_Table1.rda")
write_csv(traits, path = "./IBC-grass_Table2.csv")


# Tests for correctness -----------------------------------------------------------------------

# Leibold plot
d.t <- d %>%
    filter(SeedRain %in% c(100),
           Year %in% c(100, 200))

d.t <- d.t %>%
    filter(Biomass > 0) %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage) %>%
    mutate(p_i = Biomass / sum(Biomass)) %>%
    summarise(Shannon = 1 - sum(p_i * log(p_i)),
              Richness = n(),
              Biomass = sum(Biomass))

ggplot(d.t) +
    geom_point(aes(x = Shannon,
                   y = Biomass,
                   color = as.factor(Ninitial))) +
    geom_smooth(aes(x = Shannon,
                    y = Biomass,
                    color = as.factor(Ninitial)),
                method = "lm") +
    geom_smooth(aes(x = Shannon,
                    y = Biomass),
                method = "lm") +
    facet_grid(Model ~ Stage) +
    theme_few()

# 64 species community
d.t <- d %>%
    filter(Ninitial == 64,
           SeedRain == 100)

ggplot(d.t) +
    geom_line(aes(x = Year,
                  y = Biomass,
                  color = as.factor(SpeciesID))) +
    geom_vline(xintercept = 100, linetype = 3)

# 64 species community Richness and Shannon
d.t <- d %>%
    filter(Ninitial == 64,
           SeedRain == 100)

d.t <- d.t %>%
    filter(Biomass > 0) %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
    mutate(p_i = Biomass / sum(Biomass)) %>%
    summarise(Shannon = 1 - sum(p_i * log(p_i)),
              Richness = n(),
              Biomass = sum(Biomass))

ggplot(d.t) +
    geom_line(aes(x = Year,
                  y = Richness))

ggplot(d.t) +
    geom_line(aes(x = Year,
                  y = Shannon))
