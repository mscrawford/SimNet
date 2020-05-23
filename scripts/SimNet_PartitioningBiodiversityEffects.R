library(data.table)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(gganimate)
library(viridis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SimNet_ReadModels.R")
setwd("../tmp")

# Additive partitioning of biodiversity effects -----------------------------------------------

d <- map(.x = models,
         .f = ~
             {
                 .x <- .x %>%
                     select(Model, Ninitial, Rep, SeedRain, Stage, Year, SpeciesID, Biomass)

                 .x <- .x %>%
                     mutate(Ninitial = as.numeric(Ninitial))
             }) %>%
    bind_rows() %>%
    tbl_df() %>%
    filter(Model %in% c("adam", "lindsay"))

monocultures <- d %>%
    filter(Ninitial == 1) %>%
    select(-Ninitial, -Rep) %>%
    rename(M_i = Biomass)

d <- inner_join(d, monocultures) %>%
    filter(Ninitial > 1) %>%
    rename(Y_Oi = Biomass)

d <- d %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
    mutate(Y_O = sum(Y_Oi),
           RY_Ei = 1 / Ninitial)

d <- d %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage, Year, SpeciesID) %>%
    mutate(RY_Oi = Y_Oi / M_i,
           Y_Ei = RY_Ei * M_i)

d <- d %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
    mutate(Y_E = sum(Y_Ei))

d <- d %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage, Year, SpeciesID) %>%
    mutate(delta.Y = Y_O - Y_E,
           delta.RY_i = RY_Oi - RY_Ei)

d <- d %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
    summarise(complementarity.effect = n() * mean(delta.RY_i) * mean(M_i),
              selection.effect = n() * cov(delta.RY_i, M_i),
              net.biodiversity.effect = complementarity.effect + selection.effect)

p <- ggplot(d %>% filter(Year > 10)) +
    geom_line(aes(x = Year,
                  y = selection.effect,
                  group = Rep),
              alpha = 0.2) +
    geom_vline(xintercept = 100,
               linetype = 3) +
    facet_grid(Model ~ Ninitial,
               scales = "free",
               labeller = label_both) +
    labs(x = "Year",
         y = "Selection effect") +
    theme_few(); p

save_plot(filename = "SelectionEffect.pdf", plot = p, ncol = 6, nrow = 4)

p <- ggplot(d %>% filter(Year > 10)) +
    geom_line(aes(x = Year,
                  y = complementarity.effect,
                  group = Rep),
              alpha = 0.2) +
    geom_vline(xintercept = 100,
               linetype = 3) +
    facet_grid(Model ~ Ninitial,
               scales = "free",
               labeller = label_both) +
    labs(x = "Year",
         y = "Complementarity effect") +
    theme_few(); p

save_plot(filename = "ComplementarityEffect.pdf", plot = p, ncol = 6, nrow = 4)
