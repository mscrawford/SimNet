library(tidyverse)
library(ggthemes)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("./SimNet_ReadModels.R")
setwd("../tmp")

models <- model_runs %>%
    bind_rows() %>%
    filter(Ninitial != 64)

d <- models %>%
    group_by(Model, Ninitial, SpeciesID, Stage) %>%
    summarise(Biomass = mean(Biomass))

d.mono <- d %>%
    filter(Ninitial == 1) %>%
    rename(monoculture.biomass = Biomass) %>%
    ungroup() %>%
    select(-Ninitial)

d <- inner_join(d.mono, d) %>%
    filter(Ninitial != 1,
           Model != "IBC_grass.noNDD")

d$Ninitial <- droplevels(d$Ninitial)

p <- ggplot(d %>% filter(Ninitial == 32)) +
    geom_point(aes(x = monoculture.biomass,
                   y = Biomass,
                   color = SpeciesID),
               size = 3) +
    facet_grid(cols = vars(Stage),
               rows = vars(Model),
               scales = "free") +
    labs(x = "Monoculture biomass",
         y = "Biomass in 32-species mixture") +
    guides(color = FALSE) +
    theme_bw(18) +
    theme(aspect.ratio = 1); p

cowplot::save_plot(p, filename = "Competitiveness.png", ncol = 6, nrow = 2, base_asp = 1)
