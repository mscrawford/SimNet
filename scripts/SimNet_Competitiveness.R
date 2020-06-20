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

d <- d %>%
    filter(Ninitial == 32) %>%
    group_by(Model) %>%
    nest() %>%
    mutate(plot = map(.x = data,
                      .f = ~ {
                          ggplot(.x) +
                              geom_point(aes(x = monoculture.biomass,
                                             y = Biomass,
                                             color = SpeciesID),
                                         size = 3,
                                         show.legend = FALSE) +
                              facet_grid(cols = vars(Stage)) +
                              scale_color_viridis_d() +
                              labs(x = "Monoculture biomass",
                                   y = "Biomass in mixture") +
                              theme_bw(16) +
                              theme(aspect.ratio = 0.618)
                      })
           )

p <- cowplot::plot_grid(plotlist = d$plot,
                        labels = c("Adam", "Björn", "IBC-grass", "Lindsay", "PPA", "TROLL"),
                        ncol = 2,
                        nrow = 3)

cowplot::save_plot(p, filename = "Competitiveness.png", ncol = 4, nrow = 3, base_asp = 1.5)
