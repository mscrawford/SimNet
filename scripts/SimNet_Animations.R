library(data.table)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(gganimate)
library(viridis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SimNet_ReadModels.R")

# Shannon diversity vs. Biomass, entire datasets -----------------------------------------

plot_fun <- function(df, title)
{
    p <- ggplot(df,
                aes(x = Shannon,
                    y = total.biomass,
                    color = Stage,
                    group = Rep)) +
        geom_point() +
        facet_grid(Ninitial ~ Model,
        # facet_grid(Ninitial ~ SeedRain,
                   labeller = "label_both") +
        transition_time(Year) +
        shadow_wake(wake_length = 0.1) +
        labs(title = paste(title, ' | Year:  {frame_time}', sep = ""),
             x = "Shannon diversity",
             y = "Total Biomass") +
        theme_bw(20) +
        theme(aspect.ratio = 1)

    return(p)
}

d <- rbind(adam, lindsay, troll, PPA, IBC_grass.NDD, IBC_grass.noNDD, succ) %>%
    filter(Ninitial %in% c(1, 2, 4, 8, 16, 32, 64))

d <- d %>%
    filter(Biomass > 0) %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
    mutate(p_i = Biomass / sum(Biomass)) %>%
    summarise(total.biomass = sum(Biomass),
              Shannon = 1 - sum(p_i * log(p_i)),
              Richness = n())

p.allModels <- plot_fun(d, "All models")
anim_save(filename = "~/Desktop/allModels.gif", animation = p.allModels, width = 1400, height = 1400)

# p.adam <- plot_fun(d %>% filter(Model == "adam"), "Adam")
# anim_save(filename = "~/Desktop/adam.gif", animation = p.adam, width = 1400, height = 1400)
#
# p.lindsay <- plot_fun(d %>% filter(Model == "lindsay"), "Lindsay")
# anim_save(filename = "~/Desktop/lindsay.gif", animation = p.lindsay, width = 1400, height = 1400)
#
# p.troll <- plot_fun(d %>% filter(Model == "TROLL"), "TROLL")
# anim_save(filename = "~/Desktop/troll.gif", animation = p.troll, width = 1400, height = 1400)
#
# p.PPA <- plot_fun(d %>% filter(Model == "PPA"), "PPA")
# anim_save(filename = "~/Desktop/PPA.gif", animation = p.PPA, width = 1400, height = 1400)
#
# p.IBC_grass.NDD <- plot_fun(d %>% filter(Model == "IBC_grass.NDD"), "IBC_grass.NDD")
# anim_save(filename = "~/Desktop/IBC_grass.gif", animation = p.IBC_grass.NDD, width = 1400, height = 1400)
#
# p.IBC_grass.noNDD <- plot_fun(d %>% filter(Model == "IBC_grass.noNDD"), "IBC_grass.noNDD")
# anim_save(filename = "~/Desktop/IBC_grass.noNDD.gif", animation = p.IBC_grass.noNDD, width = 1400, height = 1400)
#
# p.succ <- plot_fun(d %>% filter(Model == "succ"), "succ")
# anim_save(filename = "~/Desktop/succ.gif", animation = p.succ, width = 1400, height = 1400)
