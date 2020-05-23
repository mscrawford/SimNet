library(data.table)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(gganimate)
library(viridis)

source("~/Documents/Dropbox/Work/Projects/SimNet/Modelling_Experiment/scripts/SimNet_ReadModels.R")

# ---------------------------------------------------------------------------------------------

# Biomass over time for each species in 64 species communities, with SeedRain gradient --------

d <- bind_rows(adam, lindsay, troll, PPA) %>%
    filter(Ninitial %in% c(64))

p <- ggplot(d) +
    geom_line(aes(x = Year,
                  y = Biomass,
                  color = as.factor(SpeciesID),
                  group = as.factor(SpeciesID))) +
    geom_hline(yintercept = 0,
               linetype = 3) +
    geom_vline(xintercept = 100,
               linetype = 3) +
    facet_grid(Model ~ SeedRain,
               scales = "free",
               labeller = label_both) +
    labs(x = "Year",
         y = "Biomass") +
    theme_few() +
    theme(legend.background = element_blank(),
          legend.position = "none")

save_plot(p,
          filename = "64_speciesCommunities_overTime.pdf",
          ncol = 2,
          nrow = 4)


# Community biomass or Shannon over time ------------------------------------------------------

d <- bind_rows(adam, lindsay, troll, PPA) %>%
    filter(SeedRain %in% c(100))

d <- d %>%
    filter(Biomass > 0) %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
    mutate(p_i = Biomass / sum(Biomass)) %>%
    summarise(total.biomass = sum(Biomass),
              Shannon = 1 - sum(p_i * log(p_i)),
              Richness = n())

p <- ggplot(d) +
    geom_line(aes(x = Year,
                  y = Shannon,
                  group = Rep),
              alpha = 0.15) +
    geom_vline(xintercept = 100,
               linetype = 3) +
    facet_grid(Model ~ Ninitial, scales = "free") +
    theme_few()

save_plot(p,
          filename = "Shannon_perSeedRain_overTime.pdf",
          ncol = 7,
          nrow = 5,
          base_aspect_ratio = 1)


# TROLL 64 species equilibrium  ---------------------------------------------------------------

d <- troll %>%
    filter(Ninitial == 64,
           SeedRain %in% c(0, 100))

p <- ggplot(d,
            aes(x = Year,
                y = Biomass,
                color = as.factor(SpeciesID))) +
    geom_line() +
    geom_vline(xintercept = 100,
               linetype = 3) +
    facet_grid(~ SeedRain) +
    theme_few() +
    theme(legend.background = element_blank(),
          legend.position = "none"); p

save_plot(p, filename = "~/Desktop/TROLL_64_species_BiomassOverTime.pdf",
          ncol = 2, nrow = 1)



# TROLL, relationship between stand productivity and biomass for monocultures -----------------

d <- troll %>%
    filter(Ninitial == 1,
           SeedRain %in% c(0, 100),
           Year %in% c(100))

p <- ggplot(d,
            aes(x = Biomass,
                y = Productivity,
                color = as.factor(SeedRain))) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_few() +
    theme(aspect.ratio = 1); p

save_plot(p, filename = "~/Desktop/TROLL_64_Productivity_vs_Monoculture.pdf",
          ncol = 1, nrow = 1)


# Adam and Lindsay rank-abundance distributions -----------------------------------------------

d <- adam %>%
    filter(Ninitial %in% c(64),
           SeedRain %in% c(100),
           Year %in% c(200))

d <- d %>%
    group_by(Ninitial, SeedRain, Rep, Year) %>%
    arrange(desc(Biomass)) %>%
    mutate(rank = row_number())

a <- ggplot(d,
            aes(x = rank,
                y = log10(Biomass),
                group = Rep)) +
    geom_line(alpha = 1) +
    geom_point(alpha = 1) +
    geom_hline(yintercept = 0, linetype = 3) +
    facet_grid(. ~ Stage) +
    theme_few() +
    theme(aspect.ratio = 1)

b <- ggplot(d %>% filter(Biomass > 1),
            aes(x = rank,
                y = Biomass,
                group = Rep)) +
    geom_line(alpha = 1) +
    geom_point(alpha = 1) +
    geom_hline(yintercept = 0, linetype = 3) +
    facet_grid(. ~ Stage) +
    theme_few() +
    theme(aspect.ratio = 1)

p <- plot_grid(a, b,
               nrow = 1, ncol = 2,
               labels = c("Unfiltered", "> 1 g")); p

save_plot(p, filename = "~/Desktop/Adam_rankAbdunance.pdf",
          nrow = 1, ncol = 2,
          base_asp = 1)


# PPA row of equal biomasses ------------------------------------------------------------------

d <- PPA %>%
    filter(SeedRain %in% c(100),
           Year %in% c(100, 200))

d <- d %>%
    filter(Biomass > 0) %>%
    group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
    mutate(p_i = Biomass / sum(Biomass)) %>%
    summarise(total.biomass = sum(Biomass),
              Shannon = 1 - sum(p_i * log(p_i)),
              Richness = n())

ggplot(d) +
    geom_point(aes(x = Shannon,
                   y = total.biomass)) +
    facet_grid(. ~ Stage)

d <- d %>%
    filter(Stage == "disassembly")

d <- d %>%
    filter(Ninitial > 1,
           total.biomass > 550,
           total.biomass < 700,
           Shannon < 1.5)

shannons <- d

p <- ggplot(d) +
    geom_point(aes(x = Shannon,
                   y = total.biomass,
                   color = as.factor(Ninitial))) +
    facet_grid(. ~ Stage) +
    ylim(c(0, 1000))

save_plot(p,
          filename = "~/Desktop/PRI2CO_Communities.pdf")

d.filter <- d %>%
    ungroup() %>%
    select(Model, Ninitial, Rep, SeedRain)

d <- semi_join(PPA, d.filter)

p <- ggplot(d %>% filter(Ninitial == 32)) +
    geom_line(aes(x = Year * 10,
                  y = Biomass,
                  color = as.factor(SpeciesID),
                  group = paste(as.factor(SpeciesID), "_",
                                as.factor(Ninitial), "_",
                                as.factor(Rep)))) +
    geom_vline(xintercept = 1000, linetype = 3) +
    facet_wrap(. ~ Rep, labeller = label_both) +
    labs(x = "Year",
         y = "Biomass") +
    theme_few() +
    theme(legend.background = element_blank(),
          legend.position = "none"); p

save_plot(p,
          filename = "~/Desktop/PPA_TotalBiomassLine.pdf",
          ncol = 5,
          nrow = 4)

d <- d %>% filter(Year == 200,
                  Biomass > 300)

d <- inner_join(d, PPA_traits) %>%
    distinct(SpeciesID, PC1score, PC2score)
