library(data.table)
library(tidyverse)
library(ggthemes)
library(cowplot)
library(gganimate)
library(viridis)

setwd("/Users/Theodore/Documents/Dropbox/Work/Projects/SimNet/Modelling_Experiment/data")


# ---------------------------------------------------------------------------------------------
# Load model results
PPA <- fread("PPA_Table1.csv") %>%
    filter(Year %in% seq(0, 2000, 10)) %>%
    mutate(Year = Year / 10) %>%
    select(-`F`, -N, -BasalArea, -Productivity)

PPA_traits <- fread("PPA_Table2.csv")

PPA <- PPA %>%
    mutate(Model = as.character(Model),
           Ninitial = as.factor(Ninitial),
           Rep = as.numeric(Rep),
           SeedRain = as.factor(SeedRain),
           SpeciesID = as.character(SpeciesID),
           Stage = as.factor(Stage),
           Year = as.numeric(Year),
           Biomass = as.numeric(Biomass))

PPA <- PPA %>%
    mutate(Stage = recode(Stage,
                          assembly = "metacommunity",
                          disassembly = "isolation"))

PPA <- PPA %>%
    filter(SeedRain %in% c(100),
           Year %in% c(100, 200))

PPA <- PPA %>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))

PPA <- PPA %>%
    filter(Biomass > 0)

PPA_traits <- PPA_traits %>%
    mutate(SpeciesID = as.character(SpeciesID))

PPA <- inner_join(PPA, PPA_traits) %>% tbl_df()


# ---------------------------------------------------------------------------------------------
# Communities with 32 species

# For the 32-species communities, it seems like there are two clusters without seed addition

PPA.diversity <- PPA.diversity %>%
    filter(Ninitial == 32,
           Stage == "isolation")

clusters <- kmeans(PPA.diversity %>%
                       select(Shannon),
                   2)

PPA.diversity <- cbind(PPA.diversity,
                       as.data.frame(clusters$cluster)) %>%
    tbl_df() %>%
    rename(cluster = `clusters$cluster`)

p <- ggplot(PPA.diversity) +
    geom_histogram(aes(x = Shannon,
                       fill = as.factor(cluster))) +
    facet_grid(rows = vars(Ninitial),
               cols = vars(Stage),
               scales = "free_y") +
    theme_bw() +
    theme(aspect.ratio = 1)

save_plot(filename = "~/Desktop/kmeansClusters_32.pdf", plot = p, nrow = 1, ncol = 2)


# exemplary communities -----------------------------------------------------------------------

PPA.diversity <- PPA.diversity %>%
    group_by(cluster) %>%
    sample_n(20)

PPA.diversity <- PPA.diversity %>%
    select(Ninitial, SeedRain, Rep, Stage, cluster)

PPA.communities <- inner_join(PPA.diversity, PPA)

blah <- PPA.communities %>%
    group_by(Ninitial, SeedRain, Rep, Stage, cluster) %>%
    filter(Biomass > 0) %>%
    mutate(p_i = Biomass / sum(Biomass)) %>%
    summarise(total.biomass = sum(Biomass),
              Shannon = 1 - sum(p_i * log(p_i)),
              Richness = n())


a <- ggplot(PPA.communities %>% filter(cluster == 1),
            aes(x = PC1score,
                y = PC2score,
                size = Biomass,
                color = as.factor(SpeciesID))) +
    geom_point() +
    facet_wrap(facets = vars(Rep), nrow = 4) +
    theme_bw() +
    theme(aspect.ratio = 1)

b <- ggplot(PPA.communities %>% filter(cluster == 2),
            aes(x = PC1score,
                y = PC2score,
                size = Biomass,
                color = as.factor(SpeciesID))) +
    geom_point() +
    facet_wrap(facets = vars(Rep), nrow = 4) +
    theme_bw() +
    theme(aspect.ratio = 1)

p <- cowplot::plot_grid(a, b, nrow = 2, labels = unique(PPA.communities$cluster))

cowplot::save_plot(p, filename = "~/Desktop/Clusters.pdf", ncol = 5, nrow = 8)

# No clustering, communities with metacommunity and isolated
sample_reps <- sample(1:64, 10)

PPA.sampled <- PPA %>%
    filter(Ninitial == 32,
           Rep %in% sample_reps)

p <- ggplot(PPA.sampled) +
    geom_point(aes(x = PC1score,
                   y = PC2score,
                   size = Biomass,
                   color = SpeciesID)) +
    facet_grid(cols = vars(Stage), rows = vars(Rep)) +
    theme_bw() +
    theme(aspect.ratio = 1)

save_plot(p, filename = "~/Desktop/CommunitiesWithStage.pdf", ncol = 2, nrow = 10)



# ---------------------------------------------------------------------------------------------
# kmeans with 8 species mixtures --------------------------------------------------------------

# Load model results
PPA <- fread("PPA_Table1.csv") %>%
    filter(Year %in% seq(0, 2000, 10)) %>%
    mutate(Year = Year / 10) %>%
    select(-`F`, -N, -BasalArea, -Productivity)

PPA_traits <- fread("PPA_Table2.csv")

PPA <- PPA %>%
    mutate(Model = as.character(Model),
           Ninitial = as.factor(Ninitial),
           Rep = as.numeric(Rep),
           SeedRain = as.factor(SeedRain),
           SpeciesID = as.character(SpeciesID),
           Stage = as.factor(Stage),
           Year = as.numeric(Year),
           Biomass = as.numeric(Biomass))

PPA <- PPA %>%
    mutate(Stage = recode(Stage,
                          assembly = "metacommunity",
                          disassembly = "isolation"))

PPA <- PPA %>%
    filter(SeedRain %in% c(100),
           Year %in% c(100, 200))

PPA <- PPA %>%
    ungroup() %>% # There shouldn't be groups anyways
    mutate(Biomass = scales::rescale(Biomass,
                                     to = c(0, 100)))

PPA <- PPA %>%
    filter(Biomass > 0)

PPA_traits <- PPA_traits %>%
    mutate(SpeciesID = as.character(SpeciesID))

PPA <- inner_join(PPA, PPA_traits) %>%
    tbl_df()

PPA.diversity <- PPA %>%
    group_by(Ninitial, SeedRain, Rep, Stage) %>%
    filter(Biomass > 0) %>%
    mutate(p_i = Biomass / sum(Biomass)) %>%
    summarise(total.biomass = sum(Biomass),
              Shannon = 1 - sum(p_i * log(p_i)),
              Richness = n()) %>%
    ungroup()

PPA.diversity <- PPA.diversity %>%
    filter(Ninitial == 8,
           Stage == "isolation")

clusters <- kmeans(PPA.diversity %>%
                       select(Shannon),
                   2)

PPA.diversity <- cbind(PPA.diversity,
                       as.data.frame(clusters$cluster)) %>%
    tbl_df() %>%
    rename(cluster = `clusters$cluster`)

p <- ggplot(PPA.diversity) +
    geom_histogram(aes(x = Shannon,
                       fill = as.factor(cluster))) +
    facet_grid(rows = vars(Ninitial),
               cols = vars(Stage),
               scales = "free_y") +
    theme_bw() +
    theme(aspect.ratio = 1)

save_plot(filename = "~/Desktop/kmeansClusters_8.pdf", plot = p, nrow = 1, ncol = 2)

PPA.diversity <- PPA.diversity %>%
    group_by(cluster) %>%
    sample_n(8)

PPA.diversity <- PPA.diversity %>%
    select(Ninitial, SeedRain, Rep, Stage, cluster)

PPA.communities <- inner_join(PPA.diversity, PPA)

a <- ggplot(PPA.communities %>% filter(cluster == 1),
            aes(x = PC1score,
                y = PC2score,
                size = Biomass,
                color = as.factor(SpeciesID))) +
    geom_point() +
    lims(x = c(min(PPA.communities$PC1score),
               max(PPA.communities$PC1score)),
         y = c(min(PPA.communities$PC2score),
               max(PPA.communities$PC2score))) +
    scale_size(limits = c(0.001, 100)) +
    facet_wrap(facets = vars(Rep), ncol = 4) +
    theme_bw() +
    theme(aspect.ratio = 1)

b <- ggplot(PPA.communities %>% filter(cluster == 2),
            aes(x = PC1score,
                y = PC2score,
                size = Biomass,
                color = as.factor(SpeciesID))) +
    geom_point() +
    lims(x = c(min(PPA.communities$PC1score),
               max(PPA.communities$PC1score)),
         y = c(min(PPA.communities$PC2score),
               max(PPA.communities$PC2score))) +
    scale_size(limits = c(0.001, 100)) +
    facet_wrap(facets = vars(Rep), ncol = 4) +
    theme_bw() +
    theme(aspect.ratio = 1)

p <- cowplot::plot_grid(a, b, nrow = 2, labels = unique(PPA.communities$cluster))

cowplot::save_plot(p, filename = "~/Desktop/Clusters_8_species.pdf", ncol = 4, nrow = 4)

