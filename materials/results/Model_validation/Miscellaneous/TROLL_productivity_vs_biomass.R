library(data.table)
library(tidyverse)
library(scales)

setwd("/Users/Theodore/Documents/Dropbox/Work/Projects/SimNet/scripts")
setwd("../data/raw")

troll <- readRDS("TROLL_Table1.rds") %>%
    select(-LineNumber) %>%
    mutate(SpeciesID = as.factor(SpeciesID),
           SpeciesID = as.numeric(SpeciesID))

troll <- troll %>%
    mutate(Model = as.character(Model),
           Ninitial = as.factor(Ninitial),
           Rep = as.numeric(Rep),
           SeedRain = as.factor(SeedRain),
           SpeciesID = as.character(SpeciesID),
           Stage = as.factor(Stage),
           Year = as.numeric(Year),
           Biomass = as.numeric(Biomass),
           Productivity = as.numeric(Productivity))

troll <- troll %>%
    mutate(Stage = recode(Stage,
                          assembly = "metacommunity",
                          disassembly = "isolation"))

troll <- troll %>%
    filter(SeedRain %in% c(100))

troll$SeedRain <- droplevels(troll$SeedRain)

troll <- troll %>%
    filter(Year %in% c(100, 200))

p <- ggplot(troll,
       aes(x = Biomass,
           y = Productivity,
           color = Stage)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(facets = vars(Ninitial), scales = "free") +
    theme_bw() +
    theme(aspect.ratio = 1); p

cowplot::save_plot(p, filename = "~/Desktop/TROLL1.pdf", ncol = 3, nrow = 3)


# -------------------------------------------------------------------------

library(data.table)
library(tidyverse)
library(scales)

setwd("/Users/Theodore/Documents/Dropbox/Work/Projects/SimNet/scripts")
setwd("../data/raw")

troll <- readRDS("TROLL_Table1.rds") %>%
    select(-LineNumber) %>%
    mutate(SpeciesID = as.factor(SpeciesID),
           SpeciesID = as.numeric(SpeciesID))

troll <- troll %>%
    mutate(Model = as.character(Model),
           Ninitial = as.factor(Ninitial),
           Rep = as.numeric(Rep),
           SeedRain = as.factor(SeedRain),
           SpeciesID = as.character(SpeciesID),
           Stage = as.factor(Stage),
           Year = as.numeric(Year),
           Biomass = as.numeric(Biomass),
           Productivity = as.numeric(Productivity))

troll <- troll %>%
    mutate(Stage = recode(Stage,
                          assembly = "metacommunity",
                          disassembly = "isolation"))

troll <- troll %>%
    filter(SeedRain %in% c(0, 100))

troll$SeedRain <- droplevels(troll$SeedRain)

troll <- troll %>%
    filter(Year %in% c(100, 200),
           Ninitial == 1)

p <- ggplot(troll,
            aes(x = Biomass,
                y = Productivity,
                color = SeedRain)) +
    geom_point() +
    geom_smooth(method = "lm") +
    facet_wrap(facets = vars(Year), scales = "free") +
    theme_bw() +
    theme(aspect.ratio = 1); p

cowplot::save_plot(p, filename = "~/Desktop/TROLL2.pdf", ncol = 1, nrow = 2)
