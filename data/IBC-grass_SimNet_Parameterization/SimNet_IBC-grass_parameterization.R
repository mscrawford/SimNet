library(FD)
library(tidyverse)


# ---------------------------------------------------------------------------------------------
# Compose SimNet communities using FD packages

d <- read.csv("~/Desktop/selectWeiss.txt", sep = " ")

rownames(d) <- d$ID

d <- d %>%
    select(-Species, -ID)

max.FD.so.far <- 0
max.community <- NA

for (i in 1:1000)
{
    tryCatch (
        {
            d.t <- sample_n(d, 64)
            this.FD <- dbFD(d.t)$FRic
            cat("i: ", i, "  |  FRic: ", this.FD, "\n", sep = "")
            if (this.FD > max.FD.so.far)
            {
                max.FD.so.far <- this.FD
                max.community <- d.t
            }
        }, error = function(e){ cat("ERROR: ", conditionMessage(e), "\n" )}
    )
}

new.community <- inner_join(d, max.community)

write.table(new.community, file = "~/Desktop/SimNetCommunity.txt", sep = " ")

# Composing SimNet communities using Kmeans
# I chose this approach because of its simplicity.

data <- read.csv("~/Desktop/selectWeiss.txt", sep = " ")

clusters <- kmeans(data %>%
                       select(AllocSeed, LMR,
                              m0, MaxMass, mSeed, Dist,
                              pEstab,
                              Gmax, SLA,
                              palat, memo,
                              RAR, growth, mThres,
                              clonal, meanSpacerLength, sdSpacerlength,
                              Resshare, AllocSpacer,
                              mSpacer),
                   64)

data <- cbind(data, clusters$cluster)

data <- data %>%
    group_by(`clusters$cluster`) %>%
    sample_n(1)

data <- data %>%
    ungroup() %>%
    select(-`clusters$cluster`)

data <- data %>%
    arrange(ID)

write.table(data, file = "~/Desktop/SimNetCommunity.txt", sep = " ", quote = F, row.names = F)


# ---------------------------------------------------------------------------------------------
# Infer appropriate seed addition levels

source("/Users/Theodore/Documents/workspace/IBC-grass/util/scripts/Utils.R")

data.dir = "/Users/Theodore/Documents/Dropbox/Work/Projects/SimNet/Modelling_Experiment/scripts/IBC-grass-SimNet_Parameterization/IBC-grass_SimNet_seedAdditionDerivation_data/out"
trait.dir = "/Users/Theodore/Documents/Dropbox/Work/Projects/SimNet/Modelling_Experiment/scripts/IBC-grass-SimNet_Parameterization/IBC-grass_SimNet_seedAdditionDerivation_data/"

parameters <- read_data(data.dir, "parameter")

populations <- read_data(data.dir, "population")

traits <- read_data(trait.dir, "SimNetCommunity") %>%
    rename(PFT = ID,
           SeedMass = mSeed)

d <- combine_data(list(populations, traits %>% select(PFT, SeedMass)), c("PFT"))

d <- inner_join(d, parameters %>% select(SimID, Stabilization))

d <- d %>%
    filter(Year == 200) %>%
    mutate(nSeeds = Repro / SeedMass,
           seedsPerPlant = nSeeds / Pop,
           reproductiveBiomassPerRecruit = Repro / Recruits,
           seedsPerRecruit = nSeeds / Recruits)

d <- d %>%
    group_by(PFT, Stabilization) %>%
    summarise_each(funs = c(mean, sd), Repro, nSeeds, seedsPerPlant, reproductiveBiomassPerRecruit, seedsPerRecruit)

yearlyReproBiomass <- d %>%
    group_by(Stabilization) %>%
    summarise(mean(Repro_fn1))

yearlyReproBiomass

# ---------------------------------------------------------------------------------------------
# Check seed addition levels viability

source("/Users/Theodore/Documents/workspace/IBC-grass/util/scripts/Utils.R")

data.dir = "/Users/Theodore/Documents/workspace/IBC-grass/data/out"

parameters <- read_data(data.dir, "parameter")

populations <- read_data(data.dir, "population")

traits <- read_data(data.dir, "SimNetCommunity") %>%
    rename(PFT = ID,
           SeedMass = mSeed)

d <- combine_data(list(parameters, populations), c("SimID"))

d <- d %>% filter(Year %in% c(100, 200))

d <- d %>%
    rename(Biomass = Shootmass) %>%
    filter(Biomass > 0) %>%
    group_by(SeedInput, Year) %>%
    mutate(p_i = Biomass / sum(Biomass)) %>%
    summarise(total.biomass = sum(Biomass),
              Shannon = 1 - sum(p_i * log(p_i)),
              Richness = n())

p <- ggplot(d) +
    geom_point(aes(x = Shannon,
                   y = total.biomass,
                   color = as.factor(SeedInput))) +
    facet_grid(~ Year); p
