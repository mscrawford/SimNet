library(FD)
library(tidyverse)

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


# ---------------------------------------------------------------------------------------------
# Compose SimNet communities using FD packages

d <- read.csv("selectWeiss.txt", sep = " ")

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

write.table(new.community, file = "SimNetCommunity.txt", sep = " ")

# Composing SimNet communities using Kmeans
# I chose this approach because of its simplicity.

data <- read.csv("selectWeiss.txt", sep = " ")

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

write.table(data, file = "SimNetCommunity.txt", sep = " ", quote = F, row.names = F)


# ---------------------------------------------------------------------------------------------
# Infer appropriate seed addition levels

library(data.table)

data.dir = "./IBC-grass_SimNet_seedRainDerivation_data/"
trait.dir = "."

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

data.dir = "./IBC-grass_SimNet_seedRainDerivation_data/"

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
