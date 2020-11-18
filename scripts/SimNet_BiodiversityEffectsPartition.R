library(assertthat)
library(rstan)
library(brms)
library(data.table)
library(tidyverse)

bootstrap_N = 2500

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SimNet_ReadModels.R")
setwd("../tmp")

Y_VAL = "relative_yield_total"

# BRMS functions
brm_mainEffect_fun <- function(d_, model, measure)
{
    cat("Model: ", model, " | Measure: ", measure, " | mainEffect\n", sep = "")

    assert_that(is.factor(d_$Stage), msg = "Stage is not a factor.")

    m_ <- brm(formula = paste(Y_VAL, " ~ -1 +
                                            Stage +
                                            Stage:", measure,
                              sep = ""),
              data = d_)

    return(m_)
}

brm_treatmentEffect_fun <- function(d_, model, measure)
{
    cat("Model: ", model, " | Measure: ", measure, " | treatmentEffect\n", sep = "")

    d_ <- d_ %>% filter(Ninitial %in% c(2, 4, 8, 16, 32)) # monocultures and 64-species plots have no variation

    assert_that(is.factor(d_$Stage), msg = "Stage is not a factor.")
    assert_that(is.factor(d_$Ninitial), msg = "Ninitial is not a factor.")

    m_ <- brm(formula = paste(Y_VAL, " ~ -1 +
                                            Ninitial:Stage +
                                            Ninitial:Stage:", measure,
                              sep = ""),
              data = d_)

    return(m_)
}

# Calculate Shannon diversity, richness, and total biomass
BEF_fun <- function(d_)
{
    d_ <- d_ %>%
        filter(Biomass > 0) %>%
        mutate(p_i = Biomass / sum(Biomass)) %>%
        summarise(biomass = sum(Biomass),
                  productivity = sum(Productivity),
                  Shannon = 1 - sum(p_i * log(p_i)),
                  Richness = n())

    return(data.frame(biomass = d_$biomass,
                      productivity = d_$productivity,
                      Shannon = d_$Shannon,
                      Richness = d_$Richness))
}

BE_partition_fun <- function(d_)
{
    d_ <- d_ %>%
        mutate(Ninitial = as.numeric(as.character(Ninitial)))

    monocultures <- d_ %>%
        filter(Ninitial == 1) %>%
        select(-Ninitial, -Rep) %>%
        rename(M_i = Biomass)

    d_ <- inner_join(d_, monocultures) %>%
        filter(Ninitial > 1) %>%
        rename(Y_Oi = Biomass)

    d_ <- d_ %>%
        group_by(Ninitial, Rep) %>%
        mutate(RY_Ei = 1 / Ninitial)

    Y_O <- d_ %>%
        group_by(Ninitial, Rep) %>%
        summarise(Y_O = sum(Y_Oi))

    d_ <- d_ %>%
        group_by(Ninitial, Rep, SpeciesID) %>%
        mutate(RY_Oi = Y_Oi / M_i,
               Y_Ei = RY_Ei * M_i)

    Y_E <- d_ %>%
        group_by(Ninitial, Rep) %>%
        summarise(Y_E = sum(Y_Ei))

    delta.Y <- inner_join(Y_O, Y_E) %>%
        mutate(delta.Y = Y_O - Y_E)

    return(data.frame(Ninitial = as.factor(delta.Y$Ninitial),
                      Rep = delta.Y$Rep,
                      relative_yield_total = delta.Y$delta.Y))
}

d <- map(.x = models,
         .f = ~
             {
                 .x <- .x %>%
                     group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
                     nest()

                 .x <- .x %>%
                     mutate(BEF = map(data, BEF_fun)) %>%
                     select(-data) %>%
                     unnest(c(BEF))
             }) %>%
    bind_rows() %>%
    as_tibble()

BEF_partition_d <- map(.x = models,
                       .f = ~
                           {
                               .x <- .x %>%
                                   select(Model, Ninitial, Rep, SeedRain, Stage, Year, SpeciesID, Biomass)

                               .x <- .x %>%
                                   group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
                                   mutate(community.biomass = sum(Biomass))

                               .x <- .x %>%
                                   ungroup() %>%
                                   mutate(standardized.community.biomass = community.biomass / max(community.biomass) * 100) %>%
                                   mutate(Biomass = (Biomass / community.biomass) * standardized.community.biomass) %>%
                                   select(-community.biomass, -standardized.community.biomass)

                               .x <- .x %>%
                                   group_by(Model, SeedRain, Stage, Year) %>%
                                   nest()

                               .x <- .x %>%
                                   mutate(relative_yield_total = map(data, BE_partition_fun)) %>%
                                   select(-data) %>%
                                   unnest(c(relative_yield_total))
                           }) %>%
    bind_rows() %>%
    as_tibble()

d <- inner_join(d, BEF_partition_d)

d <- d %>%
    group_by(Model) %>%
    nest() %>%
    expand_grid(measure = c("Shannon"))

# Run brms models
d <- d %>%
    mutate(mainEffect = pmap(.l = list(data, Model, measure),
                             .f = ~ brm_mainEffect_fun(..1, ..2, ..3)),
           treatmentEffect = pmap(.l = list(data, Model, measure),
                                  .f = ~ brm_treatmentEffect_fun(..1, ..2, ..3)))

# Save all the model results
saveRDS(d, file = "relative_yield_total_brms_models.rds")
