library(data.table)
library(tidyverse)
library(FD)

library(assertthat)

library(rstan)
library(brms)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SimNet_ReadModels.R")
setwd("../tmp")


# ---------------------------------------------------------------------------------------------
# Options -------------------------------------------------------------------------------------

calculate_FD = TRUE # Toggle whether functional diversity is incorporated (computationally expensive)

# Focal ecosystem function
# Y_VAL = "biomass"
Y_VAL = "productivity"


# Functions -----------------------------------------------------------------------------------

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

# Calculate functional diversity
FD_fun <- function(d_)
{
    d_ <- d_ %>%
        select(-Productivity)

    d_ <- d_ %>%
        filter(Biomass > 0) %>%
        arrange(SpeciesID)

    # if there is only one species remaining, it is considered maximally undispersed
    if (dim(d_)[1] == 1)
    {
        return(data.frame(FDis = 0))
    }

    d_.traits <- d_ %>%
        column_to_rownames("SpeciesID") %>%
        select(-Biomass)

    d_.community <- d_ %>%
        select(SpeciesID, Biomass) %>%
        spread(SpeciesID, Biomass)

    d_.FD <- dbFD(x = as.data.frame(d_.traits),
                  a = as.data.frame(d_.community),
                  w.abun = TRUE)

    return(data.frame(FDis = d_.FD$FDis))
}

# Calculate functional diversity on a per trait basis
FD_per_trait_fun <- function(d_)
{
    d_ <- d_ %>%
        select(-Productivity)

    d_ <- d_ %>%
        filter(Biomass > 0) %>%
        arrange(SpeciesID)

    d_.SpeciesID <- d_ %>%
        select(SpeciesID)

    d_.traits <- d_ %>%
        select(-Biomass,
               -SpeciesID)

    d_.community <- d_ %>%
        select(SpeciesID, Biomass) %>%
        spread(SpeciesID, Biomass)

    # If the trait did not vary along this axis, retain these traits ...
    d_.FD.noVariation <- keep(.x = d_.traits,
                              .p = function(x)
                              {
                                  return(length(unique(x)) == 1)
                              })

    # And make stubs that indicate they had no functional dispersion
    d_.FD.noVariation.simplified <- map_df(.x = d_.FD.noVariation,
                                           .f = function(x)
                                           {
                                               return(data.frame(FDis = 0))
                                           },
                                           .id = "trait")

    # Remove any specific traits that don't have variation in them, as they will break dbFD
    d_.FD.filtered <- purrr::discard(.x = d_.traits,
                                     .p = function(.x)
                                     {
                                         return(length(unique(.x)) == 1)
                                     })

    # Run dbFD on each trait within each community
    d_.FD <- map(.x = d_.FD.filtered,
                 .f = function(x)
                 {
                     return(dbFD(x = as.data.frame(x) %>%
                                     cbind(d_.SpeciesID) %>%
                                     column_to_rownames("SpeciesID"),
                                 a = as.data.frame(d_.community),
                                 w.abun = TRUE))
                 })

    # Combine the newly returned dbFD data.frames into one dataframe with "trait" as an identifier
    d_.FD.simplified <- map_df(.x = d_.FD,
                               .f = function(x)
                               {
                                   return(data.frame(FDis = x$FDis))
                               },
                               .id = "trait")

    return(bind_rows(d_.FD.noVariation.simplified, d_.FD.simplified))
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


# ---------------------------------------------------------------------------------------------
# The main datasets' models

d <- map(.x = models,
         .f = ~
             {
                 .x <- .x %>%
                     group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
                     nest()

                 .x <- .x %>%
                     mutate(FD = ifelse(calculate_FD, map(data, FD_fun), NA),
                            BEF = map(data, BEF_fun)) %>%
                     select(-data) %>%
                     unnest(c(FD, BEF))
             }) %>%
    bind_rows() %>%
    tbl_df()

d <- d %>%
    group_by(Model) %>%
    nest() %>%
    expand_grid(measure = c("Shannon", "Richness", "FDis"))

# When the measure is "FDis", filter data to remove Ninitial == 1
d <- d %>%
    mutate(data = ifelse(measure == "FDis",
                         map(.x = data,
                             .f = ~ .x %>% filter(Ninitial != 1)),
                         data))

# Run brms models
d <- d %>%
    mutate(mainEffect = pmap(.l = list(data, Model, measure),
                             .f = ~ brm_mainEffect_fun(..1, ..2, ..3)),
           treatmentEffect = pmap(.l = list(data, Model, measure),
                                  .f = ~ brm_treatmentEffect_fun(..1, ..2, ..3)))

# Save all the model results
saveRDS(d, file = paste0(Y_VAL, "_brms_models.rds"))


# -----------------------------------------------------------------------------------------------
# Trait-level FDis dataset ----------------------------------------------------------------------

d <- map(.x = models,
         .f = ~
             {
                 .x <- .x %>%
                     group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
                     nest()

                 .x <- .x %>%
                     mutate(FD = map(data, FD_per_trait_fun),
                            BEF = map(data, BEF_fun)) %>%
                     select(-data) %>%
                     unnest(c(FD, BEF))

                 .x <- .x %>% filter(Shannon > 1) # Only one species left, so the values will all be NA
             }) %>%
    bind_rows() %>%
    tbl_df()

d <- d %>%
    group_by(Model, trait) %>%
    nest() %>%
    expand_grid(measure = c("FDis"))

# When the measure is "FDis", filter data to remove Ninitial == 1
d <- d %>%
    mutate(data = ifelse(measure == "FDis",
                         map(.x = data,
                             .f = ~ .x %>% filter(Ninitial != 1)),
                         data))

# Run brms models
d <- d %>%
    mutate(mainEffect = pmap(.l = list(data, Model, measure),
                             .f = ~ brm_mainEffect_fun(..1, ..2, ..3)),
           treatmentEffect = pmap(.l = list(data, Model, measure),
                                  .f = ~ brm_treatmentEffect_fun(..1, ..2, ..3)))

# Save all the model results
saveRDS(d, file = "brms_models_FDis_perTrait.rds")
