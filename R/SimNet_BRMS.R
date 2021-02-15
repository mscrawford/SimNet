library(data.table)
library(tidyverse)
library(assertthat)

library(FD)

library(rstan)
library(brms)


# Options -----------------------------------------------------------------

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# -------------------------------------------------------------------------
# Functions


# BRMS functions
brm_acrossEffect_fun <- function(d_, model, measure)
{
    cat("Model: ", model,
        " | Biodiversity measure: ", measure,
        " | Function: ", Y_VAL,
        " | acrossEffect\n",
        sep = "")

    assert_that(is.factor(d_$Stage), msg = "Stage is not a factor.")

    m_ <- brm(formula = paste0(Y_VAL, " ~ -1 + Stage + Stage:", measure),
              data = d_)

    return(m_)
}


brm_withinEffect_fun <- function(d_, model, measure)
{
    cat("Model: ", model,
        " | Biodiversity measure: ", measure,
        " | Function: ", Y_VAL,
        " | withinEffect\n",
        sep = "")

    d_ <- d_ %>% filter(Ninitial %in% c(2, 4, 8, 16, 32)) # monocultures and 64-species plots have no variation

    assert_that(is.factor(d_$Stage), msg = "Stage is not a factor.")
    assert_that(is.factor(d_$Ninitial), msg = "Ninitial is not a factor.")

    m_ <- brm(formula = paste0(Y_VAL, " ~ -1 + Ninitial:Stage + Ninitial:Stage:", measure),
              data = d_)

    return(m_)
}


# Calculate Shannon diversity, richness, total biomass, and productivity
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
                  w.abun = TRUE,
                  calc.FRic = FALSE,
                  calc.FDiv = FALSE,
                  messages = FALSE)

    if (file.exists("vert.txt")) {
        file.remove("vert.txt")
    }

    return(data.frame(FDis = d_.FD$FDis))
}


# -------------------------------------------------------------------------
# Calculate biomass and productivity with the models

# Read prepared simulation datasets
if (!exists("models"))
{
    source(paste0(scripts_dir, "/SimNet_ReadModels.R"))
}

d <- map(.x = models,
         .f = ~
             {
                 .x <- .x %>%
                     group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
                     nest()

                 .x <- .x %>%
                     mutate(FD = ifelse("FDis" %in% MEASURES, map(data, FD_fun), NA),
                            BEF = map(data, BEF_fun)) %>%
                     select(-data) %>%
                     unnest(c(FD, BEF))

                 .x <- .x %>% # Standardizing all of the communities to the highest community biomass within each model
                     ungroup() %>%
                     mutate(biomass = (biomass / max(biomass)) * 100,
                            productivity = (productivity / max(productivity)) * 100)
             }) %>%
    bind_rows() %>%
    as_tibble()

d <- d %>%
    group_by(Model) %>%
    nest() %>%
    expand_grid(measure = MEASURES)

# When the measure is "FDis", filter data to remove Ninitial == 1
d <- d %>%
    mutate(data = ifelse(measure == "FDis",
                         map(.x = data,
                             .f = ~ .x %>% filter(Ninitial != 1)),
                         data))

if (SAVE_CACHE)
{
    saveRDS(d, file = paste0(tmp_dir, "/cache/", Y_VAL, "_pre_brms_models_CACHED.rds"))
}


# -------------------------------------------------------------------------
# Calculate brms models

if (READ_CACHE)
{
    assign("d", readRDS(file = paste0(tmp_dir, "/cache/", Y_VAL, "_pre_brms_models_CACHED.rds")), envir = .GlobalEnv)
}

estimates <- d %>%
    mutate(acrossEffect = pmap(.l = list(data, Model, measure),
                               .f = ~ brm_acrossEffect_fun(..1, ..2, ..3)),
           withinEffect = pmap(.l = list(data, Model, measure),
                               .f = ~ brm_withinEffect_fun(..1, ..2, ..3)))


# -------------------------------------------------------------------------
# Save the brms model results

saveRDS(estimates, file = paste0(tmp_dir, "/", Y_VAL, "_brms_models.rds"))

if (SAVE_CACHE)
{
    saveRDS(estimates, file = paste0(tmp_dir, "/cache/",  Y_VAL, "_brms_models_CACHED.rds"))
}
