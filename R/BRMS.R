library(data.table)
library(tidyverse)
library(assertthat)

library(FD)

library(rstan)
library(brms)


# -------------------------------------------------------------------------
# Options

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


# -------------------------------------------------------------------------
# Functions

# BRMS functions
brm_acrossEffect_fun <- function(d, model, x_val)
{
    cat("Model: ", model,
        " | Biodiversity measure: ", x_val,
        " | Function: ", y_val,
        " | acrossEffect\n",
        sep = "")

    assert_that(is.factor(d$Stage), msg = "Stage is not a factor.")

    m_ <- brm(formula = paste0(y_val, " ~ -1 + Stage + Stage:", x_val),
              data = d)

    return(m_)
}


brm_withinEffect_fun <- function(d, model, x_val)
{
    cat("Model: ", model,
        " | Biodiversity measure: ", x_val,
        " | Function: ", y_val,
        " | withinEffect\n",
        sep = "")

    d <- d %>% filter(Ninitial %in% c(2, 4, 8, 16, 32)) # monocultures and 64-species plots have no variation

    assert_that(is.factor(d$Stage), msg = "Stage is not a factor.")
    assert_that(is.factor(d$Ninitial), msg = "Ninitial is not a factor.")

    m_ <- brm(formula = paste0(y_val, " ~ -1 + Ninitial:Stage + Ninitial:Stage:", x_val),
              data = d)

    return(m_)
}


# Calculate Shannon diversity, richness, total biomass, and productivity
BEF_fun <- function(d)
{
    d <- d %>%
        filter(Biomass > 0) %>%
        mutate(p_i = Biomass / sum(Biomass)) %>%
        summarise(biomass = sum(Biomass),
                  productivity = sum(Productivity),
                  Shannon = 1 - sum(p_i * log(p_i)),
                  Richness = n())

    return(data.frame(biomass = d$biomass,
                      productivity = d$productivity,
                      Shannon = d$Shannon,
                      Richness = d$Richness))
}


# Calculate functional diversity
FD_fun <- function(d)
{
    d <- d %>%
        select(-Productivity)

    d <- d %>%
        filter(Biomass > 0) %>%
        arrange(SpeciesID)

    # if there is only one species remaining, it is considered maximally undispersed
    if (dim(d)[1] == 1)
    {
        return(data.frame(FDis = 0))
    }

    d.traits <- d %>%
        column_to_rownames("SpeciesID") %>%
        select(-Biomass)

    d.community <- d %>%
        select(SpeciesID, Biomass) %>%
        spread(SpeciesID, Biomass)

    d.FD <- dbFD(x = as.data.frame(d.traits),
                 a = as.data.frame(d.community),
                 w.abun = TRUE,
                 calc.FRic = FALSE,
                 calc.FDiv = FALSE,
                 messages = FALSE)

    if (file.exists("vert.txt")) {
        file.remove("vert.txt")
    }

    return(data.frame(FDis = d.FD$FDis))
}


# -------------------------------------------------------------------------
# Calculate biomass and productivity with the models

# Read prepared simulation datasets
if (!exists("models"))
{
    source(paste0(scripts_dir, "/readModels.R"))
}

d <- map(.x = models,
         .f = ~
             {
                 .x <- .x %>%
                     group_by(Model, Ninitial, Rep, SeedRain, Stage, Year) %>%
                     nest()

                 .x <- .x %>%
                     mutate(FD = ifelse("FDis" %in% X_VALS, map(data, FD_fun), NA),
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
    expand_grid(x_val = X_VALS)

# When the x_val is "FDis", filter data to remove Ninitial == 1
d <- d %>%
    mutate(data = ifelse(x_val == "FDis",
                         map(.x = data,
                             .f = ~ .x %>% filter(Ninitial != 1)),
                         data))

if (SAVE_CACHE)
{
    saveRDS(d, file = paste0(tmp_dir, "/cache/", y_val, "_pre_brms_models_CACHED.rds"))
}


# -------------------------------------------------------------------------
# Calculate brms models

if (READ_CACHE)
{
    assign("d", readRDS(file = paste0(tmp_dir, "/cache/", y_val, "_pre_brms_models_CACHED.rds")), envir = .GlobalEnv)
}

estimates <- d %>%
    mutate(acrossEffect = pmap(.l = list(data, Model, x_val),
                               .f = ~ brm_acrossEffect_fun(..1, ..2, ..3)),
           withinEffect = pmap(.l = list(data, Model, x_val),
                               .f = ~ brm_withinEffect_fun(..1, ..2, ..3)))


# -------------------------------------------------------------------------
# Save the brms model results

saveRDS(estimates, file = paste0(tmp_dir, "/", y_val, "_brms_models.rds"))

if (SAVE_CACHE)
{
    saveRDS(estimates, file = paste0(tmp_dir, "/cache/",  y_val, "_brms_models_CACHED.rds"))
}
