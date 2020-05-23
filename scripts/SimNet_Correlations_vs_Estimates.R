library(tidyverse)
library(cowplot)
library(brms)
library(rsample)
library(lmodel2)
library(infer)

bootstrap_N = 2500


# ---------------------------------------------------------------------------------------------
# Focal ecosystem function

# Y_VAL = "biomass"
Y_VAL = "productivity"


# ---------------------------------------------------------------------------------------------
# Bootstrap correlation coefficients between monoculture biomass and competitive ability (for x-axis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SimNet_ReadModels.R")
setwd("../tmp")

pearsons <- model_runs %>%
    bind_rows() %>%
    filter(Ninitial != 64,
           Model != "IBC_grass.noNDD")

pearsons.mono <- pearsons %>%
    filter(Ninitial == 1) %>%
    rename(monoculture.biomass = Biomass,
           monoculture.productivity = Productivity) %>%
    ungroup() %>%
    select(-Ninitial, -Rep)

pearsons <- inner_join(pearsons.mono, pearsons) %>%
    filter(Ninitial != 1) %>%
    select(Model, Ninitial, Rep, SpeciesID, Stage,
           Biomass, Productivity,
           monoculture.biomass, monoculture.productivity)

# Remove both factor levels Ninitial = 64 and Ninitial = 1
pearsons$Ninitial <- droplevels(pearsons$Ninitial)

bootstrap.pearsons <- pearsons %>%
    group_by(Model, Ninitial, Stage) %>%
    nest() %>%
    mutate(cors = map(.x = data,
                      .f = ~ {
                          bootstraps <- .x %>% rsample::bootstraps(times = bootstrap_N)

                          bootstraps <- bootstraps %>%
                              mutate(r = map_dbl(.x = splits,
                                                 .f = ~ {
                                                     .x <- as.data.frame(.x)
                                                     return(cor(.x[[paste0("monoculture.", Y_VAL)]], .x[[str_to_title(Y_VAL)]]))
                                                 }),
                                     id = row_number())

                          return(bootstraps)
                      })) %>%
    unnest(cors) %>%
    select(Ninitial, Stage, id, Model, r)


# ---------------------------------------------------------------------------------------------
# Bootstrap estimates of community biomass vs. realized Shannon (for y-axis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
estimates <- readRDS(file = paste0("../data/brms_models/", Y_VAL, "_brms_models.rds"))
setwd("../tmp")

estimates <- estimates %>%
    filter(measure == "Shannon",
           Model != "IBC_grass.noNDD")

# Main effects
main.effects <- estimates %>%
    select(Model, measure, mainEffect)

bootstrap.main.estimate <- main.effects %>%
    group_by(Model) %>%
    mutate(bootstrap.estimate = map(.x = mainEffect,
                                    .f = ~ {
                                        id = seq(1, bootstrap_N)

                                        estimates_ <- posterior_samples(.x, ":Shannon") %>%
                                            sample_n(bootstrap_N,
                                                     replace = TRUE)

                                        estimates_ <- cbind(estimates_, id)

                                        estimates_ <- estimates_ %>%
                                            pivot_longer(cols = -id,
                                                         names_to = "term",
                                                         values_to = "estimate") %>%
                                            mutate(Stage = ifelse(str_detect(term, "metacommunity"), "metacommunity", "isolation"),
                                                   Stage = factor(Stage, levels = c("metacommunity", "isolation"))) %>%
                                            select(Stage, id, estimate)

                                        return(estimates_)
                                    })) %>%
    unnest(bootstrap.estimate)

bootstrap.main.effect <- inner_join(bootstrap.pearsons %>% filter(Ninitial == 32),
                                    bootstrap.main.estimate) %>%
    ungroup() %>%
    select(Model, Stage, id, r, estimate)


# Treatment effects
treatment.effects <- estimates %>%
    select(Model, measure, treatmentEffect)

bootstrap.treatment.estimates <- treatment.effects %>%
    group_by(Model) %>%
    mutate(bootstrap.estimates = map(.x = treatmentEffect,
                                     .f = ~ {
                                         id = seq(1, bootstrap_N)

                                         estimates_ <- posterior_samples(.x, ":Shannon") %>%
                                             sample_n(bootstrap_N,
                                                      replace = TRUE)

                                         estimates_ <- cbind(estimates_, id)

                                         estimates_ <- estimates_ %>%
                                             pivot_longer(cols = -id,
                                                          names_to = "term",
                                                          values_to = "estimate") %>%
                                             mutate(Stage = ifelse(str_detect(term, "metacommunity"), "metacommunity", "isolation"),
                                                    Stage = factor(Stage, levels = c("metacommunity", "isolation"))) %>%
                                             mutate(Ninitial = str_extract(term, "\\d+"),
                                                    Ninitial = factor(Ninitial, levels = c("2", "4", "8", "16", "32")))

                                         estimates_ <- estimates_ %>%
                                             select(Stage, Ninitial, id, estimate)

                                         return(estimates_)
                                     })) %>%
    unnest(bootstrap.estimates)

bootstrap.treatment.effect <- inner_join(bootstrap.pearsons, bootstrap.treatment.estimates) %>%
    ungroup() %>%
    select(Stage, Ninitial, id, Model, r, estimate)


# ---------------------------------------------------------------------------------------------
# Run bootstrapped Type II Regression (reduced major axis regression)

# Main effect
bootstrap.main.effect <- bootstrap.main.effect %>%
    group_by(Stage, id) %>%
    nest() %>%
    mutate(model = map(.x = data,
                       .f = ~ lmodel2::lmodel2(.x$estimate ~ .x$r)))

bootstrap.main.effect <- bootstrap.main.effect %>%
    mutate(regression.results = map(.x = model,
                                    .f = ~ .x$regression.results[3,])) %>% # SMA Regression
    unnest(regression.results)

# Treatment effect
bootstrap.treatment.effect <- bootstrap.treatment.effect %>%
    group_by(Stage, Ninitial, id) %>%
    nest() %>%
    mutate(model = map(.x = data,
                       .f = ~ lmodel2::lmodel2(.x$estimate ~ .x$r)))

bootstrap.treatment.effect <- bootstrap.treatment.effect %>%
    mutate(regression.results = map(.x = model,
                                    .f = ~ .x$regression.results[3,])) %>% # SMA Regression
    unnest(regression.results)


# ---------------------------------------------------------------------------------------------
# Generate plots

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("../tmp")

# Main effect
main.points <- bootstrap.main.effect %>%
    ungroup() %>%
    select(Stage, data) %>%
    unnest(data)

main.effect.conf.int <- bootstrap.main.effect %>%
    ungroup() %>%
    select(Stage, id, Intercept, Slope) %>%
    group_by(Stage) %>%
    arrange(Intercept) %>%
    mutate(ci_group = dplyr::ntile(Intercept, n = 25))

main.effect.conf.int <- main.effect.conf.int %>%
    group_by(Stage, ci_group) %>%
    nest() %>%
    mutate(mean = map(.x = data,
                      .f = ~ {
                          intercept.mean = mean(.x$Intercept)
                          return(data.frame(intercept.mean))
                      })) %>%
    mutate(ci = map(.x = data,
                    .f = ~ {
                        .x %>%
                            infer::specify(response = Slope) %>%
                            infer::generate(reps = 1000) %>%
                            infer::calculate(stat = "mean") %>%
                            infer::get_confidence_interval(level = 0.95, type = "percentile")
                    })) %>%
    unnest(c(mean, ci))

a <- ggplot() +
    geom_point(data = main.points,
               mapping = aes(x = r,
                             y = estimate,
                             color = Model),
               alpha = 0.25) +
    geom_abline(data = main.effect.conf.int,
                mapping = aes(intercept = intercept.mean,
                              slope = `2.5%`,
                              group = ci_group),
                color =  "darkred",
                alpha = 0.33) +
    geom_abline(data = main.effect.conf.int,
                mapping = aes(intercept = intercept.mean,
                              slope = `97.5%`,
                              group = ci_group),
                color =  "darkred",
                alpha = 0.33) +
    facet_grid(cols = vars(Stage)) +
    labs(x = paste0("Pearson's r, monoculture ", Y_VAL, " and ", Y_VAL, " in 32-species mixture"),
         y = paste0("Estimate, mixture diversity vs. total ", Y_VAL)) +
    scale_color_viridis_d() +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme_bw(12) +
    theme(aspect.ratio = 1); a

cowplot::save_plot(a,
                   filename = paste0(Y_VAL, "_mainEffects.png"),
                   ncol = 2,
                   nrow = 1,
                   base_asp = 1.1)


# Treatment effect
treatment.points <- bootstrap.treatment.effect %>%
    ungroup() %>%
    select(Stage, Ninitial, data) %>%
    unnest(data)

treatment.effect.conf.int <- bootstrap.treatment.effect %>%
    ungroup() %>%
    select(Ninitial, Stage, id, Intercept, Slope) %>%
    group_by(Ninitial, Stage) %>%
    arrange(Intercept) %>%
    mutate(ci_group = dplyr::ntile(Intercept, n = 25))

treatment.effect.conf.int <- treatment.effect.conf.int %>%
    group_by(Ninitial, Stage, ci_group) %>%
    nest() %>%
    mutate(mean = map(.x = data,
                      .f = ~ {
                          intercept.mean = mean(.x$Intercept)
                          return(data.frame(intercept.mean))
                      })) %>%
    mutate(ci = map(.x = data,
                    .f = ~ {
                        .x %>%
                            infer::specify(response = Slope) %>%
                            infer::generate(reps = 1000) %>%
                            infer::calculate(stat = "mean") %>%
                            infer::get_confidence_interval(level = 0.95, type = "percentile")
                    })) %>%
    unnest(c(mean, ci))

b <- ggplot() +
    geom_point(data = treatment.points,
               mapping = aes(x = r,
                             y = estimate,
                             color = Model),
               alpha = 0.25) +
    geom_abline(data = treatment.effect.conf.int,
                mapping = aes(intercept = intercept.mean,
                              slope = `2.5%`,
                              group = ci_group),
                color = "darkred",
                alpha = 0.33) +
    geom_abline(data = treatment.effect.conf.int,
                mapping = aes(intercept = intercept.mean,
                              slope = `97.5%`,
                              group = ci_group),
                color = "darkred",
                alpha = 0.33) +
    facet_grid(cols = vars(Stage),
               rows = vars(Ninitial)) +
    labs(x = paste0("Pearson's r, monoculture ", Y_VAL, " and ", Y_VAL, " in mixture"),
         y = paste0("Estimate, mixture diversity vs. total ", Y_VAL)) +
    scale_color_viridis_d() +
    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
    theme_bw(14) +
    theme(aspect.ratio = 1); b

cowplot::save_plot(b,
                   filename = paste0(Y_VAL, "_treatmentEffects.png"),
                   ncol = 2,
                   nrow = 5,
                   base_asp = 1.1)
