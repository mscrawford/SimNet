library(tidyverse)
library(cowplot)
library(brms)
library(rsample)
library(lmodel2)
library(infer)

bootstrap_N = 2500


# ---------------------------------------------------------------------------------------------
# Focal ecosystem function

Y_VAL = "biomass"
# Y_VAL = "productivity"


# ---------------------------------------------------------------------------------------------
# Bootstrap correlation coefficients between monoculture biomass and competitive ability (for x-axis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SimNet_ReadModels.R")
setwd("../tmp")

pearsons <- model_runs %>%
    bind_rows() %>%
    filter(Ninitial != 64)

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
                      .f = ~
                          {
                              bootstraps <- .x %>% rsample::bootstraps(times = bootstrap_N)

                              bootstraps <- bootstraps %>%
                                  mutate(r = map_dbl(.x = splits,
                                                     .f = ~ {
                                                         .x <- as.data.frame(.x)
                                                         return(cor(.x[[paste0("monoculture.", Y_VAL)]],
                                                                    .x[[str_to_title(Y_VAL)]]))
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
# estimates <- readRDS(file = paste0("../tmp/", Y_VAL, "_brms_models.rds"))
setwd("../tmp")

# Main effects
main.effects <- estimates %>%
    select(Model, measure, mainEffect)

bootstrap.main.estimate <- main.effects %>%
    mutate(bootstrap.estimate = map(.x = mainEffect,
                                    .f = ~
                                        {
                                            id = seq(1, bootstrap_N)

                                            estimates_ <- posterior_samples(.x, paste0(":", measure)) %>%
                                                sample_n(bootstrap_N,
                                                         replace = TRUE)

                                            estimates_ <- cbind(estimates_, id)

                                            estimates_ <- estimates_ %>%
                                                pivot_longer(cols = -id,
                                                             names_to = "term",
                                                             values_to = "estimate") %>%
                                                mutate(Stage = ifelse(str_detect(term, "metacommunity"), "metacommunity", "isolation"),
                                                       Stage = factor(Stage, levels = c("metacommunity", "isolation")))

                                            estimates_ <- estimates_ %>%
                                                select(Stage, id, estimate)

                                            return(estimates_)
                                        })) %>%
    unnest(bootstrap.estimate)

bootstrap.main.effect <- inner_join(bootstrap.pearsons %>% filter(Ninitial == 32),
                                    bootstrap.main.estimate) %>%
    ungroup() %>%
    select(measure, Model, Stage, id, r, estimate)


# Treatment effects
treatment.effects <- estimates %>%
    select(measure, Model,treatmentEffect)

bootstrap.treatment.estimates <- treatment.effects %>%
    mutate(bootstrap.estimates = map(.x = treatmentEffect,
                                     .f = ~
                                         {
                                             id = seq(1, bootstrap_N)

                                             estimates_ <- posterior_samples(.x, paste0(":", measure)) %>%
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
    select(measure, Model, Stage, Ninitial, id, r, estimate)


# ---------------------------------------------------------------------------------------------
# Run bootstrapped Type II Regression (reduced major axis regression)

# Main effect
bootstrap.main.effect <- bootstrap.main.effect %>%
    group_by(measure, Stage, id) %>%
    nest() %>%
    mutate(model = map(.x = data,
                       .f = ~ lmodel2::lmodel2(.x$estimate ~ .x$r)))

bootstrap.main.effect <- bootstrap.main.effect %>%
    mutate(regression.results = map(.x = model,
                                    .f = ~ .x$regression.results[3,])) %>% # SMA Regression
    unnest(regression.results)

# Treatment effect
bootstrap.treatment.effect <- bootstrap.treatment.effect %>%
    group_by(measure, Stage, Ninitial, id) %>%
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
    select(measure, Stage, data) %>%
    unnest(data)

main.effect.conf.int <- bootstrap.main.effect %>%
    ungroup() %>%
    select(measure, Stage, id, Intercept, Slope) %>%
    group_by(measure, Stage) %>%
    arrange(Intercept) %>%
    mutate(ci_group = dplyr::ntile(Intercept, n = 25))

main.effect.conf.int <- main.effect.conf.int %>%
    group_by(measure, Stage, ci_group) %>%
    nest() %>%
    mutate(intercept.mean = map_dbl(.x = data, .f = ~ mean(.x$Intercept))) %>%
    mutate(ci = map(.x = data,
                    .f = ~
                        {
                            .x %>%
                                infer::specify(response = Slope) %>%
                                infer::generate(reps = 1000, type = "bootstrap") %>%
                                infer::calculate(stat = "mean") %>%
                                infer::get_confidence_interval(level = 0.95, type = "percentile")
                        })) %>%
    unnest(ci)

main.effect.conf.int <- main.effect.conf.int %>%
    group_by(measure, Stage) %>%
    select(-ci_group, -data) %>%
    filter(intercept.mean == min(intercept.mean) | intercept.mean == max(intercept.mean)) %>%
    mutate(X = list(seq(min(main.points$r) - 0.05, 1, length.out = 100))) %>%
    unnest(c(X)) %>%
    mutate(Y_min = `2.5%` * X + intercept.mean,
           Y_max = `97.5%` * X + intercept.mean) %>%
    select(-intercept.mean, -`2.5%`, -`97.5%`) %>%
    group_by(measure, Stage, X) %>%
    summarise(Y_min = min(Y_min), Y_max = max(Y_max))

main.effect.no.slope.lines <- main.points %>%
    group_by(measure, Stage) %>%
    summarise(mean.estimate = mean(estimate))

walk(.x = unique(estimates$measure),
     .f = ~
         {
             main.points__ <- main.points %>%
                 filter(measure == .x)

             main.effect.conf.int__ <- main.effect.conf.int %>%
                 filter(measure == .x)

             main.effect.no.slope.lines__ <- main.effect.no.slope.lines %>%
                 filter(measure == .x)

             a <- ggplot() +
                 geom_point(data = main.points__,
                            mapping = aes(x = r,
                                          y = estimate,
                                          color = Model),
                            alpha = 0.20) +
                 geom_ribbon(data = main.effect.conf.int__,
                             mapping = aes(x = X,
                                           ymin = Y_min,
                                           ymax = Y_max),
                             fill = "darkgrey",
                             alpha = 0.66) +
                 geom_hline(data = main.effect.no.slope.lines__,
                            mapping = aes(yintercept = mean.estimate),
                            linetype = 3) +
                 facet_grid(cols = vars(Stage)) +
                 labs(x = paste0("Pearson's r, monoculture ", Y_VAL, " and ", Y_VAL, " in 32-species mixture"),
                      y = paste0("Effect of realized diversity on ", Y_VAL)) +
                 scale_color_viridis_d(labels = c("C2018", "R2006", "M2009", "T2013", "R2020", "M2017")) +
                 guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                 theme_bw(11) +
                 theme(aspect.ratio = 1); a

             cowplot::save_plot(a,
                                filename = paste0(Y_VAL, "_acrossTreatmentEffects_", .x, ".png"),
                                ncol = 2,
                                nrow = 1,
                                base_asp = 1.1)

         }
)

# Treatment effect
treatment.points <- bootstrap.treatment.effect %>%
    ungroup() %>%
    select(measure, Stage, Ninitial, data) %>%
    unnest(data)

treatment.effect.conf.int <- bootstrap.treatment.effect %>%
    ungroup() %>%
    select(measure, Ninitial, Stage, id, Intercept, Slope) %>%
    group_by(measure, Ninitial, Stage) %>%
    arrange(Intercept) %>%
    mutate(ci_group = dplyr::ntile(Intercept, n = 25))

treatment.effect.conf.int <- treatment.effect.conf.int %>%
    group_by(measure, Ninitial, Stage, ci_group) %>%
    nest() %>%
    mutate(intercept.mean = map_dbl(.x = data, .f = ~ mean(.x$Intercept))) %>%
    mutate(ci = map(.x = data,
                    .f = ~ {
                        .x %>%
                            infer::specify(response = Slope) %>%
                            infer::generate(reps = 1000, type = "bootstrap") %>%
                            infer::calculate(stat = "mean") %>%
                            infer::get_confidence_interval(level = 0.95, type = "percentile")
                    })) %>%
    unnest(ci)

treatment.effect.conf.int <- treatment.effect.conf.int %>%
    group_by(measure, Stage, Ninitial,) %>%
    select(-ci_group, -data) %>%
    filter(intercept.mean == min(intercept.mean) | intercept.mean == max(intercept.mean)) %>%
    mutate(X = list(seq(min(treatment.points$r) - 0.05, 1, length.out = 100))) %>%
    unnest(c(X)) %>%
    mutate(Y_min = `2.5%` * X + intercept.mean,
           Y_max = `97.5%` * X + intercept.mean) %>%
    select(-intercept.mean, -`2.5%`, -`97.5%`) %>%
    group_by(measure, Stage, Ninitial, X) %>%
    summarise(Y_min = min(Y_min), Y_max = max(Y_max))

treatment.effect.no.slope.lines <- treatment.points %>%
    group_by(measure, Stage, Ninitial) %>%
    summarise(mean.estimate = mean(estimate))

walk(.x = unique(estimates$measure),
     .f = ~
         {

             treatment.points__ <- treatment.points %>%
                 filter(measure == .x)

             treatment.effect.conf.int__ <- treatment.effect.conf.int %>%
                 filter(measure == .x)

             treatment.effect.no.slope.lines__ <- treatment.effect.no.slope.lines %>%
                 filter(measure == .x)

             b <- ggplot() +
                 geom_point(data = treatment.points__,
                            mapping = aes(x = r,
                                          y = estimate,
                                          color = Model),
                            alpha = 0.20) +
                 geom_ribbon(data = treatment.effect.conf.int__,
                             mapping = aes(x = X,
                                           ymin = Y_min,
                                           ymax = Y_max),
                             fill = "darkgrey",
                             alpha = 0.66) +
                 geom_hline(data = treatment.effect.no.slope.lines__,
                            mapping = aes(yintercept = mean.estimate),
                            linetype = 3) +
                 facet_grid(cols = vars(Stage),
                            rows = vars(Ninitial)) +
                 labs(x = paste0("Pearson's r, monoculture ", Y_VAL, " and ", Y_VAL, " in mixture"),
                      y = paste0("Effect of realized diversity on ", Y_VAL)) +
                 scale_color_viridis_d(labels = c("C2018", "R2006", "M2009", "T2013", "R2020", "M2017")) +
                 guides(colour = guide_legend(override.aes = list(alpha = 1))) +
                 theme_bw(14) +
                 theme(aspect.ratio = 1); b

             cowplot::save_plot(b,
                                filename = paste0(Y_VAL, "_withinTreatmentEffects_", .x, ".png"),
                                ncol = 2,
                                nrow = 5,
                                base_asp = 1.1)

         }
)
