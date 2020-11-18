library(tidyverse)
library(ggthemes)
library(cowplot)
library(brms)
library(rsample)
library(lmodel2)
library(infer)

bootstrap_N = 2500


# ---------------------------------------------------------------------------------------------
# Focal ecosystem function

# Y_VAL = "biomass"
# Y_VAL = "productivity"
Y_VAL = "relative_yield_total"


# ---------------------------------------------------------------------------------------------
# Bootstrap correlation coefficients between monoculture biomass and competitive ability (for x-axis)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("SimNet_ReadModels.R")
setwd("../tmp")

pearsons <- model_runs %>%
    bind_rows()

pearsons.mono <- pearsons %>%
    filter(Ninitial == 1) %>%
    rename(monoculture.biomass = Biomass,
           monoculture.productivity = Productivity) %>%
    ungroup() %>%
    select(-Ninitial, -Rep)

pearsons <- inner_join(pearsons.mono, pearsons) %>%
    filter(Ninitial == 32) %>%
    select(Model, Rep, SpeciesID, Stage,
           Biomass, Productivity,
           monoculture.biomass, monoculture.productivity)

bootstrap.pearsons <- pearsons %>%
    group_by(Model, Stage) %>%
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
    select(Stage, id, Model, r)


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
    mutate(bootstrap.estimate = map2(.x = mainEffect,
                                     .y = measure,
                                     .f = ~
                                         {
                                             id = seq(1, bootstrap_N)

                                             estimates_ <- posterior_samples(.x, paste0(":", .y)) %>%
                                                 sample_n(bootstrap_N,
                                                          replace = TRUE)

                                             estimates_ <- cbind(estimates_, id)

                                             estimates_ <- estimates_ %>%
                                                 pivot_longer(cols = -id,
                                                              names_to = "parameter",
                                                              values_to = "estimate") %>%
                                                 mutate(Stage = ifelse(str_detect(parameter, "Withseedrain"), "With seed rain", "Without seed rain"),
                                                        Stage = as.factor(Stage))

                                             estimates_ <- estimates_ %>%
                                                 select(Stage, id, estimate)

                                             return(estimates_)
                                         })) %>%
    unnest(bootstrap.estimate)

bootstrap.main.effect <- inner_join(bootstrap.pearsons,
                                    bootstrap.main.estimate) %>%
    ungroup() %>%
    select(measure, Model, Stage, id, r, estimate)


# Treatment effects
treatment.effects <- estimates %>%
    select(measure, Model, treatmentEffect)

bootstrap.treatment.estimates <- treatment.effects %>%
    mutate(bootstrap.estimates = map2(.x = treatmentEffect,
                                      .y = measure,
                                      .f = ~
                                          {
                                              id = seq(1, bootstrap_N)

                                              estimates_ <- posterior_samples(.x, paste0(":", .y)) %>%
                                                  sample_n(bootstrap_N,
                                                           replace = TRUE)

                                              estimates_ <- cbind(estimates_, id)

                                              estimates_ <- estimates_ %>%
                                                  pivot_longer(cols = -id,
                                                               names_to = "parameter",
                                                               values_to = "estimate") %>%
                                                  mutate(Stage = ifelse(str_detect(parameter, "Withseedrain"), "With seed rain", "Without seed rain"),
                                                         Stage = as.factor(Stage)) %>%
                                                  mutate(Ninitial = str_extract(parameter, "\\d+"),
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

main.points.prediction.interval <- bootstrap.main.effect %>%
    group_by(measure, Stage) %>%
    select(measure, Stage, Intercept, Slope) %>%
    nest() %>%
    mutate(predictions = map(.x = data,
                             .f = ~ {
                                 lines <- data.frame(Slope = .x$Slope, Intercept = .x$Intercept)

                                 predictions <- expand_grid(lines, X = seq(-0.4, 1, length.out = 100)) %>%
                                     mutate(Y = Slope * X + Intercept)

                             })) %>%
    select(-data) %>%
    unnest(predictions)

main.points.prediction.interval <- main.points.prediction.interval %>%
    group_by(measure, Stage, X) %>%
    select(measure, Stage, X, Y) %>%
    nest() %>%
    mutate(quantiles = map(.x = data,
                           .f = ~ {
                               data.frame(ci_upper = quantile(.x$Y, 0.975),
                                          ci_mid = quantile(.x$Y, 0.50),
                                          ci_lower = quantile(.x$Y, 0.025))
                           })) %>%
    select(-data) %>%
    unnest(quantiles)

main.effect.slope.CI <- bootstrap.main.effect %>%
    group_by(measure, Stage) %>%
    nest() %>%
    mutate(ci = map(.x = data,
                    .f = ~
                        {
                            .x %>%
                                infer::specify(response = Slope) %>%
                                infer::generate(reps = 1000, type = "bootstrap") %>%
                                infer::calculate(stat = "mean") %>%
                                infer::get_confidence_interval(level = 0.95, type = "percentile")
                        })) %>%
    select(-data) %>%
    unnest(ci)

walk(.x = unique(estimates$measure),
     .f = ~
         {
             main.points__ <- main.points %>%
                 filter(measure == .x)

             main.points.prediction.interval__ <- main.points.prediction.interval %>%
                 filter(measure == .x)

             main.effect.slope.CI__ <- main.effect.slope.CI %>%
                 filter(measure == .x) %>%
                 mutate(lower_ci = round(lower_ci, digits = 2),
                        upper_ci = round(upper_ci, digits = 2))

             x__ <- main.points.prediction.interval__ %>%
                 ungroup() %>%
                 filter(X == max(X)) %>%
                 distinct(X) %>%
                 pull(X)

             y__ <- main.points.prediction.interval__ %>%
                 ungroup() %>%
                 filter(ci_lower == min(ci_lower)) %>%
                 distinct(ci_lower) %>%
                 pull(ci_lower)

             main.effect.slope.CI__ <- main.effect.slope.CI__ %>%
                 mutate(x_max = x__,
                        y_min = y__)

             a <- ggplot() +
                 geom_hline(yintercept = 0, color = "darkgray") +
                 geom_vline(xintercept = 0, color = "darkgray") +
                 geom_point(data = main.points__,
                            mapping = aes(x = r,
                                          y = estimate,
                                          color = Model),
                            alpha = 0.10) +
                 scale_color_brewer(palette = "Dark2") +
                 geom_ribbon(data = main.points.prediction.interval__,
                             mapping = aes(x = X,
                                           ymin = ci_lower,
                                           ymax = ci_upper),
                             fill = "darkgrey",
                             alpha = 0.66) +
                 scale_fill_brewer(palette = "Greys") +
                 geom_line(data = main.points.prediction.interval__,
                           mapping = aes(x = X,
                                         y = ci_mid),
                           color = "#3E5EF1") +
                 geom_text(data = main.effect.slope.CI__,
                           mapping = aes(x = x_max,
                                         y = y_min,
                                         label = paste0("95% CI of slope: [", lower_ci, ", ", upper_ci, "]")),
                           hjust = 1) +
                 facet_grid(cols = vars(Stage)) +
                 labs(x = "Function-dominance correlation",
                      y = "Slope of BEF relationship") +
                 guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
                 theme_few(15) +
                 theme(aspect.ratio = 0.618)

             cowplot::save_plot(a,
                                filename = paste0(Y_VAL, "_acrossTreatmentEffects_", .x, ".png"),
                                ncol = 2,
                                nrow = 1,
                                base_asp = 1.75)

         }
)


# Treatment effect
treatment.points <- bootstrap.treatment.effect %>%
    ungroup() %>%
    select(measure, Stage, Ninitial, data) %>%
    unnest(data)

treatment.points.prediction.interval <- bootstrap.treatment.effect %>%
    group_by(measure, Ninitial, Stage) %>%
    select(measure, Ninitial, Stage, Intercept, Slope) %>%
    nest() %>%
    mutate(predictions = map(.x = data,
                             .f = ~ {
                                 lines <- data.frame(Slope = .x$Slope, Intercept = .x$Intercept)

                                 predictions <- expand_grid(lines, X = seq(-0.4, 1, length.out = 100)) %>%
                                     mutate(Y = Slope * X + Intercept)
                             })) %>%
    select(-data) %>%
    unnest(predictions)

treatment.points.prediction.interval <- treatment.points.prediction.interval %>%
    select(measure, Ninitial, Stage, X, Y) %>%
    group_by(measure, Ninitial, Stage, X) %>%
    nest() %>%
    mutate(quantiles = map(.x = data,
                           .f = ~ {
                               data.frame(ci_upper = quantile(.x$Y, 0.975),
                                          ci_mid = quantile(.x$Y, 0.50),
                                          ci_lower = quantile(.x$Y, 0.025))
                           })) %>%
    select(-data) %>%
    unnest(quantiles)

treatment.effect.slope.CI <- bootstrap.treatment.effect %>%
    group_by(measure, Ninitial, Stage) %>%
    nest() %>%
    mutate(ci = map(.x = data,
                    .f = ~
                        {
                            .x %>%
                                infer::specify(response = Slope) %>%
                                infer::generate(reps = 1000, type = "bootstrap") %>%
                                infer::calculate(stat = "mean") %>%
                                infer::get_confidence_interval(level = 0.95, type = "percentile")
                        })) %>%
    select(-data) %>%
    unnest(ci)


walk(.x = unique(estimates$measure),
     .f = ~
         {

             treatment.points__ <- treatment.points %>%
                 filter(measure == .x)

             treatment.points.prediction.interval__ <- treatment.points.prediction.interval %>%
                 filter(measure == .x)

             treatment.effect.slope.CI__ <- treatment.effect.slope.CI %>%
                 filter(measure == .x) %>%
                 mutate(lower_ci = round(lower_ci, digits = 2),
                        upper_ci = round(upper_ci, digits = 2))

            x__ <- treatment.points.prediction.interval__ %>%
                ungroup() %>%
                filter(X == max(X)) %>%
                distinct(X) %>%
                pull(X)

            y__ <- treatment.points__ %>%
                ungroup() %>%
                filter(estimate == min(estimate)) %>%
                distinct(estimate) %>%
                pull(estimate)

            treatment.effect.slope.CI__ <- treatment.effect.slope.CI__ %>%
                mutate(x_max = x__,
                       y_min = y__)

             p <- ggplot() +
                 geom_hline(yintercept = 0, color = "darkgray") +
                 geom_vline(xintercept = 0, color = "darkgray") +
                 geom_point(data = treatment.points__,
                            mapping = aes(x = r,
                                          y = estimate,
                                          color = Model),
                            alpha = 0.10) +
                 geom_ribbon(data = treatment.points.prediction.interval__,
                             mapping = aes(x = X,
                                           ymin = ci_lower,
                                           ymax = ci_upper),
                             fill = "darkgrey",
                             alpha = 0.66) +
                 geom_line(data = treatment.points.prediction.interval__,
                           mapping = aes(x = X,
                                         y = ci_mid),
                           color = "#9f140f",
                           size = 1) +
                 geom_text(data = treatment.effect.slope.CI__,
                           mapping = aes(x = x_max,
                                         y = y_min,
                                         label = paste0("95% CI of slope: [", lower_ci, ", ", upper_ci, "]")),
                           hjust = 1) +
                 facet_grid(cols = vars(Stage),
                            rows = vars(Ninitial)) +
                 scale_color_brewer(palette = "Dark2") +
                 scale_fill_brewer(palette = "Greys") +
                 labs(x = "Function-dominance correlation",
                      y = "Slope of BEF relationship") +
                 guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
                 theme_few(24) +
                 theme(aspect.ratio = 0.618,
                       legend.key.size = unit(0.75, "cm"))

             cowplot::save_plot(p,
                                filename = paste0(Y_VAL, "_withinTreatmentEffects_", .x, ".png"),
                                ncol = 2,
                                nrow = 5,
                                base_asp = 1.75)

         }
)
