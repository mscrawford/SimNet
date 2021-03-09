library(tidyverse)
library(ggthemes)
library(cowplot)

library(brms)
library(rsample)
library(lmodel2)
library(infer)


# -------------------------------------------------------------------------
# Options

bootstrap_N = 2500


# -------------------------------------------------------------------------
# Bootstrap correlation coefficients between monoculture biomass and competitive ability (for x-axis)

if (!exists("model_runs"))
{
    source(paste0(scripts_dir, "/SimNet_ReadModels.R"))
}

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
                                                     .f = ~
                                                         {
                                                             .x <- as.data.frame(.x)
                                                             return(cor(.x[[paste0("monoculture.", y_val)]],
                                                                        .x[[str_to_title(y_val)]]))
                                                         }),
                                         id = row_number())

                              return(bootstraps)
                          })) %>%
    unnest(cors) %>%
    select(Stage, id, Model, r)

if (SAVE_CACHE)
{
    saveRDS(bootstrap.pearsons, file = paste0(tmp_dir, "/cache/", y_val, "_bootstrap_pearsons_CACHED.rds"))
}


# -------------------------------------------------------------------------
# Bootstrap estimates of community biomass vs. realized Shannon (for y-axis)

if (READ_CACHE)
{
    assign("estimates",          readRDS(file = paste0(tmp_dir, "/cache/", y_val, "_brms_models_CACHED.rds")),        envir = .GlobalEnv)
    assign("bootstrap.pearsons", readRDS(file = paste0(tmp_dir, "/cache/", y_val, "_bootstrap_pearsons_CACHED.rds")), envir = .GlobalEnv)
} else if (!exists("estimates")) {
    assign("estimates",          readRDS(file = paste0(model_data_dir, "/", y_val, "_brms_models.rds")), envir = .GlobalEnv)
}

# across effects
across.effects <- estimates %>%
    select(Model, x_val, acrossEffect)

bootstrap.across.estimate <- across.effects %>%
    mutate(bootstrap.estimate = map2(.x = acrossEffect,
                                     .y = x_val,
                                     .f = ~
                                         {
                                             id = seq(1, bootstrap_N)

                                             estimates <- posterior_samples(.x, paste0(":", .y)) %>%
                                                 sample_n(bootstrap_N,
                                                          replace = TRUE)

                                             estimates <- cbind(estimates, id)

                                             estimates <- estimates %>%
                                                 pivot_longer(cols = -id,
                                                              names_to = "parameter",
                                                              values_to = "estimate") %>%
                                                 mutate(Stage = ifelse(str_detect(parameter, "Withseedinflow"),
                                                                       "With seed inflow",
                                                                       "Without seed inflow"),
                                                        Stage = as.factor(Stage))

                                             estimates <- estimates %>%
                                                 select(Stage, id, estimate)

                                             return(estimates)
                                         })) %>%
    unnest(bootstrap.estimate)

bootstrap.across.effect <- inner_join(bootstrap.pearsons,
                                      bootstrap.across.estimate) %>%
    ungroup() %>%
    select(x_val, Model, Stage, id, r, estimate)

# within effects
within.effects <- estimates %>%
    select(x_val, Model, withinEffect)

bootstrap.within.estimates <- within.effects %>%
    mutate(bootstrap.estimates = map2(.x = withinEffect,
                                      .y = x_val,
                                      .f = ~
                                          {
                                              id = seq(1, bootstrap_N)

                                              estimates <- posterior_samples(.x, paste0(":", .y)) %>%
                                                  sample_n(bootstrap_N,
                                                           replace = TRUE)

                                              estimates <- cbind(estimates, id)

                                              estimates <- estimates %>%
                                                  pivot_longer(cols = -id,
                                                               names_to = "parameter",
                                                               values_to = "estimate") %>%
                                                  mutate(Stage = ifelse(str_detect(parameter, "Withseedinflow"),
                                                                        "With seed inflow",
                                                                        "Without seed inflow"),
                                                         Stage = as.factor(Stage)) %>%
                                                  mutate(Ninitial = str_extract(parameter, "\\d+"),
                                                         Ninitial = factor(Ninitial, levels = c("2", "4", "8", "16", "32")))

                                              estimates <- estimates %>%
                                                  select(Stage, Ninitial, id, estimate)

                                              return(estimates)
                                          })) %>%
    unnest(bootstrap.estimates)

bootstrap.within.effect <- inner_join(bootstrap.pearsons, bootstrap.within.estimates) %>%
    ungroup() %>%
    select(x_val, Model, Stage, Ninitial, id, r, estimate)

if (SAVE_CACHE)
{
    saveRDS(bootstrap.across.effect, file = paste0(tmp_dir, "/cache/", y_val, "_across_effect_CACHED.rds"))
    saveRDS(bootstrap.within.effect, file = paste0(tmp_dir, "/cache/", y_val, "_within_effect_CACHED.rds"))
}


# -------------------------------------------------------------------------
# Run bootstrapped Type II Regression (reduced major axis regression)

if (READ_CACHE)
{
    assign("bootstrap.across.effect", readRDS(file = paste0(tmp_dir, "/cache/", y_val, "_across_effect_CACHED.rds")), envir = .GlobalEnv)
    assign("bootstrap.within.effect", readRDS(file = paste0(tmp_dir, "/cache/", y_val, "_within_effect_CACHED.rds")), envir = .GlobalEnv)
}

# across effect
bootstrap.across.effect <- bootstrap.across.effect %>%
    group_by(x_val, Stage, id) %>%
    nest() %>%
    mutate(model = map(.x = data,
                       .f = ~ lmodel2::lmodel2(.x$estimate ~ .x$r)))

regression.across.effect <- bootstrap.across.effect %>%
    mutate(regression.results = map(.x = model,
                                    .f = ~ .x$regression.results[3,])) %>% # SMA Regression
    unnest(regression.results)

# within effect
bootstrap.within.effect <- bootstrap.within.effect %>%
    group_by(x_val, Stage, Ninitial, id) %>%
    nest() %>%
    mutate(model = map(.x = data,
                       .f = ~ lmodel2::lmodel2(.x$estimate ~ .x$r)))

regression.within.effect <- bootstrap.within.effect %>%
    mutate(regression.results = map(.x = model,
                                    .f = ~ .x$regression.results[3,])) %>% # SMA Regression
    unnest(regression.results)

if (SAVE_CACHE)
{
    saveRDS(regression.across.effect, file = paste0(tmp_dir, "/cache/", y_val, "_regression_across_effect_CACHED.rds"))
    saveRDS(regression.within.effect, file = paste0(tmp_dir, "/cache/", y_val, "_regression_within_effect_CACHED.rds"))
}


# -------------------------------------------------------------------------
# Generate plots

if (READ_CACHE)
{
    assign("estimates",                readRDS(file = paste0(tmp_dir, "/cache/", y_val, "_brms_models_CACHED.rds")),              envir = .GlobalEnv)
    assign("regression.across.effect", readRDS(file = paste0(tmp_dir, "/cache/", y_val, "_regression_across_effect_CACHED.rds")), envir = .GlobalEnv)
    assign("regression.within.effect", readRDS(file = paste0(tmp_dir, "/cache/", y_val, "_regression_within_effect_CACHED.rds")), envir = .GlobalEnv)
}

# across effect
across.points <- regression.across.effect %>%
    ungroup() %>%
    select(x_val, Stage, data) %>%
    unnest(data)

across.points.prediction.interval <- regression.across.effect %>%
    group_by(x_val, Stage) %>%
    select(x_val, Stage, Intercept, Slope) %>%
    nest() %>%
    mutate(predictions = map(.x = data,
                             .f = ~ {
                                 lines <- data.frame(Slope = .x$Slope, Intercept = .x$Intercept)

                                 predictions <- expand_grid(lines, X = seq(-0.4, 1, length.out = 100)) %>%
                                     mutate(Y = Slope * X + Intercept)
                             })) %>%
    select(-data) %>%
    unnest(predictions)

across.points.prediction.interval <- across.points.prediction.interval %>%
    group_by(x_val, Stage, X) %>%
    select(x_val, Stage, X, Y) %>%
    nest() %>%
    mutate(quantiles = map(.x = data,
                           .f = ~ {
                               data.frame(ci_upper = quantile(.x$Y, 0.975),
                                          ci_mid   = quantile(.x$Y, 0.50),
                                          ci_lower = quantile(.x$Y, 0.025))
                           })) %>%
    select(-data) %>%
    unnest(quantiles)

across.effect.slope.CI <- regression.across.effect %>%
    group_by(x_val, Stage) %>%
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

walk(.x = unique(estimates$x_val),
     .f = ~
         {
             across.points_ <- across.points %>%
                 filter(x_val == .x)

             across.points.prediction.interval_ <- across.points.prediction.interval %>%
                 filter(x_val == .x)

             across.effect.slope.CI_ <- across.effect.slope.CI %>%
                 filter(x_val == .x) %>%
                 mutate(lower_ci = round(lower_ci, digits = 2),
                        upper_ci = round(upper_ci, digits = 2))

             x_ <- across.points.prediction.interval_ %>%
                 ungroup() %>%
                 filter(X == max(X)) %>%
                 distinct(X) %>%
                 pull(X)

             y_ <- across.points.prediction.interval_ %>%
                 ungroup() %>%
                 filter(ci_lower == min(ci_lower)) %>%
                 distinct(ci_lower) %>%
                 pull(ci_lower)

             across.effect.slope.CI_ <- across.effect.slope.CI_ %>%
                 mutate(x_max = x_,
                        y_min = y_)

             across.points.corr_ <- across.points_ %>%
                 group_by(Stage) %>%
                 summarise(corr = round(cor(r, estimate), digits = 2)) %>%
                 mutate(x_max = x_,
                        y_min = y_)

             a <- ggplot() +
                 geom_hline(yintercept = 0, color = "darkgray") +
                 geom_vline(xintercept = 0, color = "darkgray") +
                 geom_point(data = across.points_,
                            mapping = aes(x = r,
                                          y = estimate,
                                          color = Model),
                            alpha = 0.10) +
                 scale_color_brewer(palette = "Dark2") +
                 geom_ribbon(data = across.points.prediction.interval_,
                             mapping = aes(x = X,
                                           ymin = ci_lower,
                                           ymax = ci_upper),
                             fill = "darkgrey",
                             alpha = 0.66) +
                 scale_fill_brewer(palette = "Greys") +
                 geom_line(data = across.points.prediction.interval_,
                           mapping = aes(x = X,
                                         y = ci_mid),
                           color = "#3E5EF1") +
                 geom_text(data = across.effect.slope.CI_,
                           mapping = aes(x = x_max,
                                         y = y_min,
                                         label = paste0("95% CI of slope: [", lower_ci, ", ", upper_ci, "]")),
                           hjust = 1) +
                 geom_text(data = across.points.corr_,
                           mapping = aes(x = x_max,
                                         y = y_min,
                                         label = paste0("Pearson's correlation: ", corr)),
                           hjust = 1,
                           nudge_y = ifelse(.x == "Shannon", 5, 1.15)) +
                 facet_grid(cols = vars(Stage)) +
                 labs(x = "Function-dominance correlation",
                      y = "Slope of BEF relationship") +
                 guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
                 theme_few(15) +
                 theme(aspect.ratio = 0.618)

             cowplot::save_plot(a,
                                filename = paste0(tmp_dir, "/", y_val, "_acrossEffects_", .x, ".png"),
                                ncol = 2,
                                nrow = 1,
                                base_asp = 1.75)

         }
)

# within effect
within.points <- regression.within.effect %>%
    ungroup() %>%
    select(x_val, Stage, Ninitial, data) %>%
    unnest(data)

within.points.prediction.interval <- regression.within.effect %>%
    group_by(x_val, Ninitial, Stage) %>%
    select(x_val, Ninitial, Stage, Intercept, Slope) %>%
    nest() %>%
    mutate(predictions = map(.x = data,
                             .f = ~
                                 {
                                     lines <- data.frame(Slope = .x$Slope, Intercept = .x$Intercept)

                                     predictions <- expand_grid(lines, X = seq(-0.4, 1, length.out = 100)) %>%
                                         mutate(Y = Slope * X + Intercept)
                                 })) %>%
    select(-data) %>%
    unnest(predictions)

within.points.prediction.interval <- within.points.prediction.interval %>%
    select(x_val, Ninitial, Stage, X, Y) %>%
    group_by(x_val, Ninitial, Stage, X) %>%
    nest() %>%
    mutate(quantiles = map(.x = data,
                           .f = ~
                               {
                                   data.frame(ci_upper = quantile(.x$Y, 0.975),
                                              ci_mid   = quantile(.x$Y, 0.50),
                                              ci_lower = quantile(.x$Y, 0.025))
                               })) %>%
    select(-data) %>%
    unnest(quantiles)

within.effect.slope.CI <- regression.within.effect %>%
    group_by(x_val, Ninitial, Stage) %>%
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

x_vals <- unique(estimates$x_val)
x_vals <- x_vals[x_vals != "Richness"] # Richness has no variation, so the plot is meaningless
walk(.x = x_vals,
     .f = ~
         {
             within.points_ <- within.points %>%
                 filter(x_val == .x)

             within.points.prediction.interval_ <- within.points.prediction.interval %>%
                 filter(x_val == .x)

             within.effect.slope.CI_ <- within.effect.slope.CI %>%
                 filter(x_val == .x) %>%
                 mutate(lower_ci = round(lower_ci, digits = 2),
                        upper_ci = round(upper_ci, digits = 2))

             x_ <- within.points.prediction.interval_ %>%
                 ungroup() %>%
                 filter(X == min(X)) %>%
                 distinct(X) %>%
                 pull(X)

             y_ <- within.points_ %>%
                 ungroup() %>%
                 filter(estimate == min(estimate)) %>%
                 distinct(estimate) %>%
                 pull(estimate)

             within.effect.slope.CI_ <- within.effect.slope.CI_ %>%
                 mutate(x_min = x_,
                        y_min = y_)

             within.points.corr_ <- within.points_ %>%
                 group_by(Stage, Ninitial) %>%
                 summarise(corr = round(cor(r, estimate), digits = 2)) %>%
                 mutate(x_min = x_,
                        y_min = y_)

             p <- ggplot() +
                 geom_hline(yintercept = 0, color = "darkgray") +
                 geom_vline(xintercept = 0, color = "darkgray") +
                 geom_point(data = within.points_,
                            mapping = aes(x = r,
                                          y = estimate,
                                          color = Model),
                            alpha = 0.10) +
                 geom_ribbon(data = within.points.prediction.interval_,
                             mapping = aes(x = X,
                                           ymin = ci_lower,
                                           ymax = ci_upper),
                             fill = "darkgrey",
                             alpha = 0.66) +
                 geom_line(data = within.points.prediction.interval_,
                           mapping = aes(x = X,
                                         y = ci_mid),
                           color = "#9f140f",
                           size = 1) +
                 geom_text(data = within.effect.slope.CI_,
                           mapping = aes(x = x_min,
                                         y = y_min,
                                         label = paste0("95% CI of slope: [", lower_ci, ", ", upper_ci, "]")),
                           hjust = 0) +
                 geom_text(data = within.points.corr_,
                           mapping = aes(x = x_min,
                                         y = y_min,
                                         label = paste0("Pearson's correlation: ", corr)),
                           hjust = 0,
                           nudge_y = 17.5) +
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
                                filename = paste0(tmp_dir, "/", y_val, "_withinEffects_", .x, ".png"),
                                ncol = 2,
                                nrow = 5,
                                base_asp = 1.75)
         }
)
