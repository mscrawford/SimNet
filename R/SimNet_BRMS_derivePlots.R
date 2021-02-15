library(rstan)
library(brms)
library(tidybayes)
library(modelr)
library(sjstats)
library(sjPlot)

library(data.table)
library(tidyverse)
library(cowplot)
library(ggthemes)


# -------------------------------------------------------------------------
# BRMS plotting function

brm.draw.fits.fun <- function(d_, acrossEffect_, withinEffect_, measure_)
{
    .extract_significances_acrossEffect <- function(model_, measure_)
    {
        hdi_ <- bayestestR::hdi(model_) %>%
            mutate(significant = !(0 >= CI_low & 0 <= CI_high))

        hdi_ <- hdi_ %>%
            filter(str_detect(Parameter, measure_)) %>%
            mutate(Stage = ifelse(str_detect(Parameter, "Withseedinflow"), "With seed inflow", "Without seed inflow"),
                   Stage = as.factor(Stage))

        hdi_ <- hdi_ %>%
            select(Stage, significant)

        return(hdi_)
    }

    .extract_significances_withinEffect <- function(model_, measure_)
    {
        hdi_ <- bayestestR::hdi(model_) %>%
            mutate(significant = !(0 >= CI_low & 0 <= CI_high))

        hdi_ <- hdi_ %>%
            filter(str_detect(Parameter, measure_)) %>%
            mutate(Stage = ifelse(str_detect(Parameter, "Withseedinflow"), "With seed inflow", "Without seed inflow"),
                   Stage = as.factor(Stage)) %>%
            mutate(Ninitial = str_extract(Parameter, "\\d+"),
                   Ninitial = factor(Ninitial, levels = c("2", "4", "8", "16", "32")))

        hdi_ <- hdi_ %>%
            select(Stage, Ninitial, significant)

        return(hdi_)
    }

    measure <- sym(measure_)

    acrossEffect <- d_ %>%
        group_by(Stage) %>%
        data_grid(!!measure) %>%
        add_fitted_draws(acrossEffect_) %>%
        inner_join(.extract_significances_acrossEffect(acrossEffect_, measure_))

    withinEffect <- d_ %>%
        filter(Ninitial %in% levels(withinEffect_$data$Ninitial)) %>%
        droplevels() %>%
        group_by(Ninitial, Stage) %>%
        data_grid(!!measure) %>%
        add_fitted_draws(withinEffect_) %>%
        inner_join(.extract_significances_withinEffect(withinEffect_, measure_))

    p <- ggplot() +
        geom_point(data = d_,
                   aes(x = !!measure,
                       y = !!sym(Y_VAL),
                       color = Ninitial),
                   alpha = 0.9) +
        stat_lineribbon(data = acrossEffect,
                        aes(x = !!measure,
                            y = .value,
                            linetype = significant),
                        alpha = 1,
                        fill = "lightgray",
                        .width = c(0.95),
                        show.legend = FALSE) +
        stat_lineribbon(data = withinEffect,
                        aes(x = !!measure,
                            y = .value,
                            color = Ninitial,
                            group = Ninitial,
                            linetype = significant),
                        alpha = 0.66,
                        fill = "lightgray",
                        .width = c(0.95),
                        show.legend = FALSE) +
        facet_grid(. ~ Stage) +
        scale_color_viridis_d() +
        scale_linetype_manual(values = c("twodash", "solid"), breaks = c(FALSE, TRUE)) +
        ylim(NA, 110) +
        labs(x = paste0("Realized ", measure),
             y = paste0("Total ", Y_VAL),
             color = "Planted species\nrichness") +
        guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3))) +
        theme_few(18) +
        theme(aspect.ratio = 0.5)

    return(p)
}


# -------------------------------------------------------------------------
# Generate fitted lines

# Read brms model data
if (READ_CACHE)
{
    assign("estimates", readRDS(file = paste0(tmp_dir, "/cache/", Y_VAL, "_brms_models_CACHED.rds")), envir = .GlobalEnv)
} else if (!exists("estimates")) {
    assign("estimates", readRDS(file = paste0(model_data_dir, "/", Y_VAL, "_brms_models.rds")), envir = .GlobalEnv)
}

estimates <- estimates %>%
    mutate(fitted_plot = pmap(.l = list(data, acrossEffect, withinEffect, measure),
                              .f = ~ brm.draw.fits.fun(..1, ..2, ..3, ..4)))

walk(.x = unique(estimates$measure),
     .f = ~
         {
             p.list <- (estimates %>% filter(measure == .x))$fitted_plot
             p.labels <- (estimates %>% filter(measure == .x))$Model

             p <- cowplot::plot_grid(plotlist = map(.x = p.list,
                                                    .f = ~ .x + theme(legend.position = "none")),
                                     labels = p.labels,
                                     ncol = 1,
                                     nrow = length(p.list))

             legend <- cowplot::get_legend(p.list[[1]])

             p.legend <- cowplot::plot_grid(p, legend, rel_widths = c(2, 0.3))

             cowplot::save_plot(p.legend,
                                filename = paste0(tmp_dir, "/", Y_VAL, "_fitted_", .x, ".png"),
                                nrow = length(p.list),
                                ncol = 2,
                                base_asp = 2)
         }
)
