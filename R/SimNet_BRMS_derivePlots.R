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

brm.draw.fits.fun <- function(d, acrossEffect, withinEffect, x_val)
{
    .extract_significances_acrossEffect <- function(model, x_val)
    {
        hdi_ <- bayestestR::hdi(model) %>%
            mutate(significant = !(0 >= CI_low & 0 <= CI_high))

        hdi_ <- hdi_ %>%
            filter(str_detect(Parameter, x_val)) %>%
            mutate(Stage = ifelse(str_detect(Parameter, "Withseedinflow"), "With seed inflow", "Without seed inflow"),
                   Stage = as.factor(Stage))

        hdi_ <- hdi_ %>%
            select(Stage, significant)

        return(hdi_)
    }

    .extract_significances_withinEffect <- function(model, x_val)
    {
        hdi_ <- bayestestR::hdi(model) %>%
            mutate(significant = !(0 >= CI_low & 0 <= CI_high))

        hdi_ <- hdi_ %>%
            filter(str_detect(Parameter, x_val)) %>%
            mutate(Stage = ifelse(str_detect(Parameter, "Withseedinflow"), "With seed inflow", "Without seed inflow"),
                   Stage = as.factor(Stage)) %>%
            mutate(Ninitial = str_extract(Parameter, "\\d+"),
                   Ninitial = factor(Ninitial, levels = c("2", "4", "8", "16", "32")))

        hdi_ <- hdi_ %>%
            select(Stage, Ninitial, significant)

        return(hdi_)
    }

    acrossEffect <- d %>%
        group_by(Stage) %>%
        data_grid(!!sym(x_val)) %>%
        add_fitted_draws(acrossEffect) %>%
        inner_join(.extract_significances_acrossEffect(acrossEffect, x_val))

    withinEffect <- d %>%
        filter(Ninitial %in% levels(withinEffect$data$Ninitial)) %>%
        droplevels() %>%
        group_by(Ninitial, Stage) %>%
        data_grid(!!sym(x_val)) %>%
        add_fitted_draws(withinEffect) %>%
        inner_join(.extract_significances_withinEffect(withinEffect, x_val))

    p <- ggplot() +
        geom_point(data = d,
                   aes(x = !!sym(x_val),
                       y = !!sym(y_val),
                       color = Ninitial),
                   alpha = 0.9) +
        stat_lineribbon(data = acrossEffect,
                        aes(x = !!sym(x_val),
                            y = .value,
                            linetype = significant),
                        alpha = 1,
                        fill = "lightgray",
                        .width = c(0.95),
                        show.legend = FALSE) +
        stat_lineribbon(data = withinEffect,
                        aes(x = !!sym(x_val),
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
        labs(x = paste0("Realized ", x_val),
             y = paste0("Total ", y_val),
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
    assign("estimates", readRDS(file = paste0(tmp_dir, "/cache/", y_val, "_brms_models_CACHED.rds")), envir = .GlobalEnv)
} else if (!exists("estimates")) {
    assign("estimates", readRDS(file = paste0(model_data_dir, "/", y_val, "_brms_models.rds")), envir = .GlobalEnv)
}

estimates <- estimates %>%
    mutate(fitted_plot = pmap(.l = list(data, acrossEffect, withinEffect, x_val),
                              .f = ~ brm.draw.fits.fun(..1, ..2, ..3, ..4)))

walk(.x = unique(estimates$x_val),
     .f = ~
         {
             p.list <- (estimates %>% filter(x_val == .x))$fitted_plot
             p.labels <- (estimates %>% filter(x_val == .x))$Model

             p <- cowplot::plot_grid(plotlist = map(.x = p.list,
                                                    .f = ~ .x + theme(legend.position = "none")),
                                     labels = p.labels,
                                     ncol = 1,
                                     nrow = length(p.list))

             legend <- cowplot::get_legend(p.list[[1]])

             p.legend <- cowplot::plot_grid(p, legend, rel_widths = c(2, 0.3))

             cowplot::save_plot(p.legend,
                                filename = paste0(tmp_dir, "/", y_val, "_fitted_", .x, ".png"),
                                nrow = length(p.list),
                                ncol = 2,
                                base_asp = 2)
         }
)
