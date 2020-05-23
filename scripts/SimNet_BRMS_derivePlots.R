library(rstan)
library(brms)
library(tidybayes)
library(modelr)
library(sjstats)
library(sjPlot)

library(data.table)
library(tidyverse)
library(cowplot)

# ---------------------------------------------------------------------------------------------
# Focal ecosystem function

# Y_VAL = "biomass"
Y_VAL = "productivity"


# ---------------------------------------------------------------------------------------------
# Load BRMS model results

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
d <- readRDS(file = paste0("../data/brms_models/", Y_VAL, "_brms_models.rds"))
setwd("../tmp")


# ---------------------------------------------------------------------------------------------
# BRMS plotting function

extract_significances_mainEffect <- function(model_, measure_)
{
    hdi_ <- bayestestR::hdi(model_) %>%
        mutate(significant = !(0 >= CI_low & 0 <= CI_high))

    hdi_ <- hdi_ %>%
        filter(str_detect(Parameter, measure_)) %>%
        mutate(Stage = ifelse(str_detect(Parameter, "metacommunity"), "metacommunity", "isolation"),
               Stage = factor(Stage, levels = c("metacommunity", "isolation")))

    hdi_ <- hdi_ %>%
        select(Stage, significant)

    return(hdi_)
}

extract_significances_treatmentEffect <- function(model_, measure_)
{
    hdi_ <- bayestestR::hdi(model_) %>%
        mutate(significant = !(0 >= CI_low & 0 <= CI_high))

    hdi_ <- hdi_ %>%
        filter(str_detect(Parameter, measure_)) %>%
        mutate(Stage = ifelse(str_detect(Parameter, "metacommunity"), "metacommunity", "isolation"),
               Stage = factor(Stage, levels = c("metacommunity", "isolation"))) %>%
        mutate(Ninitial = str_extract(Parameter, "\\d+"),
               Ninitial = factor(Ninitial, levels = c("2", "4", "8", "16", "32")))

    hdi_ <- hdi_ %>%
        select(Stage, Ninitial, significant)

    return(hdi_)
}

brm.draw.fits.fun <- function(d_, mainEffect_, treatmentEffect_, measure_)
{
    measure <- sym(measure_)

    mainEffect <- d_ %>%
        group_by(Stage) %>%
        data_grid(!!measure) %>%
        add_fitted_draws(mainEffect_) %>%
        inner_join(extract_significances_mainEffect(mainEffect_, measure_))

    treatmentEffect <- d_ %>%
        filter(Ninitial %in% levels(treatmentEffect_$data$Ninitial)) %>%
        droplevels() %>%
        group_by(Ninitial, Stage) %>%
        data_grid(!!measure) %>%
        add_fitted_draws(treatmentEffect_) %>%
        inner_join(extract_significances_treatmentEffect(treatmentEffect_, measure_))

    p <- ggplot() +
        geom_point(data = d_,
                   aes(x = !!measure,
                       y = !!sym(Y_VAL),
                       color = Ninitial),
                   alpha = 0.9) +
        stat_lineribbon(data = treatmentEffect,
                        aes(x = !!measure,
                            y = .value,
                            color = Ninitial,
                            group = Ninitial,
                            linetype = significant),
                        alpha = 0.5,
                        show.legend = FALSE) +
        stat_lineribbon(data = mainEffect,
                        aes(x = !!measure,
                            y = .value,
                            linetype = significant),
                        alpha = 0.5,
                        show.legend = FALSE) +
        scale_linetype_manual(values = c("twodash", "solid")) +
        facet_grid(. ~ Stage) +
        labs(x = paste0("Realized ", measure_),
             y = paste0("Total ", Y_VAL),
             color = "Planted species\nrichness") +
        scale_color_viridis_d() +
        scale_fill_brewer(palette = "Greys") +
        theme_bw(16) +
        theme(aspect.ratio = 1)

    return(p)
}


# ---------------------------------------------------------------------------------------------
# Main plots

# Generate plots
d <- d %>%
    mutate(fitted_plot = pmap(.l = list(data, mainEffect, treatmentEffect, measure),
                              .f = ~ brm.draw.fits.fun(..1, ..2, ..3, ..4)))

d <- d %>%
    filter(Model != "IBC_grass.noNDD")

d <- d %>% filter(measure == "Shannon")

# Save plots
walk(.x = unique(d$measure),
     .f = ~
         {
             p.list <- (d %>% filter(measure == .x))$fitted_plot

             p <- cowplot::plot_grid(plotlist = lapply(p.list,
                                                       function(x) {
                                                           x + theme(legend.position="none")
                                                       }),
                                     labels = c("Adam", "Lindsay", "IBC-grass", "PPA", "TROLL", "Björn"),
                                     ncol = 2,
                                     nrow = length(p.list) / 2)

             legend <- get_legend(p.list[[1]])

             p.legend <- cowplot::plot_grid(p, legend, rel_widths = c(4, 0.66))

             cowplot::save_plot(p.legend,
                                filename = paste0(Y_VAL, "_fitted_", .x, ".png"),
                                ncol = 4,
                                nrow = length(p.list) / 2,
                                base_asp = 1.1)
         }
)


# ---------------------------------------------------------------------------------------------
# Trait plots

# Load all the model results
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
d <- readRDS(file = "brms_models_FDis_perTrait.rds")
setwd("../../tmp")

# Generate plots
d <- d %>%
    mutate(plot = pmap(.l = list(data, mainEffect, treatmentEffect, measure),
                       .f = ~ brm.draw.fits.fun(..1, ..2, ..3, ..4)))

d <- d %>%
    filter(Model != "IBC_grass.noNDD")

# Save trait plots
pwalk(.l = expand_grid(unique(d$measure),
                       unique(d$Model)),
      .f = ~
          {
              p.list <- d %>% filter(measure == .x, Model == .y)

              p <- cowplot::plot_grid(plotlist = p.list$plot,
                                      ncol = 1,
                                      nrow = nrow(p.list),
                                      labels = p.list$trait)

              cowplot::save_plot(p,
                                 filename = paste(.y, "_", .x, "_singleTraits.png", sep = ""),
                                 ncol = 1,
                                 nrow = nrow(p.list),
                                 limitsize = FALSE)
          }
)
