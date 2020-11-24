library(tidyverse)
library(cowplot)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

PPA <- read.table("../output/PPA_output.csv", sep = "\t", header = T) %>%
    tbl_df()


# Calculate monoculture average internal seed dispersal ---------------------------------------

# Plot to ensure that its the correct dataset
p <- ggplot(PPA,
            aes(x = Year,
                y = N / 5,
                color = as.factor(SpeciesID))) +
    geom_line() +
    labs(x = "Year",
         y = "Yearly seed dispersal") +
    theme_bw() +
    theme(aspect.ratio = 1,
          legend.background = element_blank(),
          legend.position = "none"); p

PPA %>%
    filter(Year == max(Year)) %>%
    mutate(N = N / 5) %>% # Calculate the reproduction on a per year basis
    summarise(`Yearly average seed dispersal` = mean(N)) # Average across all monocultures

save_plot(p, filename = "PPA_MonocultureReproduction.PDF", ncol = 3, nrow = 1)


# 64 species ----------------------------------------------------------------------------------

p <- ggplot(PPA %>% filter(Rep == 10,
                           Ninitial == 32),
            aes(x = Year / 5,
                y = Productivity,
                color = as.factor(SpeciesID))) +
    facet_grid(Rep ~ Ninitial + SeedRain, labeller = label_both) +
    geom_line() +
    geom_vline(xintercept = 200, linetype = 3) +
    labs(x = "Year") +
    theme_bw() +
    theme(aspect.ratio = 1,
          legend.background = element_blank(),
          legend.position = "none"); p

save_plot(p, filename = "PPA_64_Species.PDF", ncol = 3, nrow = 1)
