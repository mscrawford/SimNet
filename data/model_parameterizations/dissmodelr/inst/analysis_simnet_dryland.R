
result_path <- "/Volumes/ExtremeSSD/simnet/results"

files <- list.files(result_path, pattern = "^res(.*)rds$")

filename_to_df <- function(x) {
  x <- basename(x)
  x <- gsub("\\.rds", "", x)
  x <- gsub("res_", "", x)
  x <- strsplit(x, "_")
  x <- sapply(x, as.numeric)
  x <- t(x)
  colnames(x) <- c("Ninitial", "Rep", "SeedRain", "Wdh")
  as.data.frame(x)
}

library(dissmodelr)
df <- make_simnet_simulation_df(Wdh = 1)
set.seed(2) # This is important - it is the same seed that is used in run_simnet
pools <- make_species_pools(df)

get_species_id <- function(Ninitial, Rep, pools) {
  pools[[paste0(formatC(Ninitial, width = 2, flag="0"), "_", formatC(Rep, width = 2, flag="0"))]]
}

library("dplyr")
library("tidyr")

read_as_result_df <- function(file, result_path, pools) {
  res <- readRDS(file.path(result_path, file))

  res$productivity[!is.finite(res$productivity)] <- 0 # We replace NA in productivity with 0

  Ninitial <- filename_to_df(file)$Ninitial
  Rep <- filename_to_df(file)$Rep
  SeedRain <- filename_to_df(file)$SeedRain
  Wdh <- filename_to_df(file)$Wdh

  biomass2 <- apply(res$biomass, 2, stats::filter, filter = rep(1/2, 2), sides = 1)
  biomass11 <- apply(res$biomass, 2, stats::filter, filter = rep(1/11, 11), sides = 1)
  biomass <- biomass2[seq(2, 400, by = 2), , drop = FALSE]
  biomass11 <- biomass11[seq(2, 400, by = 2), , drop = FALSE]
  species_ids <- get_species_id(Ninitial, Rep, pools)
  colnames(biomass) <- paste0("biomass_02_", species_ids)
  colnames(biomass11) <- paste0("biomass_11_", species_ids)

  productivity2 <- apply(res$productivity, 2, stats::filter, filter = rep(1/2, 2), sides = 1)
  productivity11 <- apply(res$productivity, 2, stats::filter, filter = rep(1/11, 11), sides = 1)

  productivity <- productivity2[seq(2, 400, by = 2), , drop = FALSE]
  productivity11 <- productivity11[seq(2, 400, by = 2), , drop = FALSE]
  colnames(productivity) <- paste0("productivity_02_", species_ids)
  colnames(productivity11) <- paste0("productivity_11_", species_ids)


  biomass <- data.frame("Stage" = rep(c("assembly", "disassembly"), each = 100),
                        "Year" = 1:200, biomass, biomass11)
  productivity <- data.frame("Stage" = rep(c("assembly", "disassembly"), each = 100),
                             "Year" = 1:200, productivity, productivity11)


  bml <- pivot_longer(biomass, cols = starts_with("bio"),
               names_to = c("Smooth", "SpeciesID"),
               names_prefix = "biomass_",
               names_sep = "_",
               values_to = "Biomass")
  bml$SpeciesID <- as.integer(bml$SpeciesID)


  prl <- pivot_longer(productivity, cols = starts_with("prod"),
                      names_to = c("Smooth", "SpeciesID"),
                      names_prefix = "productivity_",
                      names_sep = "_",
                      values_to = "Productivity")
  prl$SpeciesID <- as.integer(prl$SpeciesID)
  cool <- left_join(bml, prl, by = c("Stage", "Year", "SpeciesID", "Smooth"))
  cool <- cbind(filename_to_df(file), cool)
  cool$Model <- "bjoern"
  cool$Wdh <- Wdh
  rownames(cool) <- NULL
  cool[, c("Model", "Ninitial", "Rep", "SeedRain", "SpeciesID", "Stage", "Year", "Biomass", "Productivity", "Wdh", "Smooth")]
}

read_and_summarise <- function(files, result_path, pools) {
  subset_res <- lapply(files, read_as_result_df, result_path, pools)
  subset_res <- bind_rows(subset_res)
  subset_grouped <- group_by(subset_res, Model, Ninitial, Rep, SeedRain, SpeciesID, Stage, Year, Smooth)
  summarise(subset_grouped, Biomass = mean(Biomass), Productivity = mean(Productivity))
}

# remove the _00.rds files - they come from an earlier test run
files <- files[-grep("_00\\.rds$", files)]
file_groups <- apply(filename_to_df(files)[, c("Ninitial", "Rep", "SeedRain")], 1, paste, collapse = "_")

full_res_averaged <- tapply(files, file_groups, read_and_summarise, result_path, pools)
full_res_averaged <- bind_rows(full_res_averaged)

ful_res_averaged_unsmoothed <- subset(full_res_averaged, Smooth == "02")
ful_res_averaged_smoothed <- subset(full_res_averaged, Smooth == "11")

ful_res_averaged_smoothed$Biomass[ful_res_averaged_smoothed$Year %in% c(1:5, 101:105)] <- NA
ful_res_averaged_smoothed$Productivity[ful_res_averaged_smoothed$Year %in% c(1:5, 101:105)] <- NA

saveRDS(ful_res_averaged_unsmoothed, file.path(result_path, "bjoern_averaged_NAreplaced0.rds"))
saveRDS(ful_res_averaged_smoothed, file.path(result_path, "bjoern_averaged_smooth_NAreplaced0.rds"))
