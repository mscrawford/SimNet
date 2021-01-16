#' Run simnet experiment
#'
#' @param wdh Number of replicates
#' @param base_path Base path for writing files
#' @param mean_annual_seedrain_monoculture_m2 Number of seeds
#' @param test Number of simulations to test (negative values mean use all)
#' @param mcseed Random number seed for simulation
#' @return Data frame with plant parameter values
#' @export
run_simnet <- function(wdh, base_path = "~/simnet",
                       mean_annual_seedrain_monoculture_m2 = 19402,
                       test = -1,
                       mcseed = 3) {
  base_path <- path.expand(base_path)
  simulations <- make_simnet_simulation_df(Wdh = wdh)
#  return(split(simulations[, c("Ninitial", "Rep", "SeedRain", "Wdh")], factor(1:nrow(simulations)))[1:5])
  set.seed(2)
  pools <- make_species_pools(simulations)
  species_64 <- read_plantpars(system.file("input/plantPars_64.txt", package = "dissmodelr"))
  dir.create(base_path, showWarnings = FALSE)
  input_path <- file.path(base_path, "input")
  result_path <- file.path(base_path, "results")
  dir.create(input_path, showWarnings = FALSE)
  dir.create(result_path, showWarnings = FALSE)
  # Write species parameter files
  for (i in names(pools)) {
    write_plantpars(species_64[pools[[i]], ],
                    file.path(input_path, paste0("plantPars_", i, ".txt")))
  }

  # Now do this in parallel
#  library("future")
#  library("future.apply")

  run_single_simulation <- function(x) {
    j <-     paste0(formatC(x[["Ninitial"]], width = 2, flag = "0"),
                    "_",
                    formatC(x[["Rep"]], width = 2, flag = "0"))
    plant_file <- file.path(input_path, paste0("plantPars_", j, ".txt"))
    seed_multiplier <- x[["SeedRain"]]/100
    file_base <- paste0(formatC(x[["Ninitial"]], width = 2, flag = "0"),
                          "_",
                          formatC(x[["Rep"]], width = 2, flag = "0"),
                          "_",
                          formatC(x[["SeedRain"]], width = 2, flag = "0"),
                          "_",
                          formatC(x[["Wdh"]], width = 2, flag = "0"))
    expParFile <- file.path(base_path, paste0("ExpPars_", file_base, ".txt"))
    writeLines(make_configuration(repetitions = 1,
                                  plantParsFile = plant_file,
                                  numyears = 1200,
                                  seedrainIntensity = round(mean_annual_seedrain_monoculture_m2 * seed_multiplier) + 30,
                                  seedrainYears = 600,
                                  outputFilename = paste0("out_biomass_", file_base, ".txt")),
               expParFile)
    res <- dissmodelr:::simulate(expParFile)
    i <- seq(3, 1200, by = 3)
    # Extract key values
    res <- list("biomass" = res$biomass[i,, drop = FALSE],
                 "productivity" = res$productivity[i,, drop = FALSE])
    # Save as .rds
    saveRDS(res, file.path(result_path, paste0("res_", file_base, ".rds")))
    unlink(expParFile)
    rm(res)
  }
  # plan(multiprocess)
  # future_apply(simulations[1:10, c("Ninitial", "Rep", "SeedRain", "Wdh")],
  #              MARGIN = 1L, FUN = run_single_simulation,
  #              future.seed = 0xBEEF)
  # library("parallel")
  RNGkind("L'Ecuyer-CMRG")
  set.seed(3)
  parallel::mc.reset.stream()
  simulation_list <- split(simulations[, c("Ninitial", "Rep", "SeedRain", "Wdh")], factor(1:nrow(simulations)))
  if (test > 0) simulation_list <- simulation_list[1:test]
  parallel::mclapply(simulation_list, run_single_simulation, mc.cores = parallel::detectCores())
}
