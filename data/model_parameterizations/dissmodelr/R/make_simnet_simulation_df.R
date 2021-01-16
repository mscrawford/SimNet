#' Make simulation data frame for simnet
#'
#' @param Model Name of model
#' @param Ninitial Initial number of species
#' @param Rep Replicate of species pool
#' @param SeedRain Percentage of reference seed rain
#' @param Stage Simulation stage
#' @param Wdh Stochastic replicate for given species pool
#' @return Data frame with simulation parameter values
#' @export
make_simnet_simulation_df <- function(Model = "bjoern",
                                           Ninitial = as.integer(2^(0:6)),
                                           Rep = 1:64,
                                           SeedRain = c(0, 5, 10, 50, 100, 1000),
                                           Stage = c("assembly"),
                                           Wdh = 5) {

  simulations <- expand.grid("Model" = Model,
                             "Ninitial" = Ninitial,
                             "Rep" = Rep,
                             "SeedRain" = SeedRain,
                             "Stage" = Stage,
                             "Wdh" = seq(Wdh), stringsAsFactors = FALSE)
  simulations[-which(simulations$Ninitial == 64 & simulations$Rep > 1),]
}
