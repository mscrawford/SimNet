#' Create succulent model configuration text file
#'
#' @param repetitions Number of repetitions to do. UNUSED
#' @param gridWidth Number of grid columns
#' @param gridHeight Number of grid rows
#' @param cellwidth Length of cells (cm)
#' @param numyears Number of years to simulate
#' @param timeStepsPerYear Number of time steps per year to simulate
#' @param timestepLength Duration of individual time step (days)
#' @param soilDepth Depth of soil layer (m)
#' @param thetaSat Soil water saturation level (fraction)
#' @param kappa TODO
#' @param mu TODO
#' @param sd TODO
#' @param seedrainIntensity Number of external seeds (sum across all species) (m-2)
#' @param seedrainYears Number of years during with external seedrain (starts in first year)
#' @param seedMortality TODO
#' @param maxIndividualDensity TODO
#' @param sensitivityAnalysis TODO
#' @param initSoilWaterContent TODO
#' @param redistributionFactorFile TODO
#' @param dailyRainFile TODO
#' @param plantInitFile TODO
#' @param plantParsFile TODO
#' @param metacommunityFile TODO
#' @param outputFilename TODO
#' @param outputOptions TODO
#' @param reportCycle TODO
#' @param dailyEvapoFile TODO
#' @param simulationType TODO
#' @param initSeedbankDensity TODO
#' @param rocksAndPotsFile TODO
#' @param thetaCrit TODO
#' @param dewAmount TODO
#' @param dewDays TODO
#' @return Text with configuration file
#' @export
make_configuration <- function(repetitions = 1, gridWidth = 20, gridHeight = 20,
                               cellwidth = 0.25, numyears = 500, timeStepsPerYear = 365,
                               timestepLength = 1, soilDepth = 0.3, thetaSat = 0.4,
                               kappa = 0.7, mu = -1, sd = 0.8, seedrainIntensity = 30,
                               seedrainYears = 100,
                               seedMortality = 0.3, maxIndividualDensity = 499.9,
                               sensitivityAnalysis = "false", initSoilWaterContent = 0.2,
                               redistributionFactorFile = system.file("input/diversityMechanisms5/heterogeneousLandscapeSmall.txt", package = "dissmodelr"),
                               dailyRainFile = system.file("input/diversityMechanisms5/rain90mm5000.dat", package = "dissmodelr"),
                               plantInitFile = "none",
                               plantParsFile = system.file("input/diversityMechanisms5/PlantPars.dat", package = "dissmodelr"),
                               metacommunityFile = "none",
                               outputFilename = "outSim17",
                               outputOptions = "popsizeseedscoverbiomassswcdeathwateruptake",
                               reportCycle = 1825,
                               dailyEvapoFile = system.file("input/diversityMechanisms5/numeesSouthEvapo.dat", package = "dissmodelr"),
                               simulationType = "complex",
                               initSeedbankDensity = 3000,
                               rocksAndPotsFile = system.file("input/diversityMechanisms5/pots20small.txt", package = "dissmodelr"),
                               thetaCrit = 0.04, dewAmount = 0.1, dewDays = 200) {
  match.call()
  current_arg_values <- allargs()
  paste(paste(names(current_arg_values), unlist(current_arg_values), sep = " = "), collapse = " ")
}
