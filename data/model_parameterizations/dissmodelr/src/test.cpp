// pkgbuild::compile_dll()
//   devtools::document()

#include <Rcpp.h>
#include "UPlant.h"
#include "USimu.h"
using namespace Rcpp;

//' Runs dissmodel
//' @name simulate
//' @param expParFilename Name of input file
//' @return Time series of biomass and productivity of all species
//
// [[Rcpp::export]]
List simulate(std::string expParFilename) {
  Simulation testSimu;
  vector<ExperimentParameters> expPars;
  unsigned int rc;
  ExperimentParameters eP;
  eP.readExperimentParametersNew(expPars, expParFilename);
  // NumericVector result(expPars.size());

  for (size_t i = 0; i < expPars.size(); i++)
  {
    eP = expPars[i];
    // result[i] = eP.initSoilWaterContent;
    testSimu.init(eP);
    //for (size_t j = 0; j < eP.repetitions; j++)
//    size_t j = 0;
//      testSimu.initOutBiomass(eP.outputFilename, j);
//      testSimu.speciesGrowth = 0.0 * testSimu.speciesGrowth;
//      testSimu.speciesBiomass = 0.0 * testSimu.speciesBiomass;

      rc = eP.samplingDate;
      testSimu.run();
//      testSimu.gOutBiomass.close();
  }
  List result = List::create(
    Named( "biomass" ) = testSimu.speciesBiomass,
    Named( "productivity" ) = testSimu.speciesGrowth,
    Named( "seedproduction") = testSimu.speciesSeedproduction
  );
  return result;
}

// library(devtools)
// install()
// setwd("/Users/bjoern/gitlab/dissmodelr")
// setwd("/Users/bjoern/Dropbox/Projects/Simnet/dissmodel68i")
// dissmodelr:::dissmodel_base(2)
