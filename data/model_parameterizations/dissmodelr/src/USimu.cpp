//---------------------------------------------------------------------------

#include "USimu.h"
#include "UStringUtils.h"
#include <random>

//---------------------------------------------------------------------------

struct plantInitParams
{
     int speciesID;
     int x;
     int y;
     int birthday;
     double mass;
};

// wrapper around R's RNG such that we get a uniform distribution over
// [0,n) as required by the STL algorithm
// inline int randWrapper(const int n) { return floor(R::unif_rand()*n); };

// https://stackoverflow.com/questions/28033901/stdshuffle-doesnt-compile-with-stdlist
void shuffle(std::list<plant>& l)
{
  std::vector<std::list<plant>::const_iterator> v;
  for (auto it = l.cbegin(); it != l.cend(); ++it) v.push_back(it);
std::shuffle(v.begin(), v.end(), std::mt19937{ std::random_device{}() }); //, randWrapper);
  // std::shuffle(v.begin(), v.end(), std::mt19937{ floor(R::unif_rand() * 10)}); //, randWrapper);
//  std::shuffle(v.begin(), v.end(), std::mt19937{ 5 }); //, randWrapper);
  //  std::shuffle(v.begin(), v.end(), R::unif);
  std::list<plant> shuffled;
  for (auto &it : v) shuffled.splice(shuffled.end(), l, it);
  l.swap(shuffled);
}

 //gsl_rng * theRNG;
 vector<NNcell> mybarriers;

  void initPlants(list<plant>& allPlants, const vector<plantInitParams>& plantInitPars, const vector<plantParams>& plantPars, valarray<bool>& cellOccupied, int width, int height);

//----------------------------------------------------------------------------

  void Simulation::init(const ExperimentParameters & eP)
  {

       experimentPars = eP;
       repetition = 0;

       gridWidth = eP.gridWidth;
       gridHeight = eP.gridHeight;
       cellwidth = eP.cellwidth;
       cellArea = cellwidth * cellwidth;

       timeStepsPerYear = eP.timeStepsPerYear;
       numyears = eP.numyears;
       timestepLength = eP.timestepLength; // 365.0 / double(timeStepsPerYear);
       simulationDuration = numyears * timeStepsPerYear * timestepLength;

       soilDepth = eP.soilDepth;
       thetaSat = eP.thetaSat;
       thetaCrit = eP.thetaCrit;
       uptakeSlope = eP.uptakeSlope;
//       psiE = eP.psiE;
//       b = eP.b;
       initSoilWaterContent =  eP.initSoilWaterContent;

//       airTemperature = eP.airTemperature;
//       airHumidity = eP.airHumidity;
//       evapotranspiration = eP.evapotranspiration;

       initEvapotranspiration(eP);
       kappa = eP.kappa;

       mu = eP.mu;
       sd = eP.sd;
       seedMortality = eP.seedMortality;
       seedrainIntensity = eP.seedrainIntensity;
       seedrainYears = eP.seedrainYears;
       initSeedbankDensity = eP.initSeedbankDensity;

       maxIndividualDensity = eP.maxIndividualDensity;

       minLambda = 0.00001;
       waterDensity = 1.0;

       simulationType = eP.simulationType;
       startGrazingPeriod = eP.startGrazingPeriod;
       endGrazingPeriod = eP.endGrazingPeriod;
       grazingFrequency = eP.grazingFrequency;
       grazingAnnualAmount = eP.grazingAnnualAmount;

       dewAmount = eP.dewAmount;
       dewDays = eP.dewDays;



       readDailyRain(eP.dailyRainFile);
       readPlantInitFile(eP.plantInitFile);
       readPlantPars(eP.plantParsFile);

       offset = -ceil(simulationDuration);

       startEstablishmentReporting = eP.startEstablishmentReporting;
       stopEstablishmentReporting = eP.stopEstablishmentReporting;


//       double psiCrit = plantPars[0].psiCrit;
//       thetaCrit = thetaSat * pow(psiCrit / psiE, -1.0/b);

       metacommunity.resize(plantPars.size(), 1.0);       // ist zu diesem Zeitpunkt plantPars initialisiert? ja, das passiert in readPlantPars
       waterContent.resize(gridWidth*gridHeight, initSoilWaterContent);
       redistributionFactor.resize(gridWidth * gridHeight, 1.0);


       readValarray(redistributionFactor, eP.redistributionFactorFile, "const");       // TODO : How to do the casting?

//     barriers = vector<NNcell>(gridWidth * gridHeight);


     cohortCounter = 0;   // wof?r brauche ich den cohortCounter?     in initRepetition wird er auf 0 gesetzt.
     plantWaterDemand.resize(gridWidth*gridHeight, 0.0);
     maxWaterUptakeRate.resize(gridWidth*gridHeight, 0.0);
     realizedPlantWaterDemand.resize(gridWidth*gridHeight, 0.0);
     leafAreaIndex.resize(gridWidth*gridHeight, 0.0);
     potentialEvaporation.resize(gridWidth*gridHeight, 0.0);
     actualEvaporation.resize(gridWidth*gridHeight, 0.0);
     newSeeds.resize(gridWidth*gridHeight*plantPars.size(), 0); //plantPars.size());
     seedbank.resize(gridWidth*gridHeight*plantPars.size(), 0); //plantPars.size());
//     seedbankNew.resize(gridWidth*gridHeight*plantPars.size(), 0.0);
     plantDensity.resize(gridWidth*gridHeight*plantPars.size(), 0);
      totalDeaths.resize(plantPars.size(), 0.0);
     droughtDeaths.resize(plantPars.size(), 0.0);

     newSeedsOfPastYears.resize(5, newSeeds);


//     dewfallSequence.resize(numyears * timeStepsPerYear, 0.0);

     initDispersalKernel(mu, sd);           // TODO : das k?nnte man davon abh?ngig machen, ob ?berhaupt dispersal stattfindet

     globalDispersal = eP.globalDispersal;

     rocksAndPots.resize(gridWidth * gridHeight); //=
     readVector(rocksAndPots, eP.rocksAndPotsFile, "const");
     initRootPositions(plantPars,  gridHeight, gridWidth, rocksAndPots);
     vector<int> leafBarriers (gridWidth * gridHeight, 0);
     initLeafPositions(plantPars,  gridHeight, gridWidth, leafBarriers);

     adjustRedistributionMap();
     initGerminationDates(plantPars);
     initDispersalDates(plantPars);

     /*
    ofstream tstout;
    tstout.open("rootpositions.txt");
     for (size_t i = 0; i < rootPositions.size(); ++i) {
         for (size_t j = 0; j < rootPositions[i].size(); ++j) {
            tstout << rootPositions[i][j] << ", ";
         }
         tstout << "\n";
     }
     tstout.close();


    tstout.open("germinationdates.txt");
     for (size_t i = 0; i < germinationDates.size(); ++i) {
            tstout << germinationDates[i] << ", ";
     }
     tstout.close();
*/
     initLeafFraction(plantPars);

     if ( simulationType == "simple" ) {
         growth.assign(numyears*timeStepsPerYear, 0.0);        // TODO : was passiert hier? Warum muss assign verwendet werden?
         survivalProbability.assign(numyears*timeStepsPerYear, 0.0); // TODO : was passiert hier? Warum muss assign verwendet werden?
     }
     Rcpp::NumericMatrix g(numyears, plantPars.size());
     speciesGrowth = g;
     Rcpp::NumericMatrix b(numyears, plantPars.size());
     speciesBiomass = b;
     Rcpp::IntegerMatrix pp(numyears, plantPars.size());
     speciesSeedproduction = pp;
     //speciesGrowth = Rcpp::NumericMatrix::create();
  }

  //----------------------------------------------------------------------------

  void Simulation::run()
  {
     initRepetition();
     simulatedTime = 0.0;
     for (year = 0; year < numyears; year++)
     {
       R_CheckUserInterrupt();
       reportBiomass(this);
        for (timestep = 0; timestep < timeStepsPerYear; timestep++)
        {
           // shuffle(allPlants);
           checkObserver(timestep + timeStepsPerYear * year, this);
           simulatedTime += timestepLength;
           rainfall();
           dewfall();
           drainage();
           evaporationAndPlantWaterUptakeNew();

           itAllPlants = allPlants.begin();
           while (itAllPlants != allPlants.end())
           {
              (*itAllPlants).production3(plantPars, timestepLength, simulatedTime);
              speciesGrowth(year, (*itAllPlants).getspeciesid()) += (*itAllPlants).getCohortSize() * (*itAllPlants).getGrowth();

                   if ( simulationType != "simple") {
                       (*itAllPlants).allocationNew(plantPars, cellArea, leafFraction, rootPositions, year, simulatedTime, cellwidth, gridWidth, gridHeight, newSeedsOfPastYears, dispersalDates, (year >= startEstablishmentReporting && year <= stopEstablishmentReporting));
                   }

              if((*itAllPlants).mortality(plantPars, plantDensity, simulatedTime, timestepLength, simulationType, totalDeaths, droughtDeaths) == true)
              {
                //speciesGrowth(year, (*itAllPlants).getspeciesid()) -= (*itAllPlants).getCohortSize() * (*itAllPlants).getGrowth(); // remove growth of dead individuals from summary
                 itAllPlants = allPlants.erase(itAllPlants);
              }
              else
              {
                  itAllPlants++; //if the plant was removed, the itAllPlants already points to the next element in the list
              }
           }
           // reportGrowth(this);

           if ( simulationType != "simple" ) {
                if ( int(simulatedTime) == dispersalDates[year]) {
                   dispersal();
                   if (year < seedrainYears) {
                       externalSeedrain();
                   }
                }

                if (int(simulatedTime) == germinationDates[year]) {
                    germination();
                }
           }
        }

        if ( simulationType != "simple" ) {
             seedbankmortality();
        }

        // In case there are no plants, execution of this simulation experiment is stopped
        if ( allPlants.size() + seedbank.sum() == 0 && seedrainIntensity < 0.0000000001 ) {
            break;
        }
     }
//     closeOutputFiles();
  };
  //----------------------------------------------------------------------------

  void Simulation::initRepetition()  // ?berfl?ssig machen durch Einbau von init in run
  {
      cohortCounter = 0;

      unsigned int seedsPerCell = floor(initSeedbankDensity * cellArea);
      seedbank = seedsPerCell;
      newSeeds = 0;
      plantWaterDemand = 0.0;
      leafAreaIndex = 0.0;
      waterContent = initSoilWaterContent;
      totalDeaths = 0.0;
      droughtDeaths = 0.0;
      cumulativePlantWaterUptake = 0.0;
      cumulativeDrainage = 0.0;
      cumulativeEvaporation = 0.0;


      initPlants();

      if (repetition == 0) offset = experimentPars.initOffset;
      else offset += experimentPars.deltaOffset;

//      offset += ceil(simulationDuration);  /* TODO : initOffset als Parameter, deltaOffset als Parameter */
      /* TODO : Opening of output data files - when are they closed?*/
//      openOutputFiles();
      repetition++;
   };


  //----------------------------------------------------------------------------

  void Simulation::initPlants()  // ebenfalls in run aufrufen
  {
//   plantDensity = valarray<unsigned int>(gridWidth*gridHeight*plantPars.size());        //transferred to init()
   plantDensity = 0;

   allPlants.clear();
   plant newPlant;
   plantInitParams2 lpip;
   for (size_t i = 0; i < plantInitPars.size(); i++)
   {
      lpip = plantInitPars[i];
      unsigned int cellID = gridWidth*lpip.row + lpip.col;
      if (rocksAndPots[cellID] < 2) {
      newPlant.setspeciesid(lpip.speciesID);
      newPlant.setCellID(cellID);
      newPlant.setCohortSize(lpip.cohortSize);
      newPlant.setmass(lpip.mass);
      newPlant.setGrowth(0.0);
      newPlant.setLastTimePositiveGrowth(0.0);
      newPlant.setDesiccating(false);
      newPlant.setBirthday(lpip.birthday);
      newPlant.setRecruitmentYear(-99);
      newPlant.setTimeToMaturity(-1.0);
//      newPlant.initPositionVector(plantPars, gridHeight, gridWidth, cellArea, true, barriers);
//      newPlant.initPositionVector(plantPars, gridHeight, gridWidth, cellArea, false, barriers);
//      newPlant.allocation2(plantPars, cellArea);
//      newPlant.initAllocation(plantPars, cellArea);
      newPlant.initAllocationNew(plantPars, cellArea, rootPositions);
      newPlant.setPlantWaterContent(0.0);
      newPlant.setWaterUptake(0.0);
      newPlant.setTrueWaterUptake(0.0);
      newPlant.setAmin(plantPars);
//      newPlant.setThetaCrit();  //FIXME
      newPlant.setSeedmass(0.0);
      newPlant.setCohortID(cohortCounter);
      allPlants.push_back(newPlant);
      // change the plantDensityValarray
      plantDensity[cellID] += lpip.cohortSize;

      cohortCounter++;
      }
   }
 };

  //----------------------------------------------------------------------------

  void Simulation::initRNG(unsigned int rngSeed)
  {
     // gsl_rng_env_setup() ;
     // gsl_rng_default_seed = rngSeed ;
     // theRNG = gsl_rng_alloc(gsl_rng_default) ;
  };

  //----------------------------------------------------------------------------

    void Simulation::initDispersalKernel(double mu, double sd)
    {
       double centerX, centerY;

       // again: check that the x - y / cell - row conversion is consisten
       if (fmod(gridWidth, 2.0) > 0.5 ) centerX = gridWidth / 2.0;
       else centerX = (gridWidth - 1.0) / 2.0;
       if (fmod(gridHeight, 2.0) > 0.5 ) centerY = gridHeight / 2.0;
       else centerY = (gridHeight - 1.0) / 2.0;

       double zeta = mu; // log(mu / cellwidth);
       double sigma = sd; // log(sd / cellwidth);

       unsigned int nsubcells = 4;
       double dz = 1.0/double(nsubcells);


//       dispersalKernel = valarray<double>(0.0, gridWidth*gridHeight);
       dispersalKernel.resize(gridWidth * gridHeight, 0.0);

       for (size_t i = 0; i < gridWidth; i++)
       {
          for (size_t j = 0; j < gridHeight; j++)
          {
            double density = 0.0;
            for (size_t m = 0; m < nsubcells; m++)
            {
               for (size_t n = 0; n < nsubcells; n++)
               {
                  double x = i + dz * m + dz/2.0;
                  double y = j + dz * n + dz/2.0;
                  double dist = distance(x, y, centerX, centerY, gridWidth, gridHeight);
                  // check this 2D factor !
                  double circumferencefactor;
                  if (dist > 0.0) circumferencefactor = 1.0/(2*M_PI*dist);
                  else circumferencefactor = 1.0;
//                  density += circumferencefactor * gsl_ran_lognormal_pdf(dist, zeta, sigma);
//                  density += circumferencefactor * gsl_ran_lognormal_pdf(dist * cellwidth, zeta, sigma);  // dist is measured in number of cells, whereas zeta and sigma are log(m), such that the function expects dist to be scaled in m.
                  density += circumferencefactor * R::dlnorm(dist * cellwidth, zeta, sigma, false);  // dist is measured in number of cells, whereas zeta and sigma are log(m), such that the function expects dist to be scaled in m.
               }
            }
            dispersalKernel[i + j*gridWidth] = density;   // check whether this is the way I do the assignemt elsewhere!
          }
       }
       dispersalKernel /= dispersalKernel.sum();  //normalising

/*
       ofstream out;
       out.open("dispersalKernel.txt");
       for (size_t i = 0; i < dispersalKernel.size(); ++i)
           out << dispersalKernel[i] << " ";

       out.close() ;
 */

/*
       double sum = dispersalKernel.sum();
       double s;

       for (size_t i = 0; i < dispersalKernel.size(); i++)      // kumulative Dispersalkernel
       {
         s = dispersalKernel[i];
         dispersalKernel[i] /= sum;
         if (dispersalKernel[i] > 1.0) dispersalKernel[i] = 1.0; // catch rounding errors
         sum -= s;
       }
*/
    };


//----------------------------------------------------------------------------

  void Simulation::germination()
  {
     unsigned int arenaSize = gridWidth*gridHeight;
     size_t numspecies = plantPars.size();
     unsigned int maxNumIndividuals = ceil(maxIndividualDensity * cellArea);
     valarray<unsigned int> seedlings(numspecies);
     unsigned int candidates;
     unsigned int freePlaces;
//     unsigned int totalEstablishedSeedlings = 0;

     for (size_t i = 0; i < arenaSize; i++)
     {
        if (rocksAndPots[i] < 2 ) { //if we are not on a rock

        seedlings = 0;
        // Germination and reduction of seedbank size
        for (size_t j = 0; j < numspecies; j++)
        {
           if (seedbank[i + j*arenaSize] > 0)
           {
             // seedlings[j] = gsl_ran_binomial(theRNG, plantPars[j].germinationProbability, seedbank[i + j*arenaSize]);
             seedlings[j] = R::rbinom(seedbank[i + j*arenaSize], plantPars[j].germinationProbability);
             seedbank[i + j*arenaSize] -= seedlings[j];
           }
           else seedlings[j] = 0;

        }
        // Establishment via lottery competition
        candidates = seedlings.sum();
        if (candidates > 0)
        {
           freePlaces = maxNumIndividuals - plantDensity[i];
           if (freePlaces > candidates || freePlaces == candidates)
           {
              // all seedlings establish
              plantDensity[i] += candidates;
              for (size_t k = 0; k < numspecies; k++)
              {
                if (seedlings[k] > 0 ) initNewCohort(k, seedlings[k], plantPars,  allPlants, gridHeight, gridWidth, cellArea, barriers, simulatedTime, i);
              }
//              totalEstablishedSeedlings += candidates; //seedlings.sum();
            }
            else
            {
               // lottery competition
               plantDensity[i] = maxNumIndividuals;
               valarray<unsigned int> succesfulCandidates(numspecies);
               succesfulCandidates = 0;

               /*
               if (freePlaces < numspecies)  // this is a rough rule to decide when to use which of the two methods to generate random numbers.
               {
                 unsigned int winner;
                 unsigned int counter;
                 for (size_t k = 0; k < freePlaces; k++)
                 {
//                   winner = floor(gsl_ran_flat(theRNG, 0.0, candidates));
                   winner = floor(R::runif(0.0, candidates));
                   counter = 0;
                   for (size_t m = 0; m < numspecies; m++)
                   {
                      counter += seedlings[m];
                      if (counter > winner)
                      {
                         seedlings[m]--;
                         succesfulCandidates[m]++;
                         break;
                      }
                   }
                   candidates--;
                 }
               }
               else
               */
               {
                 // multivariate Hypergeometric Distribution
                 unsigned int sum = candidates;
                 unsigned int x, y, n;
                 n = freePlaces;
                 for (unsigned int ii = 0; ii != numspecies; ++ii)
                 {
                    y = seedlings[ii];
//                    x = gsl_ran_hypergeometric (theRNG, n, y, sum);
                    x = R::rhyper(y, sum-y, n);
                    n -= x;
                    sum -= y;
                    succesfulCandidates[ii] = x;
                  }
                  //succesfulCandidates[numspecies-1] = n;
               }
               for (size_t k = 0; k < numspecies; k++)
               {
                  if(succesfulCandidates[k] > 0) initNewCohort(k, succesfulCandidates[k], plantPars,  allPlants, gridHeight, gridWidth, cellArea, barriers, simulatedTime, i);
               }
//               totalEstablishedSeedlings += succesfulCandidates.sum();
            }
           }
        }
     }
  };

  //----------------------------------------------------------------------------

  void Simulation::seedbankmortality()
  {
    for (size_t i = 0; i < seedbank.size(); i++)
    {
    //   if (seedbank[i] > 0) seedbank[i] -= gsl_ran_binomial(theRNG, seedMortality, seedbank[i]);
       if (seedbank[i] > 0) seedbank[i] -= R::rbinom(seedbank[i], seedMortality);
    }
  };

  //----------------------------------------------------------------------------

  void Simulation::externalSeedrain()
  {
    unsigned int arenaSize = gridWidth * gridHeight;
    valarray<double> lambda = metacommunity;
    lambda *= seedrainIntensity * cellArea / metacommunity.sum();
    for (size_t i = 0; i < lambda.size(); i++)
    {
       if (lambda[i] > minLambda)
       {
          for (size_t j = 0; j < arenaSize; j++)
             // seedbank[j + i*arenaSize] += gsl_ran_poisson(theRNG, lambda[i]);
             seedbank[j + i*arenaSize] += R::rpois(lambda[i]);
       }
    }
  };

  //----------------------------------------------------------------------------

void Simulation::dispersal()
{
  // Seed production
   newSeeds = 0;
   itAllPlants = allPlants.begin();
   while (itAllPlants != allPlants.end())
   {
      (*itAllPlants).seedproduction(plantPars, newSeeds, gridWidth*gridHeight);
      itAllPlants++;
   }

   unsigned int nrow = gridHeight;
   unsigned int ncol = gridWidth;
   size_t arenaSize = nrow*ncol;

   for(size_t i = 0; i != plantPars.size(); ++i) { // species
     unsigned int newSeedsPerSpecies = 0;
     for(size_t j = 0; j != arenaSize; ++j) { // cells
       newSeedsPerSpecies += newSeeds[j + i * arenaSize];
     }
     speciesSeedproduction(year, i) = newSeedsPerSpecies;
   }

    if (year >= startEstablishmentReporting && year <= stopEstablishmentReporting) {
      newSeedsOfPastYears.pop_back();
      newSeedsOfPastYears.push_front(newSeeds);
    }
//   this function assumes that kernel is a cumulative pdf that sums to one.
    size_t targetIndex; // = 5;
    unsigned int x;

    size_t numspecies = seedbank.size() / arenaSize;

    unsigned int totalNumberOfSeeds;
    unsigned int n;
    valarray<double> aggregateKernel(0.0, arenaSize);

    for (size_t j = 0; j < numspecies; ++j) {
       // calculate the total number of seeds to disperse of this species
       totalNumberOfSeeds = 0;
       for (size_t i = 0; i < arenaSize; ++i) {
          totalNumberOfSeeds += newSeeds[i + arenaSize*j];
       }
       if (totalNumberOfSeeds > 0) {
          if (globalDispersal) {
             n = totalNumberOfSeeds;
             double p;
             for (size_t i = 0; i < arenaSize - 1; ++i) {
                  p = 1.0 / double(arenaSize - i);
//                  x = gsl_ran_binomial(theRNG, p, n);
                  x = R::rbinom(n, p);
                  n -= x;
                  seedbank[i + arenaSize*j] += x;
              }
              seedbank[arenaSize - 1 + arenaSize*j] += n;
          }
          else {


          if (totalNumberOfSeeds < arenaSize) { // for small numbers of seeds it is more efficient to disperse them one by one
//          if (totalNumberOfSeeds < 0) {
              for (size_t i = 0; i < arenaSize; ++i) {
                  n = newSeeds[i + arenaSize*j];
                  if (n > 0) {
                      double alpha, dist, drow, dcol;
                      double zeta = mu; // log(mu / cellwidth);
                      double sigma = sd; // log(sd / cellwidth);
                      for (size_t k = 0; k < n; k++)
                      {
                         // alpha = gsl_ran_flat(theRNG, 0.0, 2.0*M_PI);
                        alpha = R::runif(0.0, 2.0*M_PI) ;
//                          dist = gsl_ran_lognormal(theRNG, zeta, sigma) / cellwidth; // dist is measured in number of cells, whereas cellwidth is measured in m, and the random lognormal-distance is also in m.
                          dist = R::rlnorm(zeta, sigma) / cellwidth;
                          drow = dist * sin(alpha);
                          dcol = dist * cos(alpha);
                          // neuen Index berechnen
                          targetIndex = getTargetIndexDD(i, drow, dcol, ncol, nrow);
                          seedbank[targetIndex + arenaSize*j] ++;
                      }
                  }
              }
          }
          else {
              aggregateKernel = 0.0;
              for (size_t i = 0; i < arenaSize; ++i) {
                  n = newSeeds[i + arenaSize*j];
                  if (n > 0) {
                      double f = n / double(totalNumberOfSeeds);
                      for (size_t k = 0; k < arenaSize; ++k) {
                          aggregateKernel[k] += f * dispersalKernel[getTargetIndex(i, k, ncol, nrow)];
                      }
                  }
              }
              // multinomial verteilte Samen
              double sum = aggregateKernel.sum();
              double s;
              double p;
              n = totalNumberOfSeeds;
              for (size_t i = 0; i < arenaSize - 1; ++i) {
                  s = aggregateKernel[i];
                  p = s;
                  p /= sum;
                  if (p > 1.0) p = 1.0; // catch rounding errors
                  sum -= s;
                  // x = gsl_ran_binomial(theRNG, p, n);
                  x = R::rbinom(n, p);
                  n -= x;
                  seedbank[i + arenaSize*j] += x;
              }
              seedbank[arenaSize - 1 + arenaSize*j] += n;
          }
          }
       }
    }
}


//----------------------------------------------------------------------------

  void Simulation::evaporationRates()
  {
//    double Mw = 0.018; //Molar mass of water kg mol-1
//    double R = 8.3145; //Universal Gas konstant J mol-1 K-1;
    updateLAInew(plantPars); // get the current values of the Leaf Area Index
    // calculate the mean evapotranspiration rate for the current time period

     int evaporationVectorSize = evapotranspirationVector.size();
     double firstDayProportion = 1.0 - fmod(simulatedTime, 1.0);
     double lastDayProportion = fmod(simulatedTime + timestepLength, 1.0);

     // this funny stuff is to make sure that the evaporation sequence produced starts at
     // the beginning of the evapotranspiration vector once the simulatedTime is larger than
     // size of the evapotranspiration vector.

     int firstDay = floor(fmod(simulatedTime, double(evaporationVectorSize)));
     int lastDay = ceil(fmod(simulatedTime, double(evaporationVectorSize)));

     evapotranspiration = firstDayProportion * evapotranspirationVector[firstDay] + lastDayProportion * evapotranspirationVector[lastDay];
     if (lastDay >= firstDay)
     {
        for (int j = firstDay + 1; j < lastDay; j++) evapotranspiration += evapotranspirationVector[j];
     }
     else
     {
       for (int j = firstDay + 1; j < evaporationVectorSize; j++) evapotranspiration += evapotranspirationVector[j];
       for (int j = 0; j < lastDay; j++) evapotranspiration += evapotranspirationVector[j];
     }

     evapotranspiration /= timestepLength;        /* TODO : diese Zahlen sind immer etwas kleiner als der entsprechende Eintrag im Vector! */
//    potentialEvaporation =  evapotranspiration * exp(-kappa * leafAreaIndex) * (exp(Mw * psiE * pow(waterContent/thetaSat, -b) /(R * airTemperature)) - airHumidity) / (1.0 - airHumidity) * timestepLength;    // potential Evaporation has the dimension mm

     potentialEvaporation =  evapotranspiration * exp(-kappa * leafAreaIndex) * 0.25 * (1.0 - cos(waterContent / thetaSat * M_PI)) * (1.0 - cos(waterContent / thetaSat * M_PI));
  }


//----------------------------------------------------------------------------

void Simulation::initNewCohort(unsigned int speciesID, unsigned int cohortSize, const vector<plantParams>& pP, list<plant>& allPlants, unsigned int nrow, unsigned int ncol, double cellArea, const vector<NNcell> & barriers, double simulatedTime, unsigned int cellID)
{
   plant newPlant;
   newPlant.setspeciesid(speciesID);
   newPlant.setCellID(cellID);
   newPlant.setCohortSize(cohortSize);
   newPlant.setmass(pP[speciesID].seedlingsize);
   newPlant.setGrowth(0.0);
   newPlant.setLastTimePositiveGrowth(simulatedTime);
   newPlant.setDesiccating(false);
   newPlant.setTimeToMaturity(-1.0);

//   newPlant.initPositionVector(pP, nrow, ncol, cellArea, true, barriers);
//   newPlant.initPositionVector(pP, nrow, ncol, cellArea, false, barriers);
//   newPlant.setMaxRootIndex(pP, cellArea); // *must* be called after initialization of initPositionVector!
//   newPlant.allocation(pP, simulatedTime);
//   newPlant.allocation2(pP, cellArea);
//   newPlant.initAllocation(pP, cellArea);
//   assert(newPlant.getmass() > 0.0);

   newPlant.initAllocationNew(pP, cellArea, rootPositions);

   newPlant.setBirthday(simulatedTime);
//   newPlant.updateRoots(pP, cellArea);
//   newPlant.updateLeaves(pP, cellArea);
   newPlant.setPlantWaterContent(0.0);
   newPlant.setWaterUptake(0.0);
   newPlant.setTrueWaterUptake(0.0);
   newPlant.setAmin(pP);
//   newPlant.setThetaCrit();  //FIXME
   newPlant.setSeedmass(0.0);
   newPlant.setCohortID(cohortCounter);
   newPlant.setRecruitmentYear(-99);

   assert(newPlant.getLeafmass() > 0.0);

   allPlants.push_back(newPlant);
   cohortCounter++;
//   plantscheck();
 //  globalCounter++;
};

//----------------------------------------------------------------------------

void Simulation::readDailyRain(string file)
{
   dailyRainfall.clear();
   double x;
   ifstream infile(file.c_str());
   //infile.open(file.c_str(), ios::nocreate);
   while (infile.peek() != EOF)
   {
     infile >> x;
     dailyRainfall.push_back(x);
   }
   infile.close();
};

//----------------------------------------------------------------------------

void Simulation::readPlantInitFile(string file)      /* TODO : Falls die letzte Zeile leer ist, dann verdopple nicht einfach den letzten Eintrag. */
{
   plantInitParams2 lpip;
   plantInitPars.clear();

   string key = "none";
   if (file.find(key, 0) == 0) {     // i.e. if the string file starts with key, then do nothing
 //      file.erase(0, key.length());        // strip the key and turn the remainder into a double.
 //      double xx = StrToFloat(file);
 //      initSeedbankDensity =  xx;
   }
   else {
       ifstream infile(file.c_str());
       //infile.open(file.c_str(), ios::nocreate);
       while (infile.peek()!= EOF) {
           infile >> lpip.speciesID >> lpip.row >> lpip.col >> lpip.birthday >> lpip.mass >> lpip.cohortSize;
           plantInitPars.push_back(lpip);
       }
       infile.close();
   }
}

//----------------------------------------------------------------------------

void Simulation::drainage()
{
   drainageAmount = 0.0;
   for (size_t i = 0; i < waterContent.size(); i++)
   {
      if (waterContent[i] > thetaSat)
      {
         drainageAmount += waterContent[i] - thetaSat;
         waterContent[i] = thetaSat;
      }
   }
   cumulativeDrainage += drainageAmount;
}

//----------------------------------------------------------------------------

 void Simulation::rainfall()
 {

     size_t index = int(simulatedTime + offset) % dailyRainfall.size(); // this is to prevent access on the dailyRainfallvektor outside of its range. The sequence is simply recycled. This assumes that we work in daily timesteps!
     rain = dailyRainfall[index];

     /*
     if (repetition > 1) {
        ofstream rainout;
        string ofn = experimentPars.outputFilename + "raintest" + convert<string>(repetition) + ".dat" ;
        rainout.open(ofn.c_str(), ios::app);
        rainout << rain << " ";
        rainout.close();
     }
*/

     waterContent += rain * redistributionFactor / (1000.0 * soilDepth); // Water redistribution; dividing by 1000 scales from mm to m, and dividing by soilDepth scales to the proportion of soil volume that this water column represents
 }

//----------------------------------------------------------------------------

 void Simulation::dewfall()
 {

  // TODO: make this timestep-independent!
    // if (gsl_ran_flat(theRNG, 0.0, 365.0) < dewDays) waterContent += dewAmount / (1000.0 * soilDepth); // Water redistribution; dividing by 1000 scales from mm to m, and dividing by soilDepth scales to the proportion of soil volume that this water column represents
    if (R::runif(0.0, 365.0) < dewDays) waterContent += dewAmount / (1000.0 * soilDepth); // Water redistribution; dividing by 1000 scales from mm to m, and dividing by soilDepth scales to the proportion of soil volume that this water column represents
    //    gsl_ran_binomial(theRNG, , )

 }
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------

  size_t getTargetIndex(size_t newCenterPoint, size_t toBeTransformed,  int ncol,  int nrow)
  {
     unsigned int referenceRow = (nrow - 1)/2;
     unsigned int referenceCol = (ncol - 1)/2;

     unsigned int newCenterRow = floor(newCenterPoint / ncol) ;
     unsigned int newCenterCol = fmod(newCenterPoint, ncol);

     int deltaRow = referenceRow - newCenterRow;
     int deltaCol = referenceCol - newCenterCol;

     int preTransformationRow = floor(toBeTransformed / ncol);
     int preTransformationCol = fmod(toBeTransformed, ncol);

     int postTransformationRow = preTransformationRow + deltaRow;
     int postTransformationCol = preTransformationCol + deltaCol;

     if (postTransformationRow >= nrow) postTransformationRow = fmod(postTransformationRow, nrow);
     if (postTransformationRow < 0)
     {
        postTransformationRow = - fmod(- postTransformationRow, nrow);
        postTransformationRow += nrow;
     }
     if (postTransformationCol >= ncol) postTransformationCol = fmod(postTransformationCol, ncol);
     if (postTransformationCol < 0)
     {
        postTransformationCol = - fmod(- postTransformationCol, ncol);
        postTransformationCol += ncol;
     }

     return(postTransformationCol + ncol*postTransformationRow);

  }

//----------------------------------------------------------------------------

  size_t getTargetIndexDD(size_t newCenterPoint, double drow, double dcol,  int ncol,  int nrow)
{

     int postTransformationRow = floor(newCenterPoint / ncol);
     postTransformationRow += floor(drow + 0.5);
     int postTransformationCol = fmod(newCenterPoint, ncol);
     postTransformationCol +=  floor(dcol + 0.5);

     if (postTransformationRow >= nrow) postTransformationRow = fmod(postTransformationRow, nrow);
     if (postTransformationRow < 0)
     {
        postTransformationRow = nrow - 1 - fmod(- postTransformationRow, nrow); // die minus 1 ist wichtig, sonst kann row == nrow werden!
     }
     if (postTransformationCol >= ncol) postTransformationCol = fmod(postTransformationCol, ncol);
     if (postTransformationCol < 0)
     {
        postTransformationCol = ncol - 1 - fmod(- postTransformationCol, ncol);
     }

     return(postTransformationCol + ncol*postTransformationRow);
}

//----------------------------------------------------------------------------

double distance(double x, double y, double centerX, double centerY, unsigned int gridWidth, unsigned int gridHeight)
{

   double dx = min((x - centerX)*(x - centerX), (x + gridWidth - centerX)*(x + gridWidth - centerX));
   double dy = min((y - centerY)*(y - centerY), (y + gridHeight - centerY)*(y +gridHeight - centerY));

   return(sqrt(dx + dy));
}

//----------------------------------------------------------------------------

void Simulation::openOutputFiles()
{

/*
        if (experimentPars.outputOptions.find("0", 0) != string::npos)
        {
          string ofn = experimentPars.outputFilename + "allPlants" + convert<string>(repetition) + ".dat";
         if (experimentPars.outputOptions.find("bin", 0) != string::npos)  out.open(ofn.c_str(), ios::binary | ios::out | ios::app);
         else out.open(ofn.c_str(), ios::app);
          outputFileMap.insert(pair<string, ofstream *>("0", & out));
        }

/**/
/*
        if (experimentPars.outputOptions.find("popsize", 0) != string::npos)
        {
          string ofn = experimentPars.outputFilename + "popsize" + convert<string>(repetition) + ".dat";
         if (experimentPars.outputOptions.find("bin", 0) != string::npos)  out.open(ofn.c_str(), ios::binary | ios::out | ios::app);
         else out.open(ofn.c_str(), ios::app);

          outputFileMap.insert(pair<string, ofstream *>("popsize", & out));
        }

/**/

//        out.open("testgags.txt", ios::app);
//        outputFileMap.insert(pair<string, ofstream *>("t", & out));

}

//----------------------------------------------------------------------------

void Simulation::closeOutputFiles()
{
/*
    map<string, ofstream *>::iterator I = outputFileMap.begin();
    while( I != outputFileMap.end())
    {
       (*(*I).second).close();
       ++I;
    }
*/
}

//----------------------------------------------------------------------------


void Simulation::initEvapotranspiration(const ExperimentParameters & eP)
{
   double x;
   evapotranspirationVector.clear();

   string file = eP.dailyEvapoFile;
   ifstream infile(file.c_str());
   //infile.open(file.c_str(), ios::nocreate);
   while (infile.peek() != EOF)
   {
     infile >> x;
     evapotranspirationVector.push_back(x);
   }
   infile.close();

/*
   for (int i = 0; i < 365; i++)
   {
     x = eP.etAM0 + eP.etAM1 * cos(2.0 * M_PI / 365.0 * (double(i) - eP.etPH1)) + eP.etAM2 * cos(4.0*M_PI / 365.0 * (double(i) - eP.etPH2));
     if (x < 0.1) x = 0.1; // minimum Evapotranspiration is 0.1 mm day-1 (I think that is reasonable)
     evapotranspirationVector.push_back(x);
   }
*/

}

//----------------------------------------------------------------------------

  void Simulation::evaluateSimulation(const ExperimentParameters & eP)
  {
     long double oldSeedbank;
     long double oldSeedbankLtd;
     long double seedbank;
     long double seedbankLtd;
     long double initSeedbank;

     unsigned int maxAge;
     long double sizeAtMaturity;
     long double germinationRate;
     long double seedlingSize;
     long double massPerSeed;
     long double seedbankMortality;

     ofstream out;
     ofstream outLtd;
     long double result;
     long double resultLtd;

     vector<size_t> germinationDates;

     vector<long double> seedbankChange;  //(simulationPeriod);
     vector<long double> seedbankChangeLtd;  //(simulationPeriod);
     long double seeds;
     long double seedsLtd;
     long double cumSurvProb;
     long double populationLimit;

     vector<long double>  sizeAtMaturityVector;
     sizeAtMaturityVector.push_back(1.0E2);
     sizeAtMaturityVector.push_back(1.0E3);
     sizeAtMaturityVector.push_back(1.0E4);
     sizeAtMaturityVector.push_back(1.0E5);
     sizeAtMaturityVector.push_back(1.0E6);

     vector<long double>  germinationRates;
     germinationRates.push_back(0.1);
     germinationRates.push_back(1.0/3.0);
     germinationRates.push_back(2.0/3.0);
     germinationRates.push_back(0.9);

     vector<unsigned int> germinationMonths(12);
     for (unsigned int i = 0; i < 12; i++)
     {
        germinationMonths[i] = i;
     }

     size_t evalLag;
     long double rootmassDensity   =   plantPars[0].maxRootDensity;  // g m-2
     long double pRoot             =   plantPars[0].pRoot;

     long double arenaSize         =   1.0E6;        // m2, i.e. 1 km2
     evalLag = 0;
     maxAge = 20;
     initSeedbank = 10000.0;


     seedlingSize = plantPars[0].seedlingsize;
     massPerSeed = 1.0 / (plantPars[0].seedsPerSeedmass);
     seedbankMortality = seedMortality;

     valarray<long double> N(0.0, maxAge);
     valarray<long double> NLtd(0.0, maxAge);
     valarray<long double> M(0.0, maxAge);
     //------------------------------------------------------------------------------
     // Berechnung der Germination dates
     size_t getGerminationDate(const vector<double> & rain, int year, int month, double timestepLength);
     for (size_t year = 0; year < (numyears-1); year++)
     {
        for (size_t gM = 0; gM < germinationMonths.size(); gM++)
        {
           germinationDates.push_back(getGerminationDate(dailyRainfall, year, germinationMonths[gM], timestepLength));
        }
     }
    //-----------------------------------------------------------------------------



    string ofn = eP.outputFilename + "unLtd.txt";
    string ofnLtd = eP.outputFilename + "Ltd.txt";
    out.open(ofn.c_str());
    outLtd.open(ofnLtd.c_str());


        ofn = eP.outputFilename + "growth.txt";
        ofstream outTst;
        outTst.open(ofn.c_str());
        for (size_t ii = 0; ii < growth.size(); ii++)
           outTst << growth[ii] << " ";
        outTst.close();

    ofn = eP.outputFilename + "surv.txt";
        outTst.open(ofn.c_str());
        for (size_t ii = 0; ii < survivalProbability.size(); ii++)
           outTst << survivalProbability[ii] << " ";
        outTst.close();


    for (size_t gD = 0; gD < germinationMonths.size(); ++gD)
    {
       for (size_t gR = 0; gR < germinationRates.size(); ++ gR)
       {
          germinationRate = germinationRates[gR];
          for (size_t sAM = 0; sAM < sizeAtMaturityVector.size(); ++sAM)
          {
             sizeAtMaturity = sizeAtMaturityVector[sAM];
             populationLimit = arenaSize * rootmassDensity / (sizeAtMaturity * pRoot);
             oldSeedbank = initSeedbank;
             seedbank = oldSeedbank;
             seedbankLtd = oldSeedbank;
             oldSeedbankLtd = oldSeedbank;

             N = 0.0;
             NLtd = 0.0;
             M = 0.0;
             result = 0.0;
             resultLtd = 0.0;
             seedbankChange.clear();
             seedbankChangeLtd.clear();
             for (unsigned int year2 = 0; year2 < (numyears-2); year2++)
             {
                N[0] = germinationRate * seedbank ;
                long double popDiff = populationLimit - NLtd.sum();
                if (popDiff > 0.0)
                {
                   if (popDiff > N[0])
                   {
                      NLtd[0] = N[0];
                   }
                   else
                   {
                      NLtd[0] = popDiff;
                   }
                }
                else
                {
                   NLtd[0] = 0.0;
                }
                for (size_t i = (M.size() - 1); i > 0 ; i--)
                {
                   M[i] = M[i-1];
                }
                M[0] = seedlingSize;
                seedbank *= (1.0 - germinationRate) * (1.0 - seedbankMortality);
                seedbankLtd *= (1.0 - germinationRate) * (1.0 - seedbankMortality);
                size_t lo = germinationDates[year2*germinationMonths.size() + gD];
                size_t up = germinationDates[(year2+1)*germinationMonths.size() + gD];
                seeds = 0.0;
                seedsLtd = 0.0;
                cumSurvProb = 1.0;
                for (size_t day = lo; day < up; day++)
                {
                   cumSurvProb *= survivalProbability[day];
                   for (size_t i = 0; i < maxAge; i++)
                   {
                      M[i] *= (100.0 + growth[day]) / 100.0;  // it is assumed that the plant has a mass of 100.0 g
                      if (M[i] > sizeAtMaturity)
                      {
                         seeds += cumSurvProb * N[i] * (M[i] - sizeAtMaturity);
                         seedsLtd += cumSurvProb * NLtd[i] * (M[i] - sizeAtMaturity);
                         M[i] = sizeAtMaturity;
                      }
                   }
                }
                seedbank += seeds / massPerSeed;
                seedbankLtd += seedsLtd/ massPerSeed;
//                result += log(seedbank/oldSeedbank);
//                resultLtd += log(seedbankLtd/oldSeedbankLtd);
                seedbankChange.push_back(seedbank/oldSeedbank);
                seedbankChangeLtd.push_back(seedbankLtd/oldSeedbankLtd);
                N *= cumSurvProb;
                NLtd *= cumSurvProb;
                for (size_t i = (N.size() - 1); i > 0; i--)
                {
                   N[i] = N[i-1];
                   NLtd[i] = NLtd[i-1];
                }
                N[0] = 0.0;
                NLtd[0] = 0.0;
                oldSeedbank = seedbank;
                oldSeedbankLtd = seedbankLtd;
             }
             for (size_t i = evalLag; i < seedbankChange.size(); i++)
             {
                result += log(seedbankChange[i]);
                resultLtd += log(seedbankChangeLtd[i]);
             }
             out << result << " ";
             outLtd << resultLtd << " ";
          } // sizeAtMaturity
       } // germination Rates
    } // germination Dates
    out.close();
    outLtd.close();
 }

//-----------------------------------------------------------------------------

    size_t getGerminationDate(const vector<double> & rain, int year, int month, double timestepLength)
   {
      size_t monthOffset[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
      size_t start = year * 365 + monthOffset[month];
      bool  rainFound = false;
      while (!rainFound)
      {
         if (rain[start] > 0.0)
         {
           rainFound = true;
         }
         else
         {
           start++;
         }
      }
      // nun muss der Index in der dreier-Z?hlweise gefunden werden:
      return(start/timestepLength+1);
   }

//-----------------------------------------------------------------------------

  void Simulation::grazing(const vector<plantParams> & pP)
  {
      if (grazingAnnualAmount > 0.00001 && grazingFrequency > 0.00001) {
          // is this a grazing day?
          bool grazingDay = false;

          unsigned int dayInYear = fmod( simulatedTime, 365.0 );
          if ( dayInYear >= startGrazingPeriod && dayInYear <= endGrazingPeriod ) {
             double probGrazingDay = grazingFrequency / (endGrazingPeriod - startGrazingPeriod) * 365.0 / timeStepsPerYear;
             probGrazingDay = min(probGrazingDay, 1.0);
             // if ( gsl_ran_flat(theRNG, 0.0, 1.0) < probGrazingDay ) grazingDay = true;
             if ( R::runif(0.0, 1.0) < probGrazingDay ) grazingDay = true;
          }

          if ( grazingDay ) {
              double grazingAmountPerEvent = grazingAnnualAmount * gridWidth * gridHeight * cellArea / grazingFrequency; // total amount of biomass in g that is grazed on average per grazing event

              // calculate total amount of currently available biomass:
              double abovegroundBiomass = 0.0;
              itAllPlants = allPlants.begin();
              while (itAllPlants != allPlants.end()) {
                  abovegroundBiomass += (*itAllPlants).getStoragemass() + (*itAllPlants).getLeafmass();
                 itAllPlants++;
              }

              if (abovegroundBiomass > 0.000001) {
                 double averageGrazedProportion = grazingAmountPerEvent / abovegroundBiomass;
                 averageGrazedProportion = min(0.9, averageGrazedProportion); // not more than 90% of biomass is grazed per event

              // grazing
                 itAllPlants = allPlants.begin();
                 while (itAllPlants != allPlants.end()) {
                     // unsigned int grazedProportion = gsl_ran_binomial (theRNG, averageGrazedProportion, 20);
                   unsigned int grazedProportion = R::rbinom(20, averageGrazedProportion);
                   (*itAllPlants).setmass((*itAllPlants).getmass() - double(grazedProportion) / 20.0 * ((*itAllPlants).getLeafmass() + (*itAllPlants).getStoragemass()));
                     (*itAllPlants).setStoragemass(double(20 - grazedProportion) / 20.0 * (*itAllPlants).getStoragemass());
                     (*itAllPlants).setLeafmass(double(20 - grazedProportion) / 20.0 * (*itAllPlants).getLeafmass());
                     (*itAllPlants).setPlantWaterContent(double(20 - grazedProportion) / 20.0 * (*itAllPlants).getPlantWaterContent());
                     (*itAllPlants).setAmin(pP);
                     itAllPlants++;
                 }
              }
          }
      }
  }

//-----------------------------------------------------------------------------

void Simulation::initLeafFraction(const vector<plantParams> & pP) {

    // Calculate maximum mass of leaf-structural biomass.
    double a = pP[0].leafStructureAllometryFactor;
    double b = pP[0].leafStructureAllometryPower;

    double maxMass = 0.0;
    for (size_t i = 0; i < pP.size(); ++i) {
        if (pP[i].maxSize * pP[i].pLeaf > maxMass) maxMass = pP[i].maxSize * pP[i].pLeaf;
    }

//    vector <double> masses;
    leafFraction.clear();
//    deltaLeafFraction.clear();
//    masses.clear();
//    masses.push_back(0.0);
//    deltaLeafFraction.push_back(1.0);
    double x;
//    double dx;
//    x = optimalLeaffraction (pP[0].seedlingsize, a, b);
//    leafFraction.push_back(x);


    double minLogMass = floor(log(pP[0].seedlingsize * 0.001));
    if (maxMass < pP[0].seedlingsize) maxMass = pP[0].seedlingsize;
    double maxLogMass = ceil(log(maxMass));


    for (double M = minLogMass; M < maxLogMass + 2.0; ++M) {
        x =  optimalLeaffraction (exp(M), a, b);
        leafFraction.push_back(x);
    }

/*
    ofstream out;
    out.open("initLeafFraction.txt");
    for (size_t i = 0; i < leafFraction.size(); ++i) {
       out << leafFraction[i] << " ";
    }
    out.close();
 */
}

//------------------------------------------------------
   double f1 (double x, double M, double a, double b) {
       double opt = (x * M - a * pow((1.0 - x) * M, b)) / (x * M);
//       return(opt);
       return(opt * opt);
   }

//------------------------------------------------------

   double deriv1 (double x, double M, double a, double b) {
        double h1 = a * pow((1.0 - x) * M, b);
        return( 2 * (x * M - h1) * ( M + h1 * b / (1.0 - x)));
   }

//------------------------------------------------------

double optimalLeaffraction (double M, double a, double b){
//        double lowerBound = 0.0, upperBound = 1.0;
        double x = 0.5;
        double dx = 0.5;
        double tol = 0.00000000001;
//        double diff;
        int maxiter = 50;
        int iter = 0;
        bool converged = false;

        while ( !converged && iter < maxiter ) {
            ++iter;
            dx *= 0.5;
            if ( f1(x, M, a, b) < tol ) {
                converged = true;
            }
            else {
                 if (deriv1(x, M, a, b) > 0) {
                      x -= dx;
                 }
                 else {
                      x += dx;
                 }
            }
        }
        return(x);
}


//void Simulation::initRootPositions(const vector<plantParams> & pP)
  void Simulation::initRootPositions(const vector<plantParams> & pP,  unsigned int nrow, unsigned int ncol, const vector<int> & rocksAndPots){
    // determine maximum area covered by roots of a single plant:
    double maxRootArea = 0.0;
    for (size_t i = 0; i < pP.size(); ++i) {
        if (pP[i].maxSize * pP[i].pRoot / pP[i].maxRootDensity > maxRootArea) maxRootArea = pP[i].maxSize * pP[i].pRoot / pP[i].maxRootDensity;
    }
    unsigned int maxNumCells = ceil(maxRootArea / cellArea);
    initPositionVectorNew(rootPositions, nrow, ncol, maxNumCells, rocksAndPots);

}

  void Simulation::initLeafPositions(const vector<plantParams> & pP,  unsigned int nrow, unsigned int ncol, const vector<int> & rocksAndPots){

    // determine maximum area covered by roots of a single plant:
    double maxLeafArea = 0.0;
    for (size_t i = 0; i < pP.size(); ++i) {
        if (pP[i].maxSize * pP[i].pLeaf / pP[i].maxLeafDensity > maxLeafArea) maxLeafArea = pP[i].maxSize * pP[i].pLeaf / pP[i].maxLeafDensity;
    }
    unsigned int maxNumCells = ceil(maxLeafArea / cellArea);
    initPositionVectorNew(leafPositions, nrow, ncol, maxNumCells, rocksAndPots);
}



void Simulation::initPositionVectorNew( vector< vector<unsigned int> > & positionVector, unsigned int nrow, unsigned int ncol, unsigned int maxNumCells, const vector<int> & rocksAndPots) {
  valarray<bool> selected(false, rocksAndPots.size());
  vector<unsigned int> cellsAccepted;
  vector<unsigned int> candidates;
  unsigned int newAcceptedCell;
  positionVector.clear();

  for (size_t i = 0; i < rocksAndPots.size(); ++i) { // Loop over all cells in the model arena
       selected = false;         // initially, no cells are selected
       cellsAccepted.clear();
       if (rocksAndPots[i] == 0) { // if the current cell is not a rock
           cellsAccepted.push_back(i);
           selected[i] = true;
           candidates.clear();
           if (maxNumCells > 1) getCandidatesNew(candidates, selected, i, nrow, ncol, rocksAndPots);
           while (cellsAccepted.size() < maxNumCells && candidates.size() > 0) // noch ein weiteres Abbruchkriterium finden
           {
               candidates.clear();
               for(size_t j = 0; j < cellsAccepted.size(); j++) {
                   getCandidatesNew(candidates, selected, cellsAccepted[j], nrow, ncol, rocksAndPots);
               }
               if (candidates.size() > 0) {
                   newAcceptedCell = getNearestCandidateNew(candidates, i, nrow, ncol);
                   selected[newAcceptedCell] = true;
                   cellsAccepted.push_back(newAcceptedCell);
               }
           }
       }
       else {
          if (rocksAndPots[i] == 1) cellsAccepted.push_back(i);
          if (rocksAndPots[i] == 2) cellsAccepted.push_back(rocksAndPots.size());
       }
       positionVector.push_back(cellsAccepted);
  }
}


void getCandidatesNew(vector<unsigned int> & candidates, const valarray<bool> & selected, unsigned int cellID, unsigned int nrow, unsigned int ncol, const vector<int> & barriers)
{
   vector<unsigned int> fourNeighbors(4);
   unsigned int up, le, ri, lo;
   if (cellID < ncol) up = cellID + (nrow - 1)*ncol; // if cellID is in the first row
   else up = cellID - ncol;

   if ( floor(fmod(cellID, ncol)) == 0 ) le = cellID + ncol - 1;
   else le = cellID - 1;

   if ( floor(fmod(cellID + 1, ncol)) == 0) ri = cellID - ncol + 1;
   else ri = cellID + 1;

   if (cellID + ncol + 1 > ncol*ncol) lo = cellID - (nrow - 1)*ncol;
   else lo = cellID + ncol;


   fourNeighbors[0] = up;
   fourNeighbors[1] = le;
   fourNeighbors[2] = ri;
   fourNeighbors[3] = lo;

   for (size_t i = 0; i < fourNeighbors.size(); ++i)
   {
     if ( selected[fourNeighbors[i]] == false && barriers[fourNeighbors[i]] == 0 )
     {
       candidates.push_back(fourNeighbors[i]);
     }
   }
};

unsigned int getNearestCandidateNew(vector<unsigned int> & candidates, unsigned int cellID, unsigned int nrow, unsigned int ncol)
{
  unsigned int minDistance;
  unsigned int newDistance;
  unsigned int nearestCandidate;
  minDistance = getSquareDistance(cellID, candidates[0], nrow, ncol);
  nearestCandidate = candidates[0];
  for (size_t i = 1; i < candidates.size(); i++)
  {
     newDistance = getSquareDistance(cellID, candidates[i], nrow, ncol);
     if (newDistance < minDistance)
     {
        minDistance = newDistance;
        nearestCandidate = candidates[i];
     }
  }
  return(nearestCandidate);
};


double getSquareDistanceNew(unsigned int cellID1, unsigned int cellID2, unsigned int nrow, unsigned int ncol)
{
   double row1, row2, col1, col2;
   double drow, dcol;
   double maxrow, minrow, maxcol, mincol;

   row1 = floor(double(cellID1) / double(ncol));
   row2 = floor(double(cellID2) / double(ncol));
   col1 = fmod(cellID1, ncol);
   col2 = fmod(cellID2, ncol);

   maxrow = max(row1, row2);
   minrow = min(row1, row2);
   maxcol = max(col1, col2);
   mincol = min(col1, col2);

   drow = min((row1 - row2)*(row1 - row2), (minrow + nrow - maxrow)*(minrow + nrow - maxrow));
   dcol = min((col1 - col2)*(col1 - col2), (mincol + ncol - maxcol)*(mincol + ncol - maxcol));

   return(drow + dcol);

}
//----------------------------------------------------------------------------

   void Simulation::evaporationAndPlantWaterUptakeNew()
   {

      /* TODO : Diese Konstanten sollten am Anfang der Simulation berechnet werden */
      double MMToM3M3 = 0.001 / soilDepth; // Umrechnungsfaktor von mm zu m3 m-3
      double M3M3ToMM = 1.0 / MMToM3M3; // 1000.0 * soilDepth; // check this
      double MMToG = 1000.0 * cellArea;
      double GToMM = 0.001 / cellArea;

      evaporationRates();

      // Calculation of maximum water uptake rates based on soil water content
      for (size_t i = 0; i < maxWaterUptakeRate.size(); ++i) {
          if (waterContent[i] > thetaCrit)
               maxWaterUptakeRate[i] = exp( uptakeSlope * (0.01 / (thetaSat - thetaCrit) - 0.01 /(waterContent[i] - thetaCrit)) );
          else maxWaterUptakeRate[i] = 0.0;
      }

      plantWaterDemand = 0.0; // in welche Funktion kann man das einbauen=
      itAllPlants = allPlants.begin();
      while (itAllPlants != allPlants.end())
      {
         (*itAllPlants).waterorderNew(plantWaterDemand, plantPars, maxWaterUptakeRate, soilDepth, timestepLength, thetaCrit, thetaSat, evapotranspiration, cellArea, rootPositions);
         itAllPlants++;
      }

      for (size_t i = 0; i < plantWaterDemand.size(); i++)
      {
           realizedPlantWaterDemand[i] = plantWaterDemand[i] * GToMM; // realizedPlantWaterDemand is in mm

         // 2. Partition available water among evaporation and water uptake (i.e. transpiration)
           double availWater = waterContent[i] - thetaCrit; // das kann negativ sein!

           if (availWater > (realizedPlantWaterDemand[i] + potentialEvaporation[i]) * MMToM3M3)
           {
              // jeder bekommt was er will; realizedPlantWaterDemand does not change.
              waterContent[i] -= (realizedPlantWaterDemand[i] + potentialEvaporation[i]) * MMToM3M3;
              actualEvaporation[i] = potentialEvaporation[i];
           }
           else
           {
              if (availWater > 0.0)  // some water uptake for plants.
              {
                 // how much water did evaporate after thetaCrit was reached?
                 double evapo2 = potentialEvaporation[i] * (1.0 - availWater * M3M3ToMM / (realizedPlantWaterDemand[i] + potentialEvaporation[i])) ; /* TODO : Fix this */

                 realizedPlantWaterDemand[i] = M3M3ToMM * availWater * (realizedPlantWaterDemand[i] /(realizedPlantWaterDemand[i] + potentialEvaporation[i]));
                 actualEvaporation[i] = M3M3ToMM * availWater - realizedPlantWaterDemand[i];
                 waterContent[i] -= availWater;
                 availWater = thetaCrit;
                 if (availWater > evapo2 * MMToM3M3)
                 {
                    actualEvaporation[i] += evapo2;
                    waterContent[i] -= evapo2 * MMToM3M3;
                 }
                 else
                 {
                   waterContent[i] = 0.0;
                   actualEvaporation[i] += availWater * M3M3ToMM;
                 }
              }
              else // no water uptake for plants.
              {
                realizedPlantWaterDemand[i] = 0.0;
                availWater = waterContent[i];
                if (availWater > potentialEvaporation[i] * MMToM3M3)
                {
                  actualEvaporation[i] = potentialEvaporation[i];
                  waterContent[i] -= actualEvaporation[i] * MMToM3M3;
                }
                else
                {
                   waterContent[i] = 0.0;
                   actualEvaporation[i] = max(0.0, availWater * M3M3ToMM); // eigentlich ist diese Maximumsbildung ?berfl?ssig, availWater sollte nicht kleiner Null werden.
                }
              }
           }
      }
    realizedPlantWaterDemand *= MMToG;

    cumulativePlantWaterUptake += realizedPlantWaterDemand.sum();
    cumulativeEvaporation += actualEvaporation.sum();

    itAllPlants = allPlants.begin();
    while (itAllPlants != allPlants.end())
    {
       (*itAllPlants).wateruptakeNew(plantWaterDemand, plantPars, realizedPlantWaterDemand, cellArea, rootPositions, maxWaterUptakeRate);
       itAllPlants++;
    }
   }

//----------------------------------------------------------------------------

void Simulation::updateLAInew(const vector<plantParams>& pP)
{
   list<plant>::iterator itPlants;
   itPlants = allPlants.begin();
   leafAreaIndex = 0.0;

   plantParams lpp;
   double numberOfLeafCells;
   double LAIInLastCell;
   double cohortSize;
   double plantLeafmass;
   unsigned int availableLeafCells;
   unsigned int plantCellID;

   while(itPlants != allPlants.end())
   {
       lpp = pP[(*itPlants).getSpeciesID()];
       plantCellID = (*itPlants).getCellID();
       plantLeafmass = (*itPlants).getLeafmass();
       numberOfLeafCells = ceil(plantLeafmass / (cellArea * lpp.maxLeafDensity));
       cohortSize = (*itPlants).getCohortSize();

       availableLeafCells = (leafPositions[plantCellID]).size();
       if (numberOfLeafCells > availableLeafCells) {
           for (size_t i = 0; i < availableLeafCells; ++i) {
//               assert(leafPositions[plantCellID][i] >= 0);
               assert(leafPositions[plantCellID][i] < leafAreaIndex.size());

               leafAreaIndex[leafPositions[plantCellID][i]] += lpp.LAIopt * cohortSize;
           }
       }
       else {
           for (size_t i = 0; i < (numberOfLeafCells - 1); ++i) {
//               assert(leafPositions[plantCellID][i] >= 0);
               assert(leafPositions[plantCellID][i] < leafAreaIndex.size());
               leafAreaIndex[leafPositions[plantCellID][i]] += lpp.LAIopt * cohortSize;
           }
           LAIInLastCell = (plantLeafmass - ((numberOfLeafCells - 1) * cellArea * lpp.maxLeafDensity)) /  lpp.maxLeafDensity * lpp.LAIopt;
           if (numberOfLeafCells < 1) {
           assert(numberOfLeafCells > 0);
           }
//           assert(plantCellID >= 0);
           assert(plantCellID < leafPositions.size());
//           assert(leafPositions[plantCellID][(numberOfLeafCells - 1)] >= 0);
           assert(leafPositions[plantCellID][(numberOfLeafCells - 1)] < leafAreaIndex.size());

           leafAreaIndex[leafPositions[plantCellID][(numberOfLeafCells - 1)]] = LAIInLastCell * cohortSize;
       }
       itPlants++;
   }
}

//----------------------------------------------------------------------------

void Simulation::adjustRedistributionMap()
{
    // Purpose: Rainwater that falls on rocks is redistributed, proportional to the redistributionFactors.
    for (size_t i = 0; i < redistributionFactor.size(); ++i) {
        if (rocksAndPots[i] == 2) redistributionFactor[i] = 0.0;
    }
    double sum = redistributionFactor.sum() / double(redistributionFactor.size());
    assert(sum > 0.0);
    redistributionFactor /= sum;
}


//-----------------------------------------------------------------------------

    size_t getGerminationDateNew(const vector<double> & rain, int year, unsigned int germinationDate, unsigned int rainPeriod, double rainfallThreshold)
   {
      size_t maxRainIndex = rain.size();
//      size_t monthOffset[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
//      size_t start = year * 365 + monthOffset[month];
      size_t start = year * 365 + germinationDate;
      double rainSum;
      bool  rainFound = false;
      while (!rainFound && start < maxRainIndex)
      {
         rainSum = 0.0;
         for (size_t i = 0; i < rainPeriod; ++i) {
            rainSum += rain[(start + i) % maxRainIndex];
         }
         if (rainSum > rainfallThreshold)
         {
           rainFound = true;
         }
         else
         {
           start++;
         }
      }

      if (rainFound == false) start =  year * 365 + germinationDate;
      return(start);
   }

//-----------------------------------------------------------------------------


void Simulation::initGerminationDates(const vector<plantParams> & pP)
{
   germinationDates.clear();
   unsigned int germinationDate = pP[0].germinationDate;
   unsigned int rainPeriod = pP[0].rainPeriod;
   double rainfallThreshold = pP[0].rainfallThreshold;
   unsigned int date;

   for (size_t i = 0; i < numyears; ++i) {
       date = getGerminationDateNew(dailyRainfall, i, germinationDate, rainPeriod, rainfallThreshold);
       germinationDates.push_back(date);
   }
}

//-----------------------------------------------------------------------------

void Simulation::initDispersalDates(const vector<plantParams> & pP)
{
   dispersalDates.clear();
   unsigned int dispersalDate = pP[0].dispersalDate;
   unsigned int date;
   for (size_t i = 0; i < numyears; ++i) {
       date = i * 365 + dispersalDate;
       dispersalDates.push_back(date);
   }
}

void readVector(vector<int>& x, string file, string key)
{

  if (file.find(key, 0) == 0)      // i.e. if the string file starts with key
  {
    file.erase(0, key.length());        // strip the key and turn the remainder into a double.
    int xx = StrToInt(file);
    for (size_t i = 0; i < x.size(); ++i) x[i] = xx;
  }
  else
  {
    size_t c = 0;
    size_t xSize = x.size();
    ifstream infile(file.c_str());
    // infile.open(file.c_str(), ios::nocreate);
    while (infile.peek()!= EOF && c < xSize)
    {
      infile >> x[c];
      c++;
    }
  }
}

void readValarray(valarray<double>& x, string file, string key)
{

  if (file.find(key, 0) == 0)      // i.e. if the string file starts with key
  {
    file.erase(0, key.length());        // strip the key and turn the remainder into a double.
    double xx = StrToFloat(file);
    x =  xx;
  }
  else
  {
    size_t c = 0;
    size_t xSize = x.size();
    ifstream infile(file.c_str());
    //infile.open(file.c_str(), ios::nocreate);
    while (infile.peek()!= EOF && c < xSize)
    {
      infile >> x[c];
      c++;
    }
  }
}

// void checkObserver(int year){};

void checkObserver(int simuTimestep, Simulation* simu)    // TODO : check population size spikes after germination
{
/*
  Simulation& testSimu(*simu);
  const list<plant> & plants = testSimu.getPlantList();
  list<plant>::const_iterator itPlants;
  const vector<plantParams> & pP = testSimu.getPlantPars();
  const valarray<double> & soilWaterContent = testSimu.getSoilWaterContent();
  valarray<int> individualsPerSpecies(pP.size());
  individualsPerSpecies = 0;
  if (rc == 0) {
    rc = eP.reportCycle + eP.samplingDate;
    //       const map<string, ofstream *> & outputFileMap = testSimu.getOutputFileMap();
    double tt = testSimu.getCurrentYear() * testSimu.getTimeStepsPerYear() + testSimu.getCurrentTimestep();
    if (eP.outputOptions.find("0", 0) != string::npos)
    {
      itPlants = plants.begin();
      while (itPlants != plants.end()) {
        (*itPlants).writePlantShort(& gOutPlantShort,  tt);
        itPlants++;
      }
    }

    if (eP.outputOptions.find("long", 0) != string::npos)
    {
      itPlants = plants.begin();
      while (itPlants != plants.end()) {
        (*itPlants).writePlant(& gOutPlantLong, tt);
        itPlants++;
      }
    }


    if (eP.outputOptions.find("popsize", 0) != string::npos)
    {
      individualsPerSpecies = 0;
      valarray<int> bigIndividualsPerSpecies(pP.size());
      bigIndividualsPerSpecies = 0;
      valarray<int> matureIndividualsPerSpecies(pP.size());
      matureIndividualsPerSpecies = 0;
      valarray<int> recruits(pP.size());
      recruits = 0;
      valarray<double> meanAge(pP.size());
      meanAge = 0.0;
      valarray<double> meanTimeToMaturity(pP.size());
      meanTimeToMaturity = 0.0;

      int year = testSimu.getCurrentYear();

      itPlants = plants.begin();
      while (itPlants != plants.end()) {
        plant lp = (*itPlants);
        individualsPerSpecies[lp.getspeciesid()] += lp.getCohortSize();
        if (lp.getmass() > 0.1) {
          bigIndividualsPerSpecies[lp.getspeciesid()] += lp.getCohortSize();
          meanAge[lp.getspeciesid()] += (testSimu.getSimulatedTime() - lp.getbirthday()) * lp.getCohortSize();
        }
        //               if (lp.getmass() > 0.999 * pP[lp.getspeciesid()].maxSize ) {
        if (lp.getTimeToMaturity() > 0.0 ) {
          matureIndividualsPerSpecies[lp.getspeciesid()] += lp.getCohortSize();
          meanTimeToMaturity[lp.getspeciesid()] += lp.getTimeToMaturity() * lp.getCohortSize();
        }

        if (lp.getRecruitmentYear() == (year - 1))
          recruits[lp.getspeciesid()] += lp.getCohortSize();
        itPlants++;
      }
      for (size_t i = 0; i < individualsPerSpecies.size(); i++) {
        gOutPopsize << individualsPerSpecies[i] << " ";
        gOutPopsizeBig << bigIndividualsPerSpecies[i] << " ";
        gOutPopsizeMature << matureIndividualsPerSpecies[i] << " ";
        gOutRecruitment << recruits[i] << " ";
        if (bigIndividualsPerSpecies[i] > 0)
          gOutMeanAge << meanAge[i] / bigIndividualsPerSpecies[i] / double(testSimu.getTimeStepsPerYear()) << " ";
        else
          gOutMeanAge << "NA ";
        if (matureIndividualsPerSpecies[i] > 0)
          gOutMeanTimeToMaturity << meanTimeToMaturity[i] / double(matureIndividualsPerSpecies[i]) << " ";
        else
          gOutMeanTimeToMaturity << "NA ";

      }
      gOutPopsize << "\n";
      gOutPopsizeBig << "\n";
      gOutPopsizeMature << "\n";
      gOutRecruitment << "\n";
      gOutMeanAge << "\n";
      gOutMeanTimeToMaturity << "\n";
    }


    if (eP.outputOptions.find("theta", 0) != string::npos)
    {
      for (size_t i = 0; i < soilWaterContent.size()-1; i++)
      {
        gOutTheta << soilWaterContent[i] << " ";
      }
      gOutTheta << soilWaterContent[soilWaterContent.size()-1] << "\n";
    }

    if (eP.outputOptions.find("seeds", 0) != string::npos)
    {

      size_t numSpecies = pP.size();
      size_t arenaSize = soilWaterContent.size();
      valarray<unsigned int> seedsPerSpecies(0.0, numSpecies);
      const valarray<unsigned int> & seedbank = testSimu.getSeedbank();
      const vector<int> & rocksAndPots = testSimu.getRocksAndPots();

      for (size_t i = 0; i < arenaSize; i++)
      {
        if (rocksAndPots[i] < 2 ) { //if we are not on a rock
          for (size_t j = 0; j < numSpecies; j++)
          {
            seedsPerSpecies[j] += seedbank[i + j*arenaSize];
          }
        }
      }

      for (size_t i = 0; i < seedsPerSpecies.size()-1; i++)
      {
        gOutSeedbanksize << seedsPerSpecies[i] << " ";
      }
      gOutSeedbanksize << seedsPerSpecies[seedsPerSpecies.size()-1] << "\n";
    }


    if (eP.outputOptions.find("seedbank", 0) != string::npos)
    {

      size_t numSpecies = pP.size();
      size_t arenaSize = soilWaterContent.size();
      const valarray<unsigned int> & seedbank = testSimu.getSeedbank();
      for (size_t i = 0; i < arenaSize; i++)
      {
        for (size_t j = 0; j < numSpecies - 1; j++)
        {
          gOutSeedbank << seedbank[i + j*arenaSize] << " ";
        }
        gOutSeedbank << seedbank[i + (numSpecies - 1)*arenaSize] << "\n";
      }
    }
    if (eP.outputOptions.find("cover", 0) != string::npos)
    {
      valarray<double> coverage(0.0, pP.size());
      itPlants = plants.begin();
      while (itPlants != plants.end()) {
        coverage[(*itPlants).getspeciesid()] += (*itPlants).getCohortSize() * (*itPlants).getLeafmass() / pP[(*itPlants).getspeciesid()].maxLeafDensity;
        itPlants++;
      }
      coverage /= (eP.gridWidth * eP.gridHeight * eP.cellwidth * eP.cellwidth);

      for (size_t i = 0; i < coverage.size()-1; i++)
      {
        gOutCoverage << coverage[i] << " ";
      }
      gOutCoverage << coverage[coverage.size()-1] << "\n";
    }

    if (eP.outputOptions.find("biomass", 0) != string::npos)
    {
      valarray<double> biomass(0.0, pP.size());
      itPlants = plants.begin();
      while (itPlants != plants.end()) {
        biomass[(*itPlants).getspeciesid()] += (*itPlants).getCohortSize() * (*itPlants).getmass();
        itPlants++;
      }

      for (size_t i = 0; i < biomass.size()-1; i++)
      {
        gOutBiomass << biomass[i] << " ";
      }
      gOutBiomass << biomass[biomass.size()-1] << "\n";
    }

    if (eP.outputOptions.find("swc", 0) != string::npos)
    {
      size_t arenaSize = soilWaterContent.size();
      double ctr = 0.0;
      double swc = 0.0;
      const vector<int> & rocksAndPots = testSimu.getRocksAndPots();
      for (size_t i = 0; i < arenaSize; i++)
      {
        if (rocksAndPots[i] < 2 ) { //if we are not on a rock
          swc += soilWaterContent[i];
          ++ctr;
        }
      }
      gOutSWC << swc / ctr << " ";
    }

    if (eP.outputOptions.find("wateruptake", 0) != string::npos)
    {
       gOutWaterUptake << testSimu.getCumulativePlantWaterUptake() / 1000.0 << " " << testSimu.getCumulativeEvaporation() * eP.cellwidth * eP.cellwidth  << " " << testSimu.getCumulativeDrainage() * eP.soilDepth * eP.cellwidth * eP.cellwidth * eP.gridWidth * eP.gridHeight << "\n";
    }

    if (eP.outputOptions.find("evaporation", 0) != string::npos)
    {
      double swc = 0.0;
      double ctr = 0.0;
      size_t arenaSize = soilWaterContent.size();
      const vector<int> & rocksAndPots = testSimu.getRocksAndPots();
      const valarray<double> & actualEvaporation = testSimu.getActualEvaporation();
      for (size_t i = 0; i < arenaSize; i++)
      {
        if (rocksAndPots[i] < 2 ) { //if we are not on a rock
          swc += actualEvaporation[i];
          ++ctr;
        }
      }
      gOutEvaporation << swc/ctr * eP.cellwidth * eP.cellwidth * eP.gridWidth * eP.gridHeight << " " << testSimu.getDrainageAmount() * eP.soilDepth * eP.cellwidth * eP.cellwidth * eP.gridWidth * eP.gridHeight<< "\n";
    }

    if (eP.outputOptions.find("swc", 0) != string::npos)
    {
      const valarray<double> & totalDeaths = testSimu.getTotalDeaths();
      const valarray<double> & droughtDeaths = testSimu.getDroughtDeaths();
      for (size_t i = 0; i < totalDeaths.size(); i++)
      {
        gOutTotalDeaths << totalDeaths[i] << " ";
        if (totalDeaths[i] > 0)
          gOutDroughtDeaths << droughtDeaths[i] / totalDeaths[i] << " ";
        else
          gOutDroughtDeaths << "NA ";
      }
      gOutTotalDeaths <<  "\n";
      gOutDroughtDeaths <<  "\n";
    }
  }
  rc--;
 */
}

void Simulation::readPlantPars(string plantParFile)
{
  bool allOK = true;
  plantPars.clear();
  plantParams pP;
  ifstream infile(plantParFile.c_str());
  //infile.open(plantParFile.c_str(), ios::nocreate);
  int start = 0;
  while ( allOK )
  {
    readPar(pP.baseMortalityRate, infile, "baseMortalityRate", start, allOK, " = ");
    readPar(pP.desiccationAlpha, infile, "desiccationAlpha", start, allOK, " = ");
    readPar(pP.desiccationGamma, infile, "desiccationGamma", start, allOK, " = ");
    readPar(pP.germinationProbability, infile, "germinationProbability", start, allOK, " = ");
    readPar(pP.HalfSaturationTheta, infile, "HalfSaturationTheta", start, allOK, " = ");
    readPar(pP.seedsPerSeedmass, infile, "seedsPerSeedmass", start, allOK, " = ");
    readPar(pP.maxAge, infile, "maxAge", start, allOK, " = ");
    readPar(pP.maxRootDensity, infile, "maxRootDensity", start, allOK, " = ");
    readPar(pP.maxSize, infile, "maxSize", start, allOK, " = ");
    readPar(pP.maxSpecificTranspirationRate, infile, "maxSpecificTranspirationRate", start, allOK, " = ");
    readPar(pP.maxUptakeRate, infile, "maxUptakeRate", start, allOK, " = ");
    readPar(pP.pRoot, infile, "pRoot", start, allOK, " = ");
    readPar(pP.pStorage, infile, "pStorage", start, allOK, " = ");
    readPar(pP.pLeaf, infile, "pLeaf", start, allOK, " = ");
    //      readPar(pP.psiCrit, infile, "psiCrit", start, allOK, " = ");
    readPar(pP.respirationRate, infile, "respirationRate", start, allOK, " = ");
    readPar(pP.rGrowth, infile, "rGrowth", start, allOK, " = ");
    readPar(pP.seedlingsize, infile, "seedlingsize", start, allOK, " = ");
    readPar(pP.succulenceFactor, infile, "succulenceFactor", start, allOK, " = ");
    readPar(pP.TranspirationSlope, infile, "TranspirationSlope", start, allOK, " = ");
    readPar(pP.waterUseEfficiency, infile, "waterUseEfficiency", start, allOK, " = ");
    readPar(pP.maxLeafDensity, infile, "maxLeafDensity", start, allOK, " = ");
    readPar(pP.leafStructureAllometryPower, infile, "leafStructureAllometryPower", start, allOK, " = ");
    readPar(pP.leafStructureAllometryFactor, infile, "leafStructureAllometryFactor", start, allOK, " = ");
    readPar(pP.LAIopt, infile, "LAIopt", start, allOK, " = ");
    readPar(pP.germinationDate, infile, "germinationDate", start, allOK, " = ");
    readPar(pP.rainPeriod, infile, "rainPeriod", start, allOK, " = ");
    readPar(pP.rainfallThreshold, infile, "rainfallThreshold", start, allOK, " = ");
    readPar(pP.dispersalDate, infile, "dispersalDate", start, allOK, " = ");
    //      readPar(pP.uptakeSlope, infile, "uptakeSlope", start, allOK, " = ");
    //      pP.uptakeSlope = 4.3;

    /*      if ( allOK ) {
     readPar(pP.maxAgeGamma, infile, "maxAgeGamma", start, allOK, " = ");
     if ( allOK == false ) {
     pP.maxAgeGamma = 5.0;
     allOK = true;
     }
    }
     */
    pP.maxAgeGamma = 5.0;


    if ( allOK )
    {
      plantPars.push_back(pP);
      start = infile.tellg();
    }
    else {
//      writeErrorMessage("error.txt", "Could not read plant parameter file", plantPars.size());
    }
  }
  infile.close();
}

//----------------------------------------------------------------------------------

template<class T>
void readPar(T& x, ifstream & in, string key, int start, bool& found, string sep)
{
  //  T x;
  bool foundit = false;
  char str[4096];          // necessary, because borland has not yet implemented getline for stl strings
  string line;
  string s = key + sep;
  size_t currentpos, nextlinepos, findpos;
  in.seekg(start);

  while (in.peek() != EOF && !foundit)
  {
    currentpos = in.tellg();
    in.getline(str, 4096); // necessary, because borland has not yet implemented getline for stl strings
    line = str;            // necessary, because borland has not yet implemented getline for stl strings
    findpos = line.find(s, 0);
    if (findpos != string::npos) //i.e. string s was found in line
    {
      foundit = true;
      nextlinepos = in.tellg();
      in.seekg(currentpos + findpos + s.length());
      in >> x;
      in.seekg(nextlinepos);
    }
  }
  if (!foundit) found = false;
}

//----------------------------------------------------------------------------------
void ExperimentParameters::readExperimentParameters(vector<ExperimentParameters> & expPars, string file)
{

  bool allOK = true;
  expPars.clear();
  ExperimentParameters eP;
  ifstream infile(file.c_str());
  //infile.open(file.c_str(), ios::nocreate);
  int start = 0;
  while ( allOK )
  {
    readPar(eP.repetitions, infile, "repetitions", start, allOK, " = ");
    readPar(eP.gridWidth, infile, "gridWidth", start, allOK, " = ");
    readPar(eP.gridHeight, infile, "gridHeight", start, allOK, " = ");
    readPar(eP.cellwidth, infile, "cellwidth", start, allOK, " = ");
    readPar(eP.numyears, infile, "numyears", start, allOK, " = ");
    readPar(eP.timeStepsPerYear, infile, "timeStepsPerYear", start, allOK, " = ");
    readPar(eP.timestepLength, infile, "timestepLength", start, allOK, " = ");

    if (allOK)
    {
      readPar(eP.initOffset, infile, "initOffset", start, allOK, " = ");
      if (allOK == false)
      {
        eP.initOffset = 0;
        allOK = true;
      }
      readPar(eP.deltaOffset, infile, "initOffset", start, allOK, " = ");
      if (allOK == false)
      {
        eP.deltaOffset = ceil(eP.numyears * eP.timeStepsPerYear * timestepLength);
        allOK = true;
      }
    }
    readPar(eP.soilDepth, infile, "soilDepth", start, allOK, " = ");
    readPar(eP.thetaSat, infile, "thetaSat", start, allOK, " = ");
    readPar(eP.thetaCrit, infile, "thetaCrit", start, allOK, " = ");
    //      readPar(eP.psiE, infile, "psiE", start, allOK, " = ");
    //      readPar(eP.b, infile, "b", start, allOK, " = ");
    //      readPar(eP.airTemperature, infile, "airTemperature", start, allOK, " = ");
    //      readPar(eP.airHumidity, infile, "airHumidity", start, allOK, " = ");
    //      readPar(eP.evapotranspiration, infile, "evapotranspiration", start, allOK, " = ");
    readPar(eP.kappa, infile, "kappa", start, allOK, " = ");
    // Seedbank and seed dispersal parameter
    readPar(eP.mu, infile, "mu", start, allOK, " = ");
    readPar(eP.sd, infile, "sd", start, allOK, " = ");
    readPar(eP.seedrainIntensity, infile, "seedrainIntensity", start, allOK, " = ");
    readPar(eP.seedMortality, infile, "seedMortality", start, allOK, " = ");
    readPar(eP.maxIndividualDensity, infile, "maxIndividualDensity", start, allOK, " = ");
    readPar(eP.sensitivityAnalysis, infile, "sensitivityAnalysis", start, allOK, " = ");
    readPar(eP.initSoilWaterContent, infile, "initSoilWaterContent", start, allOK, " = ");
    //      readPar(eP.barriersFile, infile, "barriersFile", start, allOK, " = ");
    readPar(eP.redistributionFactorFile, infile, "redistributionFactorFile", start, allOK, " = ");
    readPar(eP.dailyRainFile, infile, "dailyRainFile", start, allOK, " = ");
    readPar(eP.plantInitFile, infile, "plantInitFile", start, allOK, " = ");
    readPar(eP.plantParsFile, infile, "plantParsFile", start, allOK, " = ");
    readPar(eP.metacommunityFile, infile, "metacommunityFile", start, allOK, " = ");
    readPar(eP.outputFilename, infile, "outputFilename", start, allOK, " = ");
    readPar(eP.outputOptions, infile, "outputOptions", start, allOK, " = ");
    readPar(eP.reportCycle, infile, "reportCycle", start, allOK, " = ");
    readPar(eP.dailyEvapoFile, infile, "dailyEvapoFile", start, allOK, " = ");
    readPar(eP.rocksAndPotsFile, infile, "rocksAndPotsFile", start, allOK, " = ");
    eP.uptakeSlope = 4.3;
    //      readPar(eP.etAM0, infile, "etAM0", start, allOK, " = ");
    //      readPar(eP.etAM1, infile, "etAM1", start, allOK, " = ");
    //      readPar(eP.etAM2, infile, "etAM2", start, allOK, " = ");
    //      readPar(eP.etPH1, infile, "etPH1", start, allOK, " = ");
    //      readPar(eP.etPH2, infile, "etPH2", start, allOK, " = ");

    if ( allOK ) {
      readPar(eP.simulationType, infile, "simulationType", start, allOK, " = ");
      if ( allOK == false ) {
        eP.simulationType = "simple";
        allOK = true;
      }
    }

    if ( allOK ) {
      readPar(eP.initSeedbankDensity, infile, "initSeedbankDensity", start, allOK, " = ");
      if ( allOK == false ) {
        eP.initSeedbankDensity = 0.0;
        allOK = true;
      }
    }


    if ( allOK ) {
      readPar(eP.startGrazingPeriod, infile, "startGrazingPeriod", start, allOK, " = ");
      if ( allOK == false ) {
        eP.startGrazingPeriod = 0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readPar(eP.endGrazingPeriod, infile, "endGrazingPeriod", start, allOK, " = ");
      if ( allOK == false ) {
        eP.endGrazingPeriod = 0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readPar(eP.grazingFrequency, infile, "grazingFrequency", start, allOK, " = ");
      if ( allOK == false ) {
        eP.grazingFrequency = 0.0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readPar(eP.grazingAnnualAmount, infile, "grazingAnnualAmount", start, allOK, " = ");
      if ( allOK == false ) {
        eP.grazingAnnualAmount = 0.0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readPar(eP.dewAmount, infile, "dewAmount", start, allOK, " = ");
      if ( allOK == false ) {
        eP.dewAmount = 0.0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readPar(eP.dewDays, infile, "dewDays", start, allOK, " = ");
      if ( allOK == false ) {
        eP.dewDays = 0.0;
        allOK = true;
      }
    }

    /*
     if ( allOK ) {
     readPar(eP.uptakeSlope, infile, "uptakeSlope", start, allOK, " = ");
     if ( allOK == false ) {
     eP.uptakeSlope = 4.3;
     allOK = true;
     }
     }

     */

    if ( allOK )
    {
      expPars.push_back(eP);
      start = infile.tellg();
    }
  }
  infile.close();
}

void ExperimentParameters::readExperimentParametersNew(vector<ExperimentParameters> & expPars, string file)
{

  // This function assumes that all parameters for a given run are written in one line.

  bool allOK = true;
  expPars.clear();
  ExperimentParameters eP;
  ifstream infile(file.c_str());
  //infile.open(file.c_str(), ios::nocreate);

  const int stringLength = 4096;

  char str[stringLength];          // necessary, because borland has not yet implemented getline for stl strings
  string line;

  while (infile.peek() != EOF )
  {
    infile.getline(str, stringLength); // necessary, because borland has not yet implemented getline for stl strings
    line = str;            // necessary, because borland has not yet implemented getline for stl strings

    allOK = true;
    readParNew(eP.repetitions, line, "repetitions", allOK, " = ");
    readParNew(eP.gridWidth, line, "gridWidth", allOK, " = ");
    readParNew(eP.gridHeight, line, "gridHeight", allOK, " = ");
    readParNew(eP.cellwidth, line, "cellwidth", allOK, " = ");
    readParNew(eP.numyears, line, "numyears", allOK, " = ");
    readParNew(eP.timeStepsPerYear, line, "timeStepsPerYear", allOK, " = ");
    readParNew(eP.timestepLength, line, "timestepLength", allOK, " = ");

    if (allOK)
    {
      readParNew(eP.initOffset, line, "initOffset", allOK, " = ");
      if (allOK == false)
      {
        eP.initOffset = 0;
        allOK = true;
      }
      readParNew(eP.deltaOffset, line, "deltaOffset", allOK, " = ");
      if (allOK == false)
      {
        eP.deltaOffset = ceil(eP.numyears * eP.timeStepsPerYear * eP.timestepLength);
        allOK = true;
      }
    }
    readParNew(eP.soilDepth, line, "soilDepth", allOK, " = ");
    readParNew(eP.thetaSat, line, "thetaSat", allOK, " = ");
    readParNew(eP.thetaCrit, line, "thetaCrit", allOK, " = ");
    readParNew(eP.kappa, line, "kappa", allOK, " = ");

    // Seedbank and seed dispersal parameter
    readParNew(eP.mu, line, "mu", allOK, " = ");
    readParNew(eP.sd, line, "sd", allOK, " = ");
    readParNew(eP.seedrainIntensity, line, "seedrainIntensity", allOK, " = ");
    if (allOK)
    {
      readParNew(eP.seedrainYears, line, "seedrainYears", allOK, " = ");
      if (allOK == false)
      {
        eP.seedrainYears = 100;
        allOK = true;
      }
    }
    readParNew(eP.seedMortality, line, "seedMortality", allOK, " = ");
    readParNew(eP.maxIndividualDensity, line, "maxIndividualDensity", allOK, " = ");
    readParNew(eP.sensitivityAnalysis, line, "sensitivityAnalysis", allOK, " = ");
    readParNew(eP.initSoilWaterContent, line, "initSoilWaterContent", allOK, " = ");
    readParNew(eP.redistributionFactorFile, line, "redistributionFactorFile", allOK, " = ");
    readParNew(eP.dailyRainFile, line, "dailyRainFile", allOK, " = ");
    readParNew(eP.plantInitFile, line, "plantInitFile", allOK, " = ");
    readParNew(eP.plantParsFile, line, "plantParsFile", allOK, " = ");
    readParNew(eP.metacommunityFile, line, "metacommunityFile", allOK, " = ");
    readParNew(eP.outputFilename, line, "outputFilename", allOK, " = ");
    readParNew(eP.outputOptions, line, "outputOptions", allOK, " = ");
    readParNew(eP.reportCycle, line, "reportCycle", allOK, " = ");
    readParNew(eP.dailyEvapoFile, line, "dailyEvapoFile", allOK, " = ");
    readParNew(eP.rocksAndPotsFile, line, "rocksAndPotsFile", allOK, " = ");

    if ( allOK ) {
      readParNew(eP.simulationType, line, "simulationType", allOK, " = ");
      if ( allOK == false ) {
        eP.simulationType = "simple";
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.initSeedbankDensity, line, "initSeedbankDensity", allOK, " = ");
      if ( allOK == false ) {
        eP.initSeedbankDensity = 0.0;
        allOK = true;
      }
    }


    if ( allOK ) {
      readParNew(eP.startGrazingPeriod, line, "startGrazingPeriod", allOK, " = ");
      if ( allOK == false ) {
        eP.startGrazingPeriod = 0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.endGrazingPeriod, line, "endGrazingPeriod", allOK, " = ");
      if ( allOK == false ) {
        eP.endGrazingPeriod = 0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.grazingFrequency, line, "grazingFrequency", allOK, " = ");
      if ( allOK == false ) {
        eP.grazingFrequency = 0.0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.grazingAnnualAmount, line, "grazingAnnualAmount", allOK, " = ");
      if ( allOK == false ) {
        eP.grazingAnnualAmount = 0.0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.dewAmount, line, "dewAmount", allOK, " = ");
      if ( allOK == false ) {
        eP.dewAmount = 0.0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.dewDays, line, "dewDays", allOK, " = ");
      if ( allOK == false ) {
        eP.dewDays = 0.0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.uptakeSlope, line, "uptakeSlope", allOK, " = ");
      if ( allOK == false ) {
        eP.uptakeSlope = 4.3;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.globalDispersal, line, "globalDispersal", allOK, " = ");
      if ( allOK == false ) {
        eP.globalDispersal = false;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.samplingDate, line, "samplingDate", allOK, " = ");
      if ( allOK == false ) {
        eP.samplingDate = 0;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.startEstablishmentReporting, line, "startEstablishmentReporting", allOK, " = ");
      if ( allOK == false ) {
        eP.startEstablishmentReporting = eP.numyears + 1;
        allOK = true;
      }
    }

    if ( allOK ) {
      readParNew(eP.stopEstablishmentReporting, line, "stopEstablishmentReporting", allOK, " = ");
      if ( allOK == false ) {
        eP.stopEstablishmentReporting = eP.numyears + 1;
        allOK = true;
      }
    }

    if ( allOK )
    {
      expPars.push_back(eP);
    }
  }
  infile.close();
}


template<class T>
void readParNew(T& x, const string & line, string key, bool& found, string sep)
{
  //  T x;
  bool foundit = false;
  string s = key + sep;
  size_t findpos = line.find(s, 0);
  if (findpos != string::npos) //i.e. string s was found in line
  {
    foundit = true;
    istringstream is(line);
    is.seekg(findpos + s.length());
    is >> x;
  }
  if (!foundit) found = false;
}

void Simulation::reportBiomass(Simulation* testSimu) {
  const list<plant> & plants = testSimu->getPlantList();
  const vector<plantParams> & pP = testSimu->getPlantPars();
  list<plant>::const_iterator itPlants;
  valarray<double> biomass(0.0, pP.size());
  itPlants = plants.begin();
  while (itPlants != plants.end()) {
    biomass[(*itPlants).getspeciesid()] += (*itPlants).getCohortSize() * (*itPlants).getmass();
    itPlants++;
  }

// A.save("A.txt", arma_ascii);

 for (size_t i = 0; i < biomass.size()-1; i++)
 {
 testSimu->gOutBiomass << biomass[i] << " ";
  speciesBiomass(testSimu->year, i) =  biomass[i];
 }
 testSimu->gOutBiomass << biomass[biomass.size()-1] << "\n";
 speciesBiomass(testSimu->year, biomass.size()-1) = biomass[biomass.size()-1];

}


void Simulation::reportGrowth(Simulation* testSimu) {
  const list<plant> & plants = testSimu->getPlantList();
  const vector<plantParams> & pP = testSimu->getPlantPars();
  list<plant>::const_iterator itPlants;
  itPlants = plants.begin();
  while (itPlants != plants.end()) {
    double gg = itPlants->getGrowth();
    double deltaGrowth = (*itPlants).getCohortSize() * gg;
    speciesGrowth(testSimu->year, (*itPlants).getspeciesid()) += deltaGrowth;
    itPlants++;
  }

}




void Simulation::initOutBiomass(std::string outputFilename, int j) {
    std::string ofn = outputFilename + "_biomass_" + convert<std::string>(j) + ".dat";
    gOutBiomass.open(ofn.c_str(), ios::app);
    gOutBiomass.setf(ios::fixed, ios::floatfield);
    gOutBiomass.precision(3);
}
/**/
