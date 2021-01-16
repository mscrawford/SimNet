//---------------------------------------------------------------------------

#ifndef USimuH
#define USimuH
//---------------------------------------------------------------------------
#endif

#include <list>
#include <fstream>
// #include <map>
#include "UPlant.h"
//#include "UStringUtils.h"


// sollte in einen eigenen header
double distance(double x, double y, double centerX, double centerY, unsigned int gridWidth, unsigned int gridHeight);

// sollte in einen eigenen header
template<class T>
void readPar(T& x, ifstream & in, string key, int start, bool& found, string sep);

template<class T>
void readParNew(T& x, const string & line, string key, bool& found, string sep);

//template<class T>
//void readValarray(valarray<T>& x, string file, string key);
void readValarray(valarray<double>& x, string file, string key);
void readVector(vector<int>& x, string file, string key);


//! Plant initialization parameters
/*!
    - species ID
    - position (row, col)
    - cohort size (number of individuals)
    - birthday  (must be <= 0.0)
    - mass (g)
*/
struct plantInitParams2
{
     int speciesID;
     int row;
     int col;
     unsigned int cohortSize;
     double birthday;
     double mass;
};

//! Contains all parameters necessary for an experiment      -- sollte in einen eigenen Header
class ExperimentParameters
{
  public:
  unsigned int repetitions;
  unsigned int initOffset;
  unsigned int deltaOffset;
  unsigned int gridWidth, gridHeight;  // initialized in setParametersGUI
  double cellwidth;                    // initialized in setParametersGUI
//  double cellArea;                     // initialized in setParametersGUI

//  double simulationDuration;           // initialized in setParametersGUI
//  double timestepLength;               // initialized in setParametersGUI
  unsigned int numyears;                        // initialized in setParametersGUI
  unsigned int timeStepsPerYear;                // initialized in setParametersGUI
  double timestepLength;

  // Soil parameters
  double soilDepth;                    // initialized in setParametersGUI
  double thetaSat;                     // initialized in setParametersGUI
  double thetaCrit;
  double uptakeSlope;
  //  double psiE;                         // initialized in setParametersGUI
//  double b;                            // initialized in setParametersGUI

  // Evaporation parameters
//  double airTemperature;               // initialized in setParametersGUI
//  double airHumidity;                  // initialized in setParametersGUI
//  double evapotranspiration;           // initialized in setParametersGUI
  double kappa;                        // initialized in setParametersGUI

//  double etAM0, etAM1, etAM2, etPH1, etPH2;

  // Seedbank and seed dispersal parameter
  double mu;                           // initialized in setParametersGUI
  double sd;                           // initialized in setParametersGUI
  double seedrainIntensity;            // initialized in setParametersGUI
  unsigned int seedrainYears;
  double seedMortality;                // initialized in setParametersGUI
  double initSeedbankDensity;

  double maxIndividualDensity;         // initialized in setParametersGUI

  bool sensitivityAnalysis;            // Determines how the plantParsFile is
                                       // interpreted. If sensitivityAnalysis is
                                       // true, then one set of
  bool globalDispersal;

  double initSoilWaterContent;

  unsigned int startGrazingPeriod;
  unsigned int endGrazingPeriod;
  double grazingFrequency;
  double grazingAnnualAmount;

  double dewAmount;
  double dewDays;

  string barriersFile;
  string redistributionFactorFile;
  string dailyRainFile;
  string plantInitFile;
  string plantParsFile;
  string metacommunityFile;
  string dailyEvapoFile;
  string simulationType;
  string rocksAndPotsFile;

  string outputFilename;
  string outputOptions;
  unsigned int reportCycle;
  unsigned int samplingDate;
  unsigned int startEstablishmentReporting;
  unsigned int stopEstablishmentReporting;

//  double minLambda = 0.00001;
//  double waterDensity = 1.0;

  //! Reads experiment parameters from file; stores them in a vector of ExperimentParameters
  void readExperimentParameters(vector<ExperimentParameters> & expPars, string file);
  void readExperimentParametersNew(vector<ExperimentParameters> & expPars, string file);

};



//! A simulation
/*!
   Contains the information on the spatio-temporal dimensions of the system, the environmental parameters (soil, evaporation, seedbank, maximum plant density).
*/
class Simulation
{
  private:

  long double cumulativePlantWaterUptake;
  long double cumulativeDrainage;
  long double cumulativeEvaporation;

  unsigned int startEstablishmentReporting;
  unsigned int stopEstablishmentReporting;


  ExperimentParameters experimentPars;
  unsigned int repetition;
  // System dimensions in space and time
  unsigned int gridWidth, gridHeight;  // initialized in setParametersGUI
  double simulationDuration;           // initialized in setParametersGUI
  double timestepLength;               // initialized in setParametersGUI
  int numyears;                        // initialized in setParametersGUI
  int timeStepsPerYear;                // initialized in setParametersGUI
  double cellwidth;                    // initialized in setParametersGUI
  double cellArea;                     // initialized in setParametersGUI
  double simulatedTime;                // initialized in Simulation::run2

  // Soil parameters
  double soilDepth;                    // initialized in setParametersGUI
  double thetaSat;                     // initialized in setParametersGUI
//  double psiE;                         // initialized in setParametersGUI
//  double b;                            // initialized in setParametersGUI
  double thetaCrit;
  double uptakeSlope;

  // Evaporation parameters
//  double airTemperature;               // initialized in setParametersGUI
//  double airHumidity;                  // initialized in setParametersGUI
  double evapotranspiration;           // initialized in setParametersGUI
  vector<double> evapotranspirationVector;
  double kappa;                        // initialized in setParametersGUI

  // Seedbank and seed dispersal parameter
  double mu;                           // initialized in setParametersGUI
  double sd;                           // initialized in setParametersGUI
  double seedrainIntensity;            // initialized in setParametersGUI
  unsigned int seedrainYears;
  double seedMortality;                // initialized in setParametersGUI
  valarray<double> metacommunity;      // initialized in setParametersGUI
  double initSeedbankDensity;


  double maxIndividualDensity;         // initialized in setParametersGUI

  double drainageAmount;                     // initialized and updated in Simulation::run2
  double rain;                         // initialized and updated in Simulation::run2

  // initialization Values
  double initSoilWaterContent;         // initialized in setParametersGUI

  double minLambda;                          // initialized in setParametersGUI
  double waterDensity;                       // initialized in setParametersGUI

  list<plant> allPlants;                     // initialized in setParametersGUI- via function call initPlants - not checked
  list<plant>::iterator itAllPlants;         // initialized and updated in Simulation::run2
  vector<plantParams> plantPars;             // initialized in setParametersGUI- via function call readPlantPars   - not checked
  vector<plantInitParams2> plantInitPars;    // initialized in setParametersGUI - via function call readPlantInitFile   - not checked

  // external Initialization
  valarray<unsigned int> plantDensity;       // initialized in Simulation::initPlants - repetition update in Simulation::initRepetition via function call Simulation::initPlants -- this way the valarray is newly created!!!
  vector<double> dailyRainfall;              // initialized in setParametersGUI - via function call readDailyRain
  valarray<double> redistributionFactor;     // initialized in setParametersGUI
  valarray<double> waterContent;             // initialized in Simulation::init, reset to initWaterContent in Simulation::initRepetition

  // internal Initialization
  valarray<double> potentialEvaporation;  // initialized in Simulation::init
  valarray<double> actualEvaporation;     // initialized in Simulation::init
  valarray<double> leafAreaIndex;         // initialized in Simulation::init
  valarray<double> dispersalKernel;       // initialized in Simulation::init via function call initDispersalKernel
  valarray<double> plantWaterDemand;      // g H2O, initialized in Simulation::init() - reset to zero in Simulation::initRepetition
  valarray<double> realizedPlantWaterDemand;      // g H2O, initialized in Simulation::init() - reset to zero in Simulation::initRepetition
  valarray<double> dewfallSequence;      //
  valarray<double> maxWaterUptakeRate;             // initialized in Simulation::init, reset to initWaterContent in Simulation::initRepetition
  valarray<double> totalDeaths;
  valarray<double> droughtDeaths;

  valarray<unsigned int> seedbank;        // initialized in Simulation::init() - reset to zero in Simulation::initRepetition
  valarray<unsigned int> newSeeds;        // initialized in Simulation::init() - reset to zero in Simulation::initRepetition
//  valarray<double> seedbankNew;        // initialized in Simulation::init() - reset to zero in Simulation::initRepetition
  list< valarray<unsigned int> > newSeedsOfPastYears;

  vector<NNcell> barriers;                // initialized in setParametersGUI - not correctly though!
  int year;                               // initialized in Simulation::run2  - and updated there (for loop)
  int timestep;                           // initialized in Simulation::run2  - and updated there (for loop)
  int offset;                    // initialized and updated in Simulation::initRepetition

  vector<unsigned int> germinationDates;
  vector<unsigned int> dispersalDates;
  vector<double> leafFraction;
//  vector<double> deltaLeafFraction;

  vector<int> rocksAndPots;
  vector<unsigned int> masterLeafIndexVector;

//  map<string, ofstream *> outputFileMap;
//  ofstream out;

  vector<long double> growth;
  vector<long double> survivalProbability;

  vector< vector<unsigned int> > rootPositions;
  vector< vector<unsigned int> > leafPositions;


  string simulationType;                  // initialized in Simulation::init

  unsigned int startGrazingPeriod;
  unsigned int endGrazingPeriod;
  double grazingFrequency;
  double grazingAnnualAmount;

  double dewAmount;
  double dewDays;

  bool globalDispersal;

  public:

  unsigned int cohortCounter;

  const list<plant> & getPlantList() {return allPlants;};
  const valarray<double> & getDispersalKernel() {return dispersalKernel;};
  const valarray<unsigned int> & getSeedbank() {return seedbank;};
  const vector<double> & getDailyRainfall() {return dailyRainfall;};
  const valarray<double> & getSoilWaterContent() {return waterContent;};
  const valarray<double> & getPlantWaterDemand() {return plantWaterDemand;};
  const valarray<double> & getActualEvaporation() {return actualEvaporation;};
  const valarray<double> & getRealizedPlantWaterDemand() {return realizedPlantWaterDemand;};
  const valarray<double> & getTotalDeaths() {return totalDeaths;};
  const valarray<double> & getDroughtDeaths() {return droughtDeaths;};
  const vector<plantParams> & getPlantPars() {return plantPars;};
  double getSimulatedTime() {return simulatedTime;};
//  const map<string, ofstream *> & getOutputFileMap() {return outputFileMap;};
  const string & getSimulationType() {return simulationType;};
  const vector<int> & getRocksAndPots() {return rocksAndPots;};

  int getCurrentYear () {return year;};
  int getCurrentTimestep () {return timestep;};
  int getTimeStepsPerYear() {return timeStepsPerYear;};
  unsigned int getRepetition() {return repetition;};
  double getRain() {return rain;};
  double getDrainageAmount() {return drainageAmount;};
  long double getCumulativePlantWaterUptake() {return cumulativePlantWaterUptake;};
  long double getCumulativeDrainage() {return cumulativeDrainage;};
  long double getCumulativeEvaporation() {return cumulativeEvaporation;};


//  void init();
//  void initMasterLeafIndexVector();
  void init(const ExperimentParameters & eP);
  void initRNG(unsigned int rngSeed);
//  void initRepetition(bool first);
  void initRepetition();
  void initPlants();
  void run();
  void setParameters(); //(int width, int height, double duration, double timestep, double cellsize);
  void setParametersGUI();
  void initDispersalKernel(double mu, double sd);
//  void evaporation();
  void evaporationRates();
  void evaporationAndPlantWaterUptake();
  void evaporationAndPlantWaterUptakeNew();
  void germination();
  void externalSeedrain();
  void seedbankmortality();
  void dispersal();

  void grazing(const vector<plantParams> & pP);

//  void updateLAI();
  void updateLAInew(const vector<plantParams>& pP);
  void readDailyRain(string file);
  void readPlantInitFile(string file);

//  void plantWaterUptake();
//  void plantWaterUptake2();
  void rainfall();
  void dewfall();
  void drainage();

  //! Reads plant parameters and initializes vector<plantParams> plantPars
  void readPlantPars(string plantParFile);
  void initNewCohort(unsigned int speciesID, unsigned int cohortSize, const vector<plantParams>& pP, list<plant>& allPlants, unsigned int nrow, unsigned int ncol, double cellArea, const vector<NNcell> & barriers, double simulatedTime,unsigned int cellID );

  void openOutputFiles();
  void closeOutputFiles();

  void initEvapotranspiration(const ExperimentParameters & eP);
  void evaluateSimulation(const ExperimentParameters & eP);

  void initLeafFraction(const vector<plantParams> & pP);
  void initRootPositions(const vector<plantParams> & pP,  unsigned int nrow, unsigned int ncol, const vector<int> & rocksAndPots);
  void initLeafPositions(const vector<plantParams> & pP,  unsigned int nrow, unsigned int ncol, const vector<int> & rocksAndPots);
  void initPositionVectorNew(vector< vector<unsigned int> > & positionVector, unsigned int nrow, unsigned int ncol, unsigned int maxNumCells, const vector<int> & rocksAndPots); /*!< Initialises the position vector */

  void adjustRedistributionMap();
  void initGerminationDates(const vector<plantParams> & pP);
  void initDispersalDates(const vector<plantParams> & pP);
  void reportBiomass(Simulation* testSimu);
  void reportGrowth(Simulation* testSimu);
  std::ofstream gOutBiomass;
  void initOutBiomass(std::string outputFilename, int j);
  Rcpp::NumericMatrix speciesBiomass;
  Rcpp::NumericMatrix speciesGrowth;
  Rcpp::IntegerMatrix speciesSeedproduction;

};
//int globalCounter;


// die sollten auch in einen eigenen Header
void plantscheck();
size_t getTargetIndex(size_t newCenterPoint, size_t toBeTransformed,  int ncol,  int nrow); // hier ist noch ein bug drin. da wird das Zentrum des Dispersal kernels verschoben.
size_t getTargetIndexDD(size_t newCenterPoint, double drow, double dcol,  int ncol,  int nrow);
// size_t getTargetIndexDD2(size_t newCenterPoint, double drow, double dcol,  int ncol,  int nrow);

//void initGraphics(int width, int height, int numyears);
void initColors(const vector<plantParams>& pP, bool expand);
void initPopSizeSeries();
void clearPopSizeSeries();
double f1 (double x, double M, double a, double b);
double optimalLeaffraction (double M, double a, double b);
double deriv1 (double x, double M, double a, double b);

void getCandidatesNew(vector<unsigned int> & candidates, const valarray<bool> & selected, unsigned int cellID, unsigned int nrow, unsigned int ncol, const vector<int> & barriers);
unsigned int getNearestCandidateNew(vector<unsigned int> & candidates, unsigned int cellID, unsigned int nrow, unsigned int ncol);
//double getSquareDistanceNew(unsigned int cellID1, unsigned int cellID2, unsigned int nrow, unsigned int ncol);
size_t getGerminationDateNew(const vector<double> & rain, int year, unsigned int germinationDate, unsigned int rainPeriod, double rainfallThreshold);


//! Calls the observer functions (graphics and file output). Called in Simulation::run2 at the beginning of every timestep.
void checkObserver(int simuTimestep, Simulation* simu);

