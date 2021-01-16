//---------------------------------------------------------------------------
#pragma once

#ifndef UPlantH
#define UPlantH
//---------------------------------------------------------------------------
#endif


#include <valarray>
#include <vector>
#include <list>
#include <string>
#include <math.h>
//#include <iostream.h>
//#include <fstream.h>
#include <fstream>
// #include <gsl/gsl_randist.h>
// #include <gsl/gsl_rng.h>
#include <cassert>
#include <Rcpp.h>

using namespace std;

// extern gsl_rng * theRNG ;
// extern ofstream gOutEstablishmentDistances;
double getSquareDistanceNew(unsigned int cellID1, unsigned int cellID2, unsigned int nrow, unsigned int ncol);

struct plantParams
{
    double pRoot;                          // Allocation to Roots (fraction)
    double pStorage;                       // Alloation to Storage (fraction)
    double pLeaf;                          // Allocation to Leaves (fraction)
    double maxSize;                        // Maximum size (g Carbon)
    double TranspirationSlope;             // Slope of regulation of Transpiration in response to soil water content (m3 m-3)
    double HalfSaturationTheta;            // Soil water content at which potential Transpiration is halfway between minimum and maximum potential Transpiration (m3 m-3)
//    double minSpecificTranspirationRate;   // g H20 g-1 Leaf Carbon day-1
    double maxSpecificTranspirationRate;   // g H20 g-1 Leaf Carbon day-1
    double maxUptakeRate;                  // g H20 g-1 Root Carbon day-1
    double respirationRate;                // g C g-1 C day-1
    double maxRootDensity;                 // g Root Carbon per dm2
    double waterUseEfficiency;             // g Carbon gained (gross!) per g H20 transpired
    double succulenceFactor;               // g H20 stored per g Storage Carbon
    double seedlingsize;                       // g Carbon per seedling
    int   maxAge;                         // year
    double maxAgeGamma;
    double germinationProbability;          // Probability of seed to germinate year-1
    double desiccationAlpha;
    double desiccationGamma;
    double rGrowth;
    double baseMortalityRate;
    double seedsPerSeedmass;
//    double psiCrit;
    double maxLeafDensity;
    double leafStructureAllometryPower;
    double leafStructureAllometryFactor;
    double LAIopt;

    unsigned int dispersalDate;
    unsigned int germinationDate;
    unsigned int rainPeriod;
    double rainfallThreshold;
    double uptakeSlope;
};


/* TODO : put this in BMisc.h */
/*
template<class T>
void writeErrorMessage(string file, string message, T value)
{
    ofstream out;
    out.open(file.c_str(), ios::app);
    out << "Error: " << message << value << "\n";
    out.close();
};
*/

struct NNcell{
    long int cellID;
    bool accessible;
    bool free;
    double distance;
};

void getCandidates(vector<unsigned int> & candidates, const valarray<bool> & selected, unsigned int cellID, unsigned int nrow, unsigned int ncol, const vector<NNcell> & barriers);
unsigned int getNearestCandidate(vector<unsigned int> & candidates, unsigned int cellID, unsigned int nrow, unsigned int ncol);
double getSquareDistance(unsigned int cellID1, unsigned int cellID2, unsigned int nrow, unsigned int ncol);

// warum ist diese Structur nicht auch private f?r die Klasse plant?
  struct localLAI
  {
    unsigned int cellID;
    double LAI;
  } ;


//----------------------------------------------------------------------------------
//!  The plant class.
/*!
  An instant of the plant class is a cohort, the size of which can be one.
*/

class plant
{
  private:
  struct localRootMass
  {
    unsigned int cellID;
    double mass;
    double waterDemand;
  } ;

   double timeToMaturity;
   double mass, rootmass, storagemass, leafmass, seedmass;
   double waterDemand;
   double waterUptake;
   double trueWaterUptake;
   double potentialSpecificTranspirationRate;
   int speciesID;
   double birthday;
   int recruitmentYear;
   vector<localLAI> localLeafAreaIndexVector;         // wo wird der initialisiert? In initNewCohort()

   double growth;
   double plantWaterContent;
   double maxCurrentUptakeRate;
//   double maxRootMass;
   double lastTimePositiveGrowth;
   unsigned int kohortSize;
   bool desiccating;
   double Amin;
   unsigned int maxRootCells;
   unsigned int maxLeafCells;
   unsigned int currentRootCells;
   vector<localRootMass> rootMassVector;
   unsigned int cellID;
   size_t maxRootCellsIndex;
   size_t maxLeafCellsIndex;
   unsigned int cohortID;
   double mortalityprob;

/*
   unsigned int numberOfRootCells;
   double rootmassInLastCell;
   double rootmassInFirstCells;
*/

  public:
//    vector<localRootMass> localRootMassesVector;    // wieso ist der public?
//  void initPositionVector(const vector<plantParams>& pP, unsigned int nrow, unsigned int ncol, double cellArea, bool root, const vector<NNcell> & barriers); /*!< Initialises the position vector */
  double getmass() const {return mass;}; /*!< Returns plant mass */
  int getspeciesid() const {return speciesID;}; /*!< Returns species ID */
  int getSpeciesID() {return speciesID;}; /*!< Returns species ID */
//  const int getspeciesid() {return speciesID;}; /*!< Returns species ID */

  double getbirthday() {return birthday;};
  double getwaterdemand() {return waterDemand;};
  double getwateruptake() {return waterUptake;};
  double getTrueWaterUptake() {return trueWaterUptake;};
  double gettranspiration() {return potentialSpecificTranspirationRate;};
  double getseedmass() {return seedmass;};
  const vector<localLAI> & getLAIvector() {return localLeafAreaIndexVector;};
  double getCohortSize() const {return kohortSize;};
  double getPotentialSpecificTranspirationRate() {return potentialSpecificTranspirationRate;};
  double getAmin(){return Amin;};
  double getWaterDemand(){return waterDemand;};
  size_t getMaxLeafCellsIndex(){return maxLeafCellsIndex;};
  double getGrowth() const {return growth;};
  double getMortalityprob(){return mortalityprob;};
  double getStoragemass() {return storagemass;};
  double getLeafmass() {return leafmass;};
  double getPlantWaterContent() {return plantWaterContent;};
  unsigned int getCellID() {return cellID;};
  int getRecruitmentYear() {return recruitmentYear;};
  double getTimeToMaturity() {return timeToMaturity;};

  void setPlantWaterContent(double x) {plantWaterContent = x;};
  void setCohortID(unsigned int x) {cohortID = x;};

  //! A public function.
  /*!
    Sets the water uptake by the plant. Used in plant initialization.
  */
  void setTimeToMaturity(double x) {timeToMaturity = x;};
  void setWaterUptake(double x) {waterUptake = x;};
  void setTrueWaterUptake(double x) {trueWaterUptake = x;};
  void setbirthday(int x) {birthday = x;};       // soll weg
  void setBirthday(double x) {birthday = x;};
  void setspeciesid(int x) {speciesID = x;};     // Schreibweise!
  void setmass(double m) {mass = m;};
  void setSeedmass(double m) {seedmass = m;};
  void setLastTimePositiveGrowth(double x) {lastTimePositiveGrowth = x;};
  void setDesiccating(bool d) {desiccating = d;};
  void setCellID(int x, int y, int ncol)
  {
    cellID = ncol*y + x;
  };
  void setCellID(unsigned int x) {cellID = x;};
  void setCohortSize(double x){kohortSize = x;};
  void setGrowth(double x){growth = x;};
  void setAmin(const vector<plantParams>& pP);
  void setLeafmass(double x){leafmass = x;};
  void setStoragemass(double x){storagemass = x;};
  void setRecruitmentYear(int x){recruitmentYear = x;};

  void reportDistancesToPotentialMotherPlants( const valarray<double> & newSeeds );



//------------------------------------------------------------------------------

//! Calculates the transpiration rate
/* The transpiration rate is a function of plant water storage and soil water content. Water has to pass through the
storage before it can be transpired through the leaves.

If the plant has no water in the storage, then no transpiration occurs.
\f[ A_i = \left\{
\begin{array}{r@{,\quad}l}
                  0 & \mbox{if } W_i = 0 \\
                  A_{min,j} + \frac{W_i}{W_{max,i}} \frac{A_{max} -
A_{min,j}}{1 + e^{- \beta_i (\bar{\theta_i} - \theta_{M,i})}} &
\mbox{otherwise}.
                  \end{array} \right.
\f]
*/

//  void potentialTranspiration( const vector<plantParams>& pP, const valarray<double>& theta, double ET);
  void potentialTranspirationNew( const vector<plantParams>& pP, const valarray<double>& theta, double ET, double cellArea, const vector< vector<unsigned int> > & rootPositions);
  void waterdemand2( const vector<plantParams>& pP, double deltatime);

  //------------------------------------------------------------------------------
    //! The maximum amount of water the cohort can take up during the next time step.
    /*! The total amount of water that the plant cohort could maximally take up during a time step
       is ordered from the different cells where the plant cohort has roots. The weight for the different
       cells is given by the fraction of the cohort's root mass that is in those cells.
    */

//  void waterorder(valarray<double>& plantwaterdemand,  const vector<plantParams>& pP, const valarray<double>& theta, double soilDepth, double timestepLength, double thetaCrit, double ET);
//  void wateruptake2(const valarray<double>& plantWaterDemand, const valarray<double>& realizedPlantWaterDemand);
  void waterorderNew(valarray<double>& plantwaterdemand,  const vector<plantParams>& pP, const valarray<double>& maxWaterUptakeRate, double soilDepth, double timestepLength, double thetaCrit, double thetaSat, double ET, double cellArea, const vector< vector<unsigned int> > & rootPositions);
  void wateruptakeNew(const valarray<double>& plantWaterDemand, const vector<plantParams>& pP, const valarray<double>& realizedPlantWaterDemand, double cellArea, const vector< vector<unsigned int> > & rootPositions, const valarray<double>& maxWaterUptakeRate);

//------------------------------------------------------------------------------

//! Production of photosynthates
/*!
Production is the difference between gross
photosynthesis and respiration. Gross photosynthesis is the
product of water use efficiency and transpiration, and respiration
is directly proportional to biomass: \f[
P_i = \mbox{\textsc{wue}}  \cdot T_i - r_{resp} M_i
\f]
It is assumed that the average respiration rate, \f$ r_{resp} \f$, is
independent of the plant's allocation strategy.


\f[
P_i = \left\{ \begin{array}{r@{,\quad}l}
   (1 - r_{growth}) \cdot ( \mbox{\textsc{wue}} \cdot A_i \cdot p_{L,j} - r_{resp}) \cdot
   M_i & \mbox{if } \mbox{\textsc{wue}} \cdot A_i \cdot p_{L,j} > r_{resp}\\
   ( \mbox{\textsc{wue}} \cdot A_i \cdot p_{L,j} - r_{resp}) \cdot
   M_i    & \mbox{otherwise}
         \end{array} \right.
\f]
*/

  void production3 (  const vector<plantParams>& pP, double deltatime, double time  );
  bool mortality( const vector<plantParams>& pP, valarray<unsigned int> & plantDensity, double time, double deltatime, const string & simulationType, valarray<double> & totalDeaths, valarray<double> & droughtDeaths);

  //! Allocates the photosynthetic products
//  void allocation2( const vector<plantParams>& pP, double cellArea, const vector<double>& leafFraction);
//  void initAllocation( const vector<plantParams>& pP, double cellArea);
  void initAllocationNew( const vector<plantParams>& pP, double cellArea, const vector< vector<unsigned int> > & rootPositions);
  void allocationNew( const vector<plantParams>& pP, double cellArea, const vector<double>& leafFraction, const vector< vector<unsigned int> > & rootPositions, int year, double simulatedTime, double cellwidth, unsigned int width, unsigned int height, const list< valarray<unsigned int> > & newSeedsOfPastYears, const vector<unsigned int> & dispersalDates, bool reportEstablishment);


  void seedproduction( const vector<plantParams>& pP, valarray<unsigned int>& newSeeds, unsigned int arenaSize);
  void updateLeaves2(const vector<plantParams>& pP, double cellArea);
  void updateLeavesNew(const vector<plantParams>& pP, double cellArea);
  void updateRootsNew();
/*
  void writePlant(ofstream * pout, double time) const;
  void writePlantBin(ofstream * pout, double time) const;
  void writePlantShort(ofstream * pout, double time) const;
 */
};

