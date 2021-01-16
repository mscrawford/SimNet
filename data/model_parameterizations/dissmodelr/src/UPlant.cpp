#include "UPlant.h"

double getSquareDistance(unsigned int cellID1, unsigned int cellID2, unsigned int nrow, unsigned int ncol)
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

//--------------------------------------------------------------------


void getCandidates(vector<unsigned int> & candidates, const valarray<bool> & selected, unsigned int cellID, unsigned int nrow, unsigned int ncol, const vector<NNcell> & barriers)
{
   vector<unsigned int> fourNeighbors(4);
   unsigned int up, le, ri, lo;
   bool isAccessible(unsigned int target, unsigned int source, const vector<NNcell> & barriers);

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

   for (size_t i = 0; i < fourNeighbors.size(); i++)
   {
     if (selected[fourNeighbors[i]] == false && isAccessible(fourNeighbors[i], cellID, barriers))
     {
       candidates.push_back(fourNeighbors[i]);
     }
   }

};

//--------------------------------------------------------------------

   bool isAccessible(unsigned int target, unsigned int source, const vector<NNcell> & barriers)
   {
     return(true);
   }

//--------------------------------------------------------------------


unsigned int getNearestCandidate(vector<unsigned int> & candidates, unsigned int cellID, unsigned int nrow, unsigned int ncol)
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

//------------------------------------------------------------------------------

  void plant::updateLeaves2(const vector<plantParams>& pP, double cellArea)
  {
    plantParams lpp = pP[speciesID];
    maxLeafCellsIndex = ceil(leafmass / (cellArea * lpp.maxLeafDensity));
    if (maxLeafCellsIndex - 1 < localLeafAreaIndexVector.size())
    {
       for (size_t i = 0; i < maxLeafCellsIndex - 1; i++) localLeafAreaIndexVector[i].LAI = lpp.maxLeafDensity;
       localLeafAreaIndexVector[maxLeafCellsIndex - 1].LAI = (leafmass - lpp.maxLeafDensity * (maxLeafCellsIndex - 1) * cellArea)/cellArea;
    }
    else
       for (size_t i = 0; i < localLeafAreaIndexVector.size(); i++) localLeafAreaIndexVector[i].LAI = lpp.maxLeafDensity;

  }

//------------------------------------------------------------------------------
/*
void plant::writePlant(ofstream * pout, double time) const
{
   (* pout) << time << " " << cohortID << " " << speciesID << " " << cellID  << " " << kohortSize << " " << mass << " "
            << birthday << " " <<  leafmass  << " " << desiccating << " " << plantWaterContent << " "
            << potentialSpecificTranspirationRate << " " << maxCurrentUptakeRate << " " << lastTimePositiveGrowth << " "
            << waterDemand << " " << waterUptake << " " << seedmass << "\n";
};
*/
//------------------------------------------------------------------------------

/*
void plant::writePlantShort(ofstream * pout, double time) const
{
   (* pout) << time << " " << plantWaterContent << " " << growth << " " << lastTimePositiveGrowth << "\n";
};
*/
//------------------------------------------------------------------------------

/*
void plant::writePlantBin(ofstream * pout, double time) const
{
   (* pout).write(reinterpret_cast<const char*>(&time), sizeof(time));
   (* pout).write(reinterpret_cast<const char*>(&speciesID), sizeof(speciesID));
   (* pout).write(reinterpret_cast<const char*>(&growth), sizeof(time));
   (* pout).write(reinterpret_cast<const char*>(&lastTimePositiveGrowth), sizeof(time));
}
*/
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
  bool plant::mortality( const vector<plantParams>& pP, valarray<unsigned int> & plantDensity, double time, double deltatime, const string & simulationType,  valarray<double> & totalDeaths, valarray<double> & droughtDeaths)
  {
    bool allDead = false;
    unsigned int deadIndividuals;
    plantParams lpp = pP[speciesID];
    if (mass < 0.1 * lpp.seedlingsize) allDead = true;
    else
    {
      {
         // Age dependent mortality
         double upperAgeBound = time + deltatime - birthday;
         double lowerAgeBound = time - birthday;
         double ageMortalityFactor = pow((upperAgeBound/(lpp.maxAge * 365.0)), lpp.maxAgeGamma) - pow((lowerAgeBound/(lpp.maxAge * 365.0)), lpp.maxAgeGamma);
         double desiccationMortalityFactor = 0.0;
         if (desiccating == true)
         {
            double upperBound = time + deltatime - lastTimePositiveGrowth;
            double lowerBound; // = max(time, lastTimePositiveGrowth) - lastTimePositiveGrowth;
            if (upperBound > deltatime)
            {
              lowerBound = upperBound - deltatime;
            }
            else
            {
              lowerBound = 0.0;
            }
//            double desiccationMortalityFactor = lpp.desiccationAlpha * (pow(upperBound, lpp.desiccationGamma) - pow(lowerBound, lpp.desiccationGamma));     // wieso wird deltatime addiert?  weil mindestens ein Zeitschritt lang - der kommende - unter Trockenstress steht. Und gerade der lastTimePositiveGrowth gesetzt wurde. Ne Quatsch.
            desiccationMortalityFactor = pow((upperBound/lpp.desiccationAlpha), lpp.desiccationGamma) - pow((lowerBound/lpp.desiccationAlpha), lpp.desiccationGamma);     // wieso wird deltatime addiert?  weil mindestens ein Zeitschritt lang - der kommende - unter Trockenstress steht. Und gerade der lastTimePositiveGrowth gesetzt wurde. Ne Quatsch.
            mortalityprob = 1.0 - exp(- (desiccationMortalityFactor + lpp.baseMortalityRate * deltatime + ageMortalityFactor));
         }
         else
         {
            mortalityprob = 1.0 - exp(-(lpp.baseMortalityRate * deltatime + ageMortalityFactor));
         }
         if ( simulationType != "simple" ) {
             //deadIndividuals = gsl_ran_binomial(theRNG, mortalityprob, kohortSize);
             deadIndividuals = R::rbinom(kohortSize, mortalityprob);
             kohortSize -= deadIndividuals; // gsl_ran_binomial(theRNG, mortalityprob, kohortSize); // Mortality from unspecified causes
             if (deadIndividuals > 0 && mass > 0.1) {     // only established plant mortality is recorded
             totalDeaths[speciesID] += deadIndividuals;
             droughtDeaths[speciesID] += deadIndividuals *  desiccationMortalityFactor / (desiccationMortalityFactor + ageMortalityFactor);
             }
             plantDensity[cellID] -= deadIndividuals;
             if (kohortSize < 1) allDead = true;
         }
      }
    }
    return(allDead);
  };

//------------------------------------------------------------------------------

 void plant::waterdemand2( const vector<plantParams>& pP, double deltatime)
 {
    plantParams lpp = pP[speciesID];
    double maxWaterStorage = storagemass * lpp.succulenceFactor;
    double waterStorageDeficit = maxWaterStorage - plantWaterContent;


    double refillRate = lpp.maxUptakeRate * rootmass - lpp.respirationRate * mass / lpp.waterUseEfficiency;
    double refillTime = waterStorageDeficit / refillRate;

    if (refillTime < deltatime)
    {
      waterDemand = lpp.maxUptakeRate * rootmass * refillTime;
      waterDemand += (deltatime - refillTime) * min(lpp.maxUptakeRate * rootmass, potentialSpecificTranspirationRate * leafmass);
    }
    else
      waterDemand = lpp.maxUptakeRate * rootmass * deltatime;

    waterDemand *= kohortSize;

 };

//------------------------------------------------------------------------------

  //------------------------------------------------------------------------------

    void plant::production3 (  const vector<plantParams>& pP, double deltatime, double time  )
    {
       plantParams lpp = pP[speciesID];

       double maxWaterStorage = storagemass * lpp.succulenceFactor;
       double minTranspiration = lpp.respirationRate * deltatime * mass / lpp.waterUseEfficiency;

       double waterStorageDeficit = maxWaterStorage - plantWaterContent;
       double availableWater = waterUptake + plantWaterContent;

       if (waterUptake >= minTranspiration)
       {
         waterUptake -= minTranspiration;
         if (waterUptake >= waterStorageDeficit)
         {
           plantWaterContent = maxWaterStorage;
           waterUptake -= waterStorageDeficit;
           growth = (1.0 - lpp.rGrowth) * lpp.waterUseEfficiency * waterUptake;
         }
         else
         {
            plantWaterContent += waterUptake;
            growth = 0.0;
         }
         lastTimePositiveGrowth = time + deltatime;
         desiccating = false;
       }
       else
       {
          minTranspiration -= waterUptake;
          if (plantWaterContent >= minTranspiration)
          {
             plantWaterContent -= minTranspiration;
             growth = 0.0;
             lastTimePositiveGrowth = time + deltatime;
             desiccating = false;
          }
          else
          {
            growth = lpp.waterUseEfficiency * availableWater - lpp.respirationRate * deltatime * mass;
            plantWaterContent = 0.0;
            // if the amount of available water is sufficient for the transpiration of at least one day,
            // then the calculation of the dessiccation time starts new.
            double minTranspirationOneDay = lpp.respirationRate * 1.0 * mass / lpp.waterUseEfficiency;
            if (availableWater >= minTranspirationOneDay)
            {
               double sufficientWaterTime = availableWater / minTranspirationOneDay;
               lastTimePositiveGrowth = time + sufficientWaterTime;
             }
             desiccating = true;
          }
       }
    }

  //------------------------------------------------------------------------------

  void plant::setAmin(const vector<plantParams>& pP)
  {
     plantParams lpp = pP[speciesID];
     if (leafmass > 0.0000001 && mass > 0.0000001) {
        Amin = lpp.respirationRate / (lpp.waterUseEfficiency * leafmass / mass);
     }
     else {
        Amin = lpp.respirationRate / (lpp.waterUseEfficiency * lpp.pLeaf);
     }
     if (lpp.maxSpecificTranspirationRate < Amin)
       Amin = lpp.maxSpecificTranspirationRate;
  };

  //------------------------------------------------------------------------------

  void plant::seedproduction( const vector<plantParams>& pP, valarray<unsigned int>& newSeeds, unsigned int arenaSize)
  {
    plantParams lpp = pP[speciesID];
    newSeeds[cellID + speciesID * arenaSize] += floor(lpp.seedsPerSeedmass * seedmass * kohortSize);
    seedmass = 0.0;
  };

  //------------------------------------------------------------------------------

  void plant::potentialTranspirationNew( const vector<plantParams>& pP, const valarray<double>& theta, double ET, double cellArea, const vector< vector<unsigned int> > & rootPositions)
  {
     double maxLeafDensity =  (pP[speciesID]).maxLeafDensity;
     potentialSpecificTranspirationRate = (pP[speciesID]).maxSpecificTranspirationRate;
     if (potentialSpecificTranspirationRate * maxLeafDensity / 1000.0 > ET)
     {
       potentialSpecificTranspirationRate = ET * 1000.0 / maxLeafDensity;
     }
  }


void plant::waterorderNew(valarray<double>& plantwaterdemand,  const vector<plantParams>& pP, const valarray<double>& maxWaterUptakeRate, double soilDepth, double timestepLength, double thetaCrit, double thetaSat, double ET, double cellArea, const vector< vector<unsigned int> > & rootPositions)
{
    // potential transpiration rate
     const plantParams & lpp = pP[speciesID];
     potentialSpecificTranspirationRate = lpp.maxSpecificTranspirationRate;
     if (potentialSpecificTranspirationRate * lpp.maxLeafDensity / 1000.0 > ET)
     {
       potentialSpecificTranspirationRate = ET * 1000.0 / lpp.maxLeafDensity;
     }

    // water demand
      double rootmassInFirstCells;
      double rootmassInLastCell;
      unsigned int numberOfRootCells;
      unsigned int availableRootCells;

      rootmassInFirstCells = lpp.maxRootDensity * cellArea;
      numberOfRootCells = ceil(0.999 * rootmass / rootmassInFirstCells);
      rootmassInLastCell = rootmass - (numberOfRootCells - 1) * rootmassInFirstCells; // Can this go wrong if numberOfRootCells was exactly integer and did not have to be rounded up?
      availableRootCells = (rootPositions[cellID]).size();
      // if (numberOfRootCells > availableRootCells) writeErrorMessage("error.txt", "Too many root cells in potential Transpiration calculation in parent cell ", cellID);
      maxCurrentUptakeRate = 0.0;

      for (size_t i = 0; i < (numberOfRootCells - 1); ++i) {
          assert(rootPositions[cellID][i] >= 0 && rootPositions[cellID][i] < maxWaterUptakeRate.size());
          maxCurrentUptakeRate += rootmassInFirstCells * maxWaterUptakeRate[rootPositions[cellID][i]];
      }
      maxCurrentUptakeRate += rootmassInLastCell * maxWaterUptakeRate[rootPositions[cellID][numberOfRootCells - 1]];
      maxCurrentUptakeRate /= rootmass;
      maxCurrentUptakeRate *= lpp.maxUptakeRate;

     if (maxCurrentUptakeRate > 0.0) {
         double maxWaterStorage = storagemass * lpp.succulenceFactor;
         double waterStorageDeficit = maxWaterStorage - plantWaterContent;

         double refillRate = maxCurrentUptakeRate * rootmass - lpp.respirationRate * mass / lpp.waterUseEfficiency;
         double refillTime = waterStorageDeficit / refillRate;

         if (refillTime < timestepLength)
         {
            waterDemand = maxCurrentUptakeRate * rootmass * refillTime;
            waterDemand += (timestepLength - refillTime) * min(maxCurrentUptakeRate * rootmass, potentialSpecificTranspirationRate * leafmass);
         }
         else
            waterDemand = maxCurrentUptakeRate * rootmass * timestepLength;

        waterDemand *= kohortSize;


        for (size_t i = 0; i < (numberOfRootCells - 1); ++i) {
            plantwaterdemand[rootPositions[cellID][i]] += rootmassInFirstCells * maxWaterUptakeRate[rootPositions[cellID][i]] * lpp.maxUptakeRate / rootmass / maxCurrentUptakeRate * waterDemand;
        }
        plantwaterdemand[rootPositions[cellID][numberOfRootCells - 1]] += rootmassInLastCell * maxWaterUptakeRate[rootPositions[cellID][numberOfRootCells - 1]] * lpp.maxUptakeRate / rootmass / maxCurrentUptakeRate * waterDemand;
     }
     else
         waterDemand = 0.0;
}


void plant::wateruptakeNew(const valarray<double>& plantWaterDemand, const vector<plantParams>& pP, const valarray<double>& realizedPlantWaterDemand, double cellArea, const vector< vector<unsigned int> > & rootPositions, const valarray<double>& maxWaterUptakeRate)
{
    waterUptake = 0.0;
    if (waterDemand > 0.0) {
       const plantParams & lpp = pP[speciesID];
       double localPlantWaterDemand;
       double rootmassInFirstCells;
       double rootmassInLastCell;
       unsigned int numberOfRootCells;
       double multFac;
       rootmassInFirstCells = lpp.maxRootDensity * cellArea;
       numberOfRootCells = ceil(0.999 * rootmass / rootmassInFirstCells);
       rootmassInLastCell = rootmass - (numberOfRootCells - 1) * rootmassInFirstCells; // Can this go wrong if numberOfRootCells was exactly integer and did not have to be rounded up?

       double fixFac = 1.0 / rootmass * lpp.maxUptakeRate / maxCurrentUptakeRate * waterDemand / kohortSize;

       if (numberOfRootCells > (rootPositions[cellID]).size()) {
         // writeErrorMessage("error.txt", "Too many root cells in Water uptake calculation in parent cell ", cellID);
       }
       for (size_t i = 0; i < (numberOfRootCells - 1); ++i) {
            localPlantWaterDemand = plantWaterDemand[rootPositions[cellID][i]];
            if (localPlantWaterDemand > 0.0) {
                multFac = rootmassInFirstCells * maxWaterUptakeRate[rootPositions[cellID][i]] / localPlantWaterDemand * fixFac;
                waterUptake +=  realizedPlantWaterDemand[rootPositions[cellID][i]] * multFac;
            }
       }
       localPlantWaterDemand = plantWaterDemand[rootPositions[cellID][(numberOfRootCells - 1)]];
       if (localPlantWaterDemand > 0.0) {
            multFac = rootmassInLastCell * maxWaterUptakeRate[rootPositions[cellID][numberOfRootCells - 1]] / localPlantWaterDemand * fixFac;
            waterUptake += realizedPlantWaterDemand[rootPositions[cellID][(numberOfRootCells - 1)]] * multFac;
       }
    }
   trueWaterUptake = waterUptake;              // necessary for output functions
}


void plant::initAllocationNew( const vector<plantParams>& pP, double cellArea, const vector< vector<unsigned int> > & rootPositions)
{
     plantParams lpp = pP[speciesID];
     maxRootCells = ceil(pP[speciesID].pRoot * pP[speciesID].maxSize / (pP[speciesID].maxRootDensity * cellArea));
     maxLeafCells = ceil(pP[speciesID].pLeaf * pP[speciesID].maxSize /(pP[speciesID].maxLeafDensity * cellArea) );          //FIXME -- why??
     rootmass = lpp.pRoot*mass;
     storagemass = lpp.pStorage*mass;
     leafmass = lpp.pLeaf*mass;

     growth = 0.0;
}


void plant::allocationNew( const vector<plantParams>& pP, double cellArea, const vector<double>& leafFraction, const vector< vector<unsigned int> > & rootPositions, int year, double simulatedTime, double cellwidth, unsigned int width, unsigned int height, const list< valarray<unsigned int> > & newSeedsOfPastYears, const vector<unsigned int> & dispersalDates, bool reportEstablishment)
{
     const plantParams & lpp = pP[speciesID];
     if ( growth >= 0.0)
     {
        if ((rootPositions[cellID]).size() < maxRootCells) // i.e. the plant cannot grow to its maximum size because of space constraints
        {
           double maxGrowth = 0.99 * ((rootPositions[cellID]).size() * lpp.maxRootDensity * cellArea - rootmass) / lpp.pRoot;
           if (growth > maxGrowth) growth = maxGrowth;
        }
        mass += growth;
        if ( mass > lpp.maxSize ) {
           seedmass += mass - lpp.maxSize;
           mass = lpp.maxSize;
           if (timeToMaturity < 0.0) timeToMaturity = simulatedTime - birthday;
        }
     }
     else
     {
        mass += growth;
     }
     growth = 0.0;

     if (mass > 0.1 && recruitmentYear < 0)
     {
         recruitmentYear = year;

       /*
        if (reportEstablishment)
        {
         int birthyear = int(birthday /365.0);
         if (birthyear > 0 && birthyear < dispersalDates.size()) {
           int dispersalOccurredThisYear = dispersalDates[year] < simulatedTime;
           int dispersalOccurredNotInBirthyear = dispersalDates[birthyear] >= birthday;
           int yearsToGoBack = recruitmentYear - birthyear + dispersalOccurredThisYear + dispersalOccurredNotInBirthyear;
           if (recruitmentYear == birthyear && dispersalOccurredThisYear == 1 && dispersalOccurredNotInBirthyear == 0)
               yearsToGoBack = 0; // this is a special case of dispersal, germination and establishment in the same year
           if (yearsToGoBack < newSeedsOfPastYears.size())
           {
              list< valarray<unsigned int> >::const_iterator itNewSeeds = newSeedsOfPastYears.begin();
              for (size_t i = 0; i < yearsToGoBack; ++i)
                 itNewSeeds++;

            // now we have reached the correct newSeedsEntry
           gOutEstablishmentDistances << speciesID << " " << kohortSize << " ";
              unsigned int arenaSize = width * height;
              double dist = 0.0;
              for (size_t i = 0; i < arenaSize; ++i)
              {
                 if ((*itNewSeeds)[i + arenaSize*speciesID] > 0)
                 {
                    dist = getSquareDistanceNew(i, cellID, height, width);
                    gOutEstablishmentDistances << dist << " ";
                 }
              }
              gOutEstablishmentDistances << "\n";
           }
         }
        }
        */

     }

     rootmass = lpp.pRoot * mass;
     storagemass = lpp.pStorage * mass;

     double logLeafstemmass = log(lpp.pLeaf * mass);

     double dx = exp(floor(logLeafstemmass) + 1.0) - exp(floor(logLeafstemmass));

     int leafFractionIndex = int(floor(logLeafstemmass) - floor(log(lpp.seedlingsize * 0.001)));
     if (leafFractionIndex < 0) leafFractionIndex = 0;
     if (leafFractionIndex > leafFraction.size() - 1) leafFractionIndex = leafFraction.size() - 1;

     double fLeaf = leafFraction[leafFractionIndex];

     double dy;
     if (leafFractionIndex + 1 < leafFraction.size()) {
         dy = fLeaf - leafFraction[leafFractionIndex + 1];
     }
     else {
         dy = 0.0;
     }
     fLeaf -= (lpp.pLeaf * mass - exp(floor(logLeafstemmass))) * dy / dx;
     leafmass = fLeaf * lpp.pLeaf * mass;
}
