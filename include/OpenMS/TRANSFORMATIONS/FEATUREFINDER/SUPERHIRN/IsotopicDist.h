/*
 *  IsotopicDist.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  Copyright 2006 ETHZ, IMSB, Switzerland. All rights reserved.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */


#ifndef _ISOTOPICDIST_h_
#define _ISOTOPICDIST_h_

class CentroidPeak;
class DeconvPeak;

#include <list>
#include <iostream>

using namespace std;

class IsotopicDist
{
public:

  static double	sfDetectableIsoFact;
  static double	sfIntensityCV;
  
	static void init();
	static bool getMatchingPeaks(list<CentroidPeak>::iterator, list<CentroidPeak>::iterator, int, double&, double, list<list<CentroidPeak>::iterator>&);
	static void subtractMatchingPeaks(list<list<CentroidPeak>::iterator>&, int, double, DeconvPeak&);
	static void getDistribution(double,double*&,double*&);
	static void getMassBounds(double,int,int,int,double,double&,double&);

	static bool getDebug() {return sfDebug;}
	static ostream* getDebugStream() {return sfStream;}
	
	static void setDebug(bool pDebug) {sfDebug = pDebug;}
	static void setDebugStream(ostream* pStream) {sfStream = pStream;}

protected:

	static int getIndex(double,int);

	static double	sfIsoDist10[96][20];
	static double	sfIsoDist50[96][20];
	static double	sfIsoDist90[96][20];
	static double	sfIsoMass10[96][20];
	static double	sfIsoMass50[96][20];
	static double	sfIsoMass90[96][20];
	static int		sfNrIsotopes[96];
 	static int		sfMaxMassIndex;
 	static int		sfMaxIsotopeIndex;
	static double	sfMinMass;
	static double	sfMaxMass;
	static double	sfMassStep;
	static bool		sfDebug;
	static ostream* sfStream;
};

// Returns lower and upper mass bounds for all isotopic peaks
inline void IsotopicDist::getMassBounds(
	double pMass, // m/z of monoisotopic peak
	int pCharge, // charge
	int pIdx, // index in isodist table
	int pIsotope, // 0 = mono isotopic peaks, 1 = C13 peak, ....
	double pTol, // Mass error (without calibration) of centroids
	double& pLower, // Iso dist masses as C array, monoisotopic peak is set to 0 mass
	double& pUpper) // Iso dist intensities as C array
{
	pLower = pMass + sfIsoDist10[pIdx][pIsotope]/pCharge - pTol;
	pUpper = pMass + sfIsoMass90[pIdx][pIsotope]/pCharge + pTol;
}

// Returns index in isotopic tables
inline int IsotopicDist::getIndex(
                                  double pMass, // m/z of monoisotopic peak
                                  int pCharge) // charge
{
  double diff;
  int idx;
  
  diff = (pMass*pCharge-sfMinMass)/sfMassStep;
  if (diff < 0) idx = 0;
  else if (diff < sfMaxMassIndex) idx = (int)((pMass*pCharge-sfMinMass)/sfMassStep);
  else idx = sfMaxMassIndex;
  
  return idx;
}



#endif
