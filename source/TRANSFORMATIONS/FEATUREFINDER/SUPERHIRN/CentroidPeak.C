/*
 *  CentroidPeak.cpp
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  Copyright 2006 ETHZ, IMSB, Switzerland. All rights reserved.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#include <list>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>

// TODO Refactor
// The include IsotopicDist is just used to know one constant
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/IsotopicDist.h>

#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <iomanip>

using namespace std;

int	CentroidPeak::sfCentroidWindowWidth = 5; // Centroid window width

// public methods

// Default constructor	
CentroidPeak::CentroidPeak()
{
	fMass = 0.0;
	fIntensity = 0.0;
	fIsotopIdx = 0;
	fFittedIntensity = 0.0;
	fOrgIntensity = 0.0;
}

// Constructs a centroid peak object with mass / intensity values
CentroidPeak::CentroidPeak(
                           double pMass,
                           double pIntensity)
{
  fMass = pMass;
  fIntensity = pIntensity;
  fIsotopIdx = 0;
  fFittedIntensity = 0.0;
  fOrgIntensity = pIntensity;
}

// Constructs a centroid peak object with mass / intensity values
CentroidPeak::CentroidPeak(
                           double pMass,
                           double pIntensity,
                           double iRT)
{
  fMass = pMass;
  fIntensity = pIntensity;
  fIsotopIdx = 0;
  fFittedIntensity = 0.0;
  fOrgIntensity = pIntensity;
  fTr = iRT;
}


// Copy constructor
CentroidPeak::CentroidPeak(
	const CentroidPeak& pCentroidPeak) // Object to copy
{
	fMass = pCentroidPeak.fMass;
	fIntensity = pCentroidPeak.fIntensity;
	fIsotopIdx = pCentroidPeak.fIsotopIdx;
	fFittedIntensity = pCentroidPeak.fFittedIntensity;
	fSignalToNoise = pCentroidPeak.fSignalToNoise;
	fOrgIntensity = pCentroidPeak.fOrgIntensity;
  extraPeakInfo = pCentroidPeak.extraPeakInfo;
  fTr = pCentroidPeak.fTr;
}

// Copy constructor
CentroidPeak::CentroidPeak(
                           const CentroidPeak* pCentroidPeak) // Object to copy
{
  fMass = pCentroidPeak->fMass;
  fIntensity = pCentroidPeak->fIntensity;
  fIsotopIdx = pCentroidPeak->fIsotopIdx;
  fFittedIntensity = pCentroidPeak->fFittedIntensity;
  fSignalToNoise = pCentroidPeak->fSignalToNoise;
  fOrgIntensity = pCentroidPeak->fOrgIntensity;
  extraPeakInfo = pCentroidPeak->extraPeakInfo;
  fTr = pCentroidPeak->fTr;
}

// Destructor	
CentroidPeak::~CentroidPeak()
{
}
	
// Operators

// Copies values by assignemnt = operator
CentroidPeak& CentroidPeak::operator=(
	const CentroidPeak& pCentroidPeak) // Object to be assigned
{
  fMass = pCentroidPeak.fMass;
  fIntensity = pCentroidPeak.fIntensity;
  fIsotopIdx = pCentroidPeak.fIsotopIdx;
  fFittedIntensity = pCentroidPeak.fFittedIntensity;
  fSignalToNoise = pCentroidPeak.fSignalToNoise;
  fOrgIntensity = pCentroidPeak.fOrgIntensity;
  extraPeakInfo = pCentroidPeak.extraPeakInfo;
  fTr = pCentroidPeak.fTr;
  return (*this);
}

// Allows sorting objects in order of ascending mass
bool CentroidPeak::operator<(
	const CentroidPeak& pCentroidPeak) // Object to be assigned
{
	return (fMass < pCentroidPeak.fMass);
}

// subtract intensity
void CentroidPeak::subtractIntensity(
	double pIntensity) // intensity to be subtracted
{
	if (fIntensity<0.0) return;  // do nothing for small intensities
	
	if (abs(fIntensity-pIntensity)/fIntensity > IsotopicDist::sfIntensityCV) {
		fIntensity -= pIntensity; // subtract if difference is larger than stat variation (CV)
	} else {
		fIntensity = 0.0; // if difference not stat. significant, set to zero
	}
}


// Writes data to out stream using the << operator
std::ostream& operator<<(
                         std::ostream& pOut, // output stream 
	CentroidPeak& pCentroidPeak) // 
{
	pOut << std::fixed << std::setprecision(4) << pCentroidPeak.getMass() << " " << std::fixed << std::setprecision(2) << pCentroidPeak.getIntensity();
	return pOut;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DeconvPeak class methods

// public methods

// Default constructor	
DeconvPeak::DeconvPeak()
{
	fMass = 0.0;
	fIntensity = 0.0;
	fIsotopIdx = 0;
	fCharge = 0;
	fNrIsotopes = 0;
	fScore = 0.0;
	fC13MassError = 0.0;
}

// Constructs a centroid peak object with mass / intensity values
DeconvPeak::DeconvPeak(
                       double pMass,
                       double pIntensity,
                       int pCharge,
                       int pNrIsotopes,
                       double pScore,
                       double pC13MassError
                       )
{
  fMass = pMass;
  fIntensity = pIntensity;
  fIsotopIdx = 0;
  fCharge = pCharge;
  fNrIsotopes = pNrIsotopes;
  fScore = pScore;
  fC13MassError = pC13MassError;
}

// Copy constructor
DeconvPeak::DeconvPeak(
                       const DeconvPeak& pDeconvPeak) // Object to copy
{
  fMass = pDeconvPeak.fMass;
  fIntensity = pDeconvPeak.fIntensity;
  fIsotopIdx = pDeconvPeak.fIsotopIdx;
  fCharge = pDeconvPeak.fCharge;
  fNrIsotopes = pDeconvPeak.fNrIsotopes;
  fScore = pDeconvPeak.fScore;
  fC13MassError = pDeconvPeak.fC13MassError;
  fIsotopicPeaks = pDeconvPeak.fIsotopicPeaks;
  extraPeakInfo = pDeconvPeak.extraPeakInfo;
  
}


// Copy constructor
DeconvPeak::DeconvPeak(
                           const DeconvPeak* pDeconvPeak) // Object to copy
{
  fMass = pDeconvPeak->fMass;
  fIntensity = pDeconvPeak->fIntensity;
  fIsotopIdx = pDeconvPeak->fIsotopIdx;
  fCharge = pDeconvPeak->fCharge;
  fNrIsotopes = pDeconvPeak->fNrIsotopes;
  fScore = pDeconvPeak->fScore;
  fC13MassError = pDeconvPeak->fC13MassError;
  fIsotopicPeaks = pDeconvPeak->fIsotopicPeaks;
  extraPeakInfo = pDeconvPeak->extraPeakInfo;
}


// Destructor	
DeconvPeak::~DeconvPeak()
{
}

// Operators

// Copies values by assignemnt = operator
DeconvPeak& DeconvPeak::operator=(
                                  const DeconvPeak& pDeconvPeak) // Object to be assigned
{
  fMass = pDeconvPeak.fMass;
  fIntensity = pDeconvPeak.fIntensity;
  fIsotopIdx = pDeconvPeak.fIsotopIdx;
  fCharge = pDeconvPeak.fCharge;
  fNrIsotopes = pDeconvPeak.fNrIsotopes;
  fScore = pDeconvPeak.fScore;
  fC13MassError = pDeconvPeak.fC13MassError;
  fIsotopicPeaks = pDeconvPeak.fIsotopicPeaks;
  extraPeakInfo = pDeconvPeak.extraPeakInfo;
  return (*this);
}

// Writes data to out stream using the << operator
ostream& operator<<(
	ostream& pOut, // output stream 
	DeconvPeak& pDeconvPeak) // 
{
	pOut << (CentroidPeak&) pDeconvPeak;
	pOut << " " << pDeconvPeak.getCharge() << " " << fixed << setprecision(5) << pDeconvPeak.getC13MassError();
	pOut << " " << fixed << setprecision(2) << pDeconvPeak.getScore();
	return pOut;
}


// shows the info of the peak:
void DeconvPeak::show_info(){
  printf( "\tDeconvoluted Peak: mz=%.4f,I=%.4f\n",  fMass, fIntensity);

  if( !extraPeakInfo.empty() ){
    //cout<<"\t"<<extraPeakInfo<<endl;
  }
  
  if( ! fIsotopicPeaks.empty()){ 
    printf("\t");
    vector<CentroidPeak>::iterator I = fIsotopicPeaks.begin();
    while( I != fIsotopicPeaks.end() ) {
      printf("%0.4f(%0.0f[%0.0f]) ", (*I).getMass(), (*I).getFittedIntensity(), (*I).getOrgIntensity());
      I++;
    }
    printf("\n");
  }
}

// shows the info of the peak:
void CentroidPeak::show_info(){
  printf( "\tCentroidPeak: m/z=%.3f,I=%.4f\n",  fMass, fIntensity);
  
  if( !extraPeakInfo.empty() ){
    //cout<<"\t"<<extraPeakInfo<<endl;
  }
}


