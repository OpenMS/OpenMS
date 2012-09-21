/*
 *  CentroidData.cpp
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  Copyright 2006 ETHZ, IMSB, Switzerland. All rights reserved.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>

using namespace std;

double	CentroidData::sfMassTolPpm = 10.0; // mass tolerance in ppm between isotopes
double	CentroidData::sfMassTolDa = 0.01; // mass tolerance in Da between isotopes - total mass to = mass*fMassTolPpm/1000000 + fMassTolDa
double	CentroidData::sfMinIntensity = 0.0; // peak below this values are not considered as monoisotopic peaks
double	CentroidData::sfIntensityFloor = 1.0; // intensities below this value are considered as 0


bool CentroidData::MonoIsoDebugging = false;
double CentroidData::DebugMonoIsoMassMin =  318.8; 
double CentroidData::DebugMonoIsoMassMax =  319; 


// Default constructor
CentroidData::CentroidData() 
{
	fWindowWidth = 3;
	fNoise = 0.0;
}

// Use this constructor if MS data is provided later
CentroidData::CentroidData(
	int pWindowWidth)  // width of window for centroid algo
{
	fWindowWidth = pWindowWidth;
	fNoise = 0.0;
}

// Constructor for raw data, which is then centroided. Converts it to CentroidPeak objects
CentroidData::CentroidData(
	int pWindowWidth,  // width of window for quadratic fit
	RawData &pRawData, // Profile data object 
bool centMode // if data are in centroid modus
) 
{
	fWindowWidth = pWindowWidth;
	fNoise = 0.0;
 CENTROID_DATA_MODUS = centMode;

	set(pRawData);
}

// Constructor for raw data, which is then centroided. Converts it to CentroidPeak objects
CentroidData::CentroidData(
                           int pWindowWidth,  // width of window for quadratic fit
                           RawData &pRawData, // Profile data object
                           double iRT,
                           bool centMode // if data are in centroid modus
) 
{
  fScanRetentionTime = iRT;
  fWindowWidth = pWindowWidth;
  fNoise = 0.0;
  CENTROID_DATA_MODUS = centMode;
  
  set(pRawData);
}


// Constructor for already centroided data. This may be useful if centroiding is 
// already done by MS instrument software. Converts them to CentroidPeak objects
// The mass/intensity values must be sorted with respect to increasing mass
CentroidData::CentroidData(
	int pWindowWidth,  // width of window for quadratic fit
	vector<double> &pProfileMasses, // Mass sample values in profile mode
	vector<double> &pProfileIntens  // Intensity sample values in profile mode
) 
{
	fWindowWidth = pWindowWidth;
	fNoise = 0.0;
	
	set(pProfileMasses,pProfileIntens);
}

// destructor
CentroidData::~CentroidData() 
{
	fCentroidPeaks.clear();
}

// Operators

// Writes data to out stream using the << operator
ostream& operator<<(
                    ostream& pOut, // output stream 
                    CentroidData& pCentroidData) // 
{
  list<CentroidPeak> p;
  list<CentroidPeak>::iterator pi;
  
  pCentroidData.get(p);
  for (pi=p.begin();pi!=p.end();++pi) {
    pOut << *pi << endl;
  }
  
  return pOut;
}

// Public methods

// get centroid peak list
void CentroidData::get(list<CentroidPeak>& pCentroidPeaks){
  pCentroidPeaks = fCentroidPeaks;
}


// Sets list of already centroided mass and intensity values. Converts them to CentroidPeak objects 
// These values must be sorted with respect to increasing mass
void CentroidData::set(
                       vector<double>& pCentroidMasses, // Centroided masses
                       vector<double>& pCentroidIntens  // Centroided intensities
)
{
  vector<double>::iterator mi,hi;
  
  fCentroidPeaks.clear();
  
  for (mi=pCentroidMasses.begin(),hi=pCentroidIntens.begin();mi!=pCentroidMasses.end();++mi,++hi) {
    CentroidPeak peak(*mi,*hi);
    fCentroidPeaks.push_back(peak);		
  }
  
  resetPeakGroupIter();
}

// Sets raw profile data object. Converts it to CentroidPeak objects 
void CentroidData::set(
                       RawData &pRawData // Profile data
)
{
  calcCentroids(pRawData);
  resetPeakGroupIter();
}

// Calculates the precentile of all intensities and sets fNoise to it
void CentroidData::setNoise(
                            double pPrctile // percentile, which noiselevel is set to
)
{
  vector<double> intens;
  list<CentroidPeak>::iterator pi;
  
  for (pi=fCentroidPeaks.begin();pi!=fCentroidPeaks.end();++pi) {
    intens.push_back(pi->getIntensity());
  }
  
  sort(intens.begin(),intens.end()); // ascending order
  
  int len = (int)intens.size();
  
  /////////////////////////////////////////
  // modification Lukas:
  // only do if the vector length is longer than 0!
  if( len > 0 ){
    double prt = pPrctile*len/100.0;
    int i1 = (int)prt; 
    int i2 = i1+1;
    if (i2==len) i2--;
    
    // interpolate linearly between values if prt is not integer
    fNoise = (double)((prt-i1)*intens[i1]+(1-prt+i1)*intens[i2]);
  }
}

// Removes and deletes all CentroidePeak objects with a intensity below fNoise
// setNoise has to be called first
void CentroidData::removeNoise()
{
	list<CentroidPeak>::iterator pi;
	
	for (pi=fCentroidPeaks.begin();pi!=fCentroidPeaks.end();++pi) {
		if (fNoise > pi->getIntensity()) {
			pi = fCentroidPeaks.erase(pi);
			--pi;
		}	
	}
}

// A peak group is a set of peaks (ordered by mass) with maximal spacing of 1 + eps Da. If the gap is higher than
//  1 + eps Da, then a new peak group is formed. Deisotoping only within one peak group. Start and end iterators
// are set for the peakgroup
bool CentroidData::getNextPeakGroup(
	list<CentroidPeak>::iterator& pStart, // start position of peak group
	list<CentroidPeak>::iterator& pEnd) // end position of peak group
{
	list<CentroidPeak>::iterator pi,prev;
	double eps;  
	
	pi = fPeakGroupStart;
	prev = fPeakGroupStart;
	++pi;
	for (;pi!=fCentroidPeaks.end();++pi,++prev) {
		eps = CentroidData::sfMassTolPpm*pi->getMass()/1.0e6 + CentroidData::sfMassTolDa;
		if (abs(pi->getMass()-prev->getMass())>1.0+eps) {		
			break;
		}
	}
	
	pStart = fPeakGroupStart;
	pEnd = pi;
	fPeakGroupStart = pi; // store for next call
	
	return (fPeakGroupStart!=fCentroidPeaks.end());
}

// Resets the peak group iterator to start
void	CentroidData::resetPeakGroupIter()
{
	fPeakGroupStart = fCentroidPeaks.begin();
}

// Private methods


// Calculates centroides of peaks from raw data
void CentroidData::calcCentroids(	RawData &pRawData ) // Profile data object
// Calculates centroide data from profile data
{
  
  int i,hw,j;
  double cm,toti,min_dh;
  vector<double> masses,intens;
  
  pRawData.get(masses,intens);
  fCentroidPeaks.clear();
  
  ////////////////////////////////////////////
  
  if (CENTROID_DATA_MODUS) { // if data alread centroided in mzXML file
    for (i=0;i<(int)masses.size();i++) { 
      
      double inte = intens[i];
      double mz = masses[i];
      
      if( inte >= pRawData.LOW_INTENSITY_MS_SIGNAL_THRESHOLD ){
        CentroidPeak peak(mz,inte, fScanRetentionTime);
        fCentroidPeaks.push_back(peak);
      }
      
    }
  } 
  else {
    
    ////////////////////////////////////////////
    // centroid raw data
    min_dh = pRawData.LOW_INTENSITY_MS_SIGNAL_THRESHOLD; // min height
    hw = fWindowWidth/2;
    
    for (i=2;i<(int)masses.size()-2;i++) { 
      
      // Peak must be concave in the interval [i-2 .. i+2]
      if (intens[i]>min_dh && intens[i]>intens[i-1]+min_dh && intens[i]>=intens[i+1] && intens[i-1]>intens[i-2]+min_dh && intens[i+1]>=intens[i+2]) {
        
        // centroid mass:
        cm = 0.0;
        // total intensity:
        toti = 0.0;
        
        // double Tinte = intens[i];
        double Tmz = masses[i];
        
        /*
        if( MonoIsoDebugging ){
          if( ( DebugMonoIsoMassMin <= Tmz) && ( DebugMonoIsoMassMax >= Tmz) ){
            cout<<endl<<"*Centroid: "<<Tmz<<": "<<Tinte<<endl;
          }
        }
        */
        
        for (j= -hw;j<=hw;j++) {
          double inte = intens[i-j];
          double mz = masses[i-j];

          /*
          if( MonoIsoDebugging ){
            if( ( DebugMonoIsoMassMin <= Tmz) && ( DebugMonoIsoMassMax >= Tmz) ){
              cout<<"** add: "<<mz<<": "<<inte<<endl;
            }
          }
           */
          
          cm += inte*mz;
          toti += (double) intens[i-j];
        }
        cm = cm/toti;  // Centre of gravity = centroid
        
        // take the intensity at the apex of the profile peak:
        CentroidPeak peak(cm,intens[i],fScanRetentionTime);
       
        // Lukas: take summed intensity over all peaks:
        // CentroidPeak peak( cm, toti);
        
        if( MonoIsoDebugging ){
          if( ( DebugMonoIsoMassMin <= Tmz) && ( DebugMonoIsoMassMax >= Tmz) ){
            cout<<"*final: ";
            peak.show_info();
          }
        }
        
        
        fCentroidPeaks.push_back(peak);
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////
