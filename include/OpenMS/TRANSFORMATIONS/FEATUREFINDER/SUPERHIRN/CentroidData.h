/*
 *  CentroidData.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  Copyright 2006 ETHZ, IMSB, Switzerland. All rights reserved.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

//#ifndef _CentroidData_h_
//#define _CentroidData_h_

#include <list>
#include <vector>

using namespace std;

#ifndef _USE_CENTROID_DATA_H
#define _USE_CENTROID_DATA_H


class RawData;
class CentroidPeak;
 


//using namespace std;

class CentroidData
{
public:
  
  static double	sfMassTolPpm;
  static double	sfMassTolDa;
  static double	sfMinIntensity;
  static double	sfIntensityFloor;
  
  // debugging parameters:
  static bool MonoIsoDebugging;
  static double DebugMonoIsoMassMin; 
  static double DebugMonoIsoMassMax; 
  
  

	CentroidData();
	CentroidData(int);
	CentroidData(int,RawData&, bool);
  CentroidData(int,RawData&, double, bool);
	CentroidData(int,vector<double>&,vector<double>&);
	virtual ~CentroidData();
  
	void get(list<CentroidPeak>&);
	void set(RawData&);
	void set(vector<double>&,vector<double>&);
  
	void setWidth(int pWidth) {fWindowWidth = pWidth;}
	int	 getWidth(){return fWindowWidth;}

	void setNoise(double);
	double getNoise() {return fNoise;}
  	void removeNoise();
	
	bool getNextPeakGroup(list<CentroidPeak>::iterator&,list<CentroidPeak>::iterator&);	
	void resetPeakGroupIter();
 // if data are in centroid modus:
 bool CENTROID_DATA_MODUS;
    
protected:
	
	void calcCentroids(RawData&);
  
	int		fWindowWidth;
	double	fNoise;
  double fScanRetentionTime;
	list<CentroidPeak> fCentroidPeaks;
	list<CentroidPeak>::iterator fPeakGroupStart;			
};

ostream& operator<<(ostream&, CentroidData&);


#endif
