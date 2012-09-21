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

#ifndef _USE_CENTROID_DATA_H
#define _USE_CENTROID_DATA_H

class RawData;
class CentroidPeak;

class CentroidData
{
public:
  
  static double	sfMassTolPpm;
  static double	sfMassTolDa;
  static double	sfMinIntensity;
  static double	sfIntensityFloor;
  
  // debugging parameters => used by other classes :( :( :(   
  static bool MonoIsoDebugging;
  static double DebugMonoIsoMassMin; 
  static double DebugMonoIsoMassMax; 
  
	CentroidData(int,RawData&, bool);
  CentroidData(int,RawData&, double, bool);
	virtual ~CentroidData();
  
	void get(std::list<CentroidPeak>&);
	void set(RawData&);
	void set(std::vector<double>&, std::vector<double>&);
  
	void setWidth(int pWidth) {fWindowWidth = pWidth;}
	int	 getWidth(){return fWindowWidth;}

	void setNoise(double);
	double getNoise() {return fNoise;}
  void removeNoise();
	
	bool getNextPeakGroup(std::list<CentroidPeak>::iterator&, std::list<CentroidPeak>::iterator&);	
	void resetPeakGroupIter();
  
  bool CENTROID_DATA_MODUS;
    
protected:
	
	void calcCentroids(RawData&);
  
	int	fWindowWidth;
	double fNoise;
  double fScanRetentionTime;
  std::list<CentroidPeak> fCentroidPeaks;
  std::list<CentroidPeak>::iterator fPeakGroupStart;			
};

std::ostream& operator<<(std::ostream&, CentroidData&);

#endif
