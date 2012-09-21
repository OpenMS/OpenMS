/*
 *  Deisotoper.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  Copyright 2006 ETHZ, IMSB, Switzerland. All rights reserved.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef _DEISOTOPER_H_
#define _DEISOTOPER_H_

class CentroidData;
class DeconvPeak;

#include <list>
#include <iostream>

class Deisotoper 
{
public:
  static int sfMinCharge;
  static int	sfMaxCharge;
  
	Deisotoper();
	Deisotoper(CentroidData&);
	virtual ~Deisotoper();
	
  std::list<DeconvPeak>& getDeconvPeaks() {return fDeconvPeaks;}
	
	void go(CentroidData&);
	void cleanDeconvPeaks();
	
	inline int getMinPeakGroupSize() {return fMinPeakGroupSize;}
	inline double getTheta() {return fTheta;}
	inline int getScanNumber() {return fScanNumber;}
	inline bool getShortReportFlag() {return fShortReportFlag;}

	inline void setMinPeakGroupSize(int pMinPeakGroupSize) {fMinPeakGroupSize = pMinPeakGroupSize;}
	inline void setTheta(double pTheta) {fTheta = pTheta;}
	inline void setScanNumber(int pScanNumber) {fScanNumber = pScanNumber;}
	inline void setShortReportFlag(bool pShortReportFlag) {fShortReportFlag = pShortReportFlag;}
	
protected:

  std::list<DeconvPeak> fDeconvPeaks;
 
	int		fMinPeakGroupSize;
	double	fTheta;
	int		fScanNumber;
	bool	fShortReportFlag;
};

std::ostream& operator<<(std::ostream&, Deisotoper&);

#endif
