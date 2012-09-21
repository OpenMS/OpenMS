// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Florian Zeller $
// $Authors: Markus Mueller $
// --------------------------------------------------------------------------
//
/*
 *  CentroidData.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

//#ifndef _CentroidData_h_
//#define _CentroidData_h_

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>

#include <list>
#include <vector>

#ifndef _USE_CENTROID_DATA_H
#define _USE_CENTROID_DATA_H

namespace OpenMS
{

class OPENMS_DLLAPI CentroidData
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

} // ns

#endif
