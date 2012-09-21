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
 *  Deisotoper.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef _DEISOTOPER_H_
#define _DEISOTOPER_H_

#include <OpenMS/CONCEPT/Types.h>

#include <list>
#include <iostream>

namespace OpenMS
{

class CentroidData;
class DeconvPeak;

class OPENMS_DLLAPI Deisotoper 
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

} // ns

#endif
