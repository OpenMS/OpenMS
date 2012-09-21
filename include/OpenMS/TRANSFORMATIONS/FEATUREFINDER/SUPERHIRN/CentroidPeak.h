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
 *  CentroidPeak.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef _CENTROIDPEAK_H_
#define _CENTROIDPEAK_H_

#include <OpenMS/CONCEPT/Types.h>

#include <ostream>
#include <cmath>
#include <vector>

namespace OpenMS
{

class OPENMS_DLLAPI CentroidPeak{
public:
  
  static int	sfCentroidWindowWidth;
  
  CentroidPeak();
  CentroidPeak(double,double);
  CentroidPeak(double,double, double);
  CentroidPeak(const CentroidPeak&);
  // Copy constructor
  CentroidPeak(const CentroidPeak*);
  
  CentroidPeak& operator=(const CentroidPeak&);
  
  bool operator<(const CentroidPeak&);
  
  virtual ~CentroidPeak();
  
  inline double getMass(){return fMass;}
  inline double getIntensity(){return fIntensity;}
  inline int	  getIsotopIdx(){return fIsotopIdx;}
  inline double getSignalToNoise(){return fSignalToNoise;}	
  inline double getFittedIntensity(){return fFittedIntensity;}
  inline double getOrgIntensity(){return fOrgIntensity;}
  inline std::string getExtraPeakInfo( ){ return extraPeakInfo;};
  inline double getRetentionTime( ){ return fTr;};
  
  inline void setMass(double pMass){fMass = pMass;}
  inline void setIntensity(double pIntensity){fIntensity = pIntensity;}
  inline void setIsotopIdx(double pIsotopIdx){fIsotopIdx = (int)pIsotopIdx;}
  inline void setSignalToNoise(double in){fSignalToNoise = in;}	
  inline void setFittedIntensity(double pFittedIntensity){fFittedIntensity = pFittedIntensity;}
  inline void setOrgIntensity(double pOrgIntensity){fOrgIntensity = pOrgIntensity;}
  inline void setExtraPeakInfo( std::string in){ extraPeakInfo = in;};
  inline void setRetentionTime(double in ){ fTr = in;};

  void show_info();
  void subtractIntensity(double);
  
protected:
  
  int		fIsotopIdx;
  double	fMass;
  double	fIntensity;
  double	fFittedIntensity;
  double	fOrgIntensity;
  double fTr;
  
  double fSignalToNoise;
  
  std::string extraPeakInfo;
};

std::ostream& operator<<(std::ostream&, CentroidPeak&);

// Class for deconvoluted isotopic patterns 
class DeconvPeak : public CentroidPeak{
public:
  
  DeconvPeak();
  DeconvPeak(double,double,int,int,double,double);
  DeconvPeak(const DeconvPeak&);
  DeconvPeak(const DeconvPeak*);
  
  DeconvPeak& operator=(const DeconvPeak&);
  
  virtual ~DeconvPeak();
  
  inline int getCharge(){return fCharge;}
  inline int getNrIsotopes(){return fNrIsotopes;}
  inline double getC13MassError(){return fC13MassError;}
  inline double getScore(){return fScore;}
  // shows the info of the peak:
  void show_info();
  
  inline std::vector<CentroidPeak> getIsotopicPeaks() {return fIsotopicPeaks;}
  
  inline void setCharge(int pCharge){fCharge = pCharge;}	
  inline void setC13MassError(double pC13MassError){fC13MassError = pC13MassError;}
  inline void setNrIsotopes(int pNrIsotopes){fNrIsotopes = pNrIsotopes;}
  inline void setScore(double pScore){fScore = pScore;}
  inline void setIsotopicPeaks(std::vector<CentroidPeak> pIsotopicPeaks) {
    fIsotopicPeaks = pIsotopicPeaks;
  }
  
protected:
  
  int						fCharge;
  int						fNrIsotopes;
  double					fC13MassError;
  double					fScore;
  std::vector<CentroidPeak>	fIsotopicPeaks;
};

std::ostream& operator<<(std::ostream&, DeconvPeak&);

} // ns

#endif
