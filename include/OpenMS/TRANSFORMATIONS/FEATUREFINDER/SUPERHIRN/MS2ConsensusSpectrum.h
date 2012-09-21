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
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//


#ifndef _MS2_CONSENSUS_SPECTRUM_H
#define _MS2_CONSENSUS_SPECTRUM_H

#include <OpenMS/CONCEPT/Types.h>

#include <string>
#include <map>

namespace OpenMS
{

class OPENMS_DLLAPI MS2ConsensusSpectrum{

    
    ////////////////////////////////////////////////
    // declaration of the private members:

private:
  
  // stores the MS2 fragments:
  std::multimap<double, MS2Fragment> MS2FragmentPeaks;
  
  
protected:

    ////////////////////////////////////////////////
    // declaration of the public members:

    
  double startTR;
  double endTR;
  int z;
  int apexScan;
  int startScan;
  int endScan;
  
  
public:
  
  double precursorMZ;
  double TR;

     
  // mass to charge tolerance for MS2 trace level:
  static double MS2_MZ_TOLERANCE;
  
  // class destructor
  ~MS2ConsensusSpectrum();
  
  // class constructor
  MS2ConsensusSpectrum();
  MS2ConsensusSpectrum( MS2Fragment* );
  MS2ConsensusSpectrum(double iPrecursorMZ, double iTR, int iChrg, int iApexScan);

  // class copy constructor
  MS2ConsensusSpectrum(const MS2ConsensusSpectrum&);
  // class copy constructor
  MS2ConsensusSpectrum(const MS2ConsensusSpectrum*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  MS2ConsensusSpectrum& operator=(const MS2ConsensusSpectrum&);
  bool operator==(const MS2ConsensusSpectrum&);
  MS2ConsensusSpectrum& operator<=(const MS2ConsensusSpectrum&);
  MS2ConsensusSpectrum& operator>=(const MS2ConsensusSpectrum&);
  MS2ConsensusSpectrum& operator<(const MS2ConsensusSpectrum&);
  MS2ConsensusSpectrum& operator>(const MS2ConsensusSpectrum&);
  
  
  //////////////////////////////////////////////////////
  // copmute the similarity of the elution shape of the
  // MS2 fragment to this MS2 consensus spectrum
  double getLCElutionPeakSimilarity( MS2Fragment* );
    
  
  // add a MS2 fragment:
  void addMS2Fragment( MS2Fragment* );
  
  // process the stored fragments:
  void processConsenusSpectraFragments();
  
  // compute MS2 parameters
  void computeMS2SpectrumParameters();
  
  // remove outlier fragments based on their:
  // MS2Fragment::OutlierAttribute = ...
  // 1: retention time
  // 2: precursor mass
  // etc.
  void removeOutlierFragments();
  
  // remove H2O loss region of the MS2 spectra
  void removeWaterLossRegion( );
  // show MS2 spectrum info:
  void show_info( );

  
  // find a corresponding MS2 fragment
  MS2Fragment* findMS2Fragment( double );

  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
  
  // precursor mass:
  double getPrecursorMZ(){return precursorMZ;};

  // TR:
  double getTR(){return TR;};
  // start TR
  double getStartTR(){return startTR;};
  // end TR
  double getEndTR(){return endTR;};
  
  
  // set / get  the charge state of the precursor MZ:
  void setPrecursorChrg( int IN ){ z = IN;};
  int getPrecursorChrg(){ return z;};
  // apex scan:
  int getApexScan(){return apexScan;};
  // start scan
  int getStartScan(){return startScan;};
  // end scan
  int getEndScan(){return endScan;};
  // get the number of consensus fragments:
  int getNbMS2Fragments(){return (int) MS2FragmentPeaks.size();};
  
  // get the MS2 fragments list iterator:
  std::multimap<double, MS2Fragment>::iterator getMS2FragmentPeakStart(){return MS2FragmentPeaks.begin();};
  std::multimap<double, MS2Fragment>::iterator getMS2FragmentPeakEnd(){return MS2FragmentPeaks.end();};
  std::multimap<double, MS2Fragment>* getMS2FragmentMap(){return &MS2FragmentPeaks;};

  
};

} // ns

#endif

    
