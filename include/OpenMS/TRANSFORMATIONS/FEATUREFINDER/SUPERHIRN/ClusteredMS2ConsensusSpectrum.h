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


#ifndef _CLUSTERED_MS2_CONSENSUS_SPECTRUM_H
#define _CLUSTERED_MS2_CONSENSUS_SPECTRUM_H

#include <OpenMS/CONCEPT/Types.h>

#include <vector>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2ConsensusSpectrum.h>

namespace OpenMS
{

class OPENMS_DLLAPI ClusteredMS2ConsensusSpectrum : public MS2ConsensusSpectrum{

  
   
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
  // stores the individual MS/MS spectra:
  std::vector<int> MS2Scans;
  
public:
  
  using MS2ConsensusSpectrum::operator=;
  
  // class destructor
  ~ClusteredMS2ConsensusSpectrum();
  
  // class constructor
  ClusteredMS2ConsensusSpectrum();
  ClusteredMS2ConsensusSpectrum( MS2Fragment*);
  ClusteredMS2ConsensusSpectrum( MS2ConsensusSpectrum* );
  ClusteredMS2ConsensusSpectrum(double, double, int, int);
  
  // class copy constructor
  ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum&);
  // class copy constructor
  ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  ClusteredMS2ConsensusSpectrum& operator=(const ClusteredMS2ConsensusSpectrum&);
  bool operator==(const ClusteredMS2ConsensusSpectrum&);
  ClusteredMS2ConsensusSpectrum& operator<=(const ClusteredMS2ConsensusSpectrum&);
  ClusteredMS2ConsensusSpectrum& operator>=(const ClusteredMS2ConsensusSpectrum&);
  ClusteredMS2ConsensusSpectrum& operator<(const ClusteredMS2ConsensusSpectrum&);
  ClusteredMS2ConsensusSpectrum& operator>(const ClusteredMS2ConsensusSpectrum&);
  
  
  
  
  
  
  
  //////////////////////////////////////////////////
  // trace the fragments across MS/MS scans runs:
  void constructClusteredConsenusSpectraFragments(MS2ConsensusSpectrum*);

  // add a MS2 Consensus Spectrum:
  void addMS2ConsensusSpectrum( MS2ConsensusSpectrum* );
  
  //////////////////////////////////////////////////
  // extracts fragments from a MS/MS spectra and inserts
  // them into the Clustered MS/MS spectrum:
  void extractFragmentsFromSpectra( MS2ConsensusSpectrum*);
  // merge a MS2 fragment into the target MS2 fragment:
  void mergeMS2Fragments(MS2Fragment* , MS2Fragment* );
  
    
  // plot all the consensus MS2 spectrum in one plot:
  void plotCombinedSpectra( );

  
  // remove outlier fragments based on their:
  // MS2Fragment::OutlierAttribute = ...
  // 1: retention time
  // 2: precursor mass
  // etc.
  void removeOutlierFragments();
  
  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class
 
  int getNumberOfSpectraScan(){return MS2Scans.size();};
  std::vector<int>::iterator getSpectraScanNumberStart(){return MS2Scans.begin();};
  std::vector<int>::iterator getSpectraScanNumberEnd(){return MS2Scans.end();};
    
};

} // ns

#endif

    
