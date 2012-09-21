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

#include <vector>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2ConsensusSpectrum.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ClusteredMS2ConsensusSpectrum.h>

namespace OpenMS
{

using namespace std;


////////////////////////////////////////////////
// constructor for the object ClusteredMS2ConsensusSpectrum:
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum( MS2ConsensusSpectrum* in):MS2ConsensusSpectrum(in){
  
  precursorMZ = in->getPrecursorMZ();
  TR = in->getTR();
  z = in->getPrecursorChrg();
  apexScan = in->getApexScan();
  
  this->addMS2ConsensusSpectrum( in );
}

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(double iPrecursorMZ, double iTR, int iChrg, int iApexScan) : MS2ConsensusSpectrum(iPrecursorMZ, iTR, iChrg, iApexScan ){}


////////////////////////////////////////////////
// constructor for the object ClusteredMS2ConsensusSpectrum:
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum( MS2Fragment* in):MS2ConsensusSpectrum(in){
  // store the scan numbers:
  MS2Scans.push_back( in->getApexScan() );
}


//////////////////////////////////////////////////
// class desctructor of ClusteredMS2ConsensusSpectrum
ClusteredMS2ConsensusSpectrum::~ClusteredMS2ConsensusSpectrum(){
  // ClusteredSpectra.clear();
}

//////////////////////////////////////////////////
// class copy constructor of ClusteredMS2ConsensusSpectrum
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum& tmp):MS2ConsensusSpectrum(tmp) {
  MS2Scans = tmp.MS2Scans;
}

//////////////////////////////////////////////////
// class copy constructor of ClusteredMS2ConsensusSpectrum
ClusteredMS2ConsensusSpectrum::ClusteredMS2ConsensusSpectrum(const ClusteredMS2ConsensusSpectrum* tmp):MS2ConsensusSpectrum(tmp){
  MS2Scans = tmp->MS2Scans;
}

//////////////////////////////////////////////////
// extracts fragments from a MS/MS spectra and inserts
// them into the Clustered MS/MS spectrum:
void ClusteredMS2ConsensusSpectrum::extractFragmentsFromSpectra( MS2ConsensusSpectrum* in  ){
  
  
  // go through the MS2 fragments and find the common one:
  multimap<double, MS2Fragment>::iterator P = in->getMS2FragmentPeakStart();
  while( P != in->getMS2FragmentPeakEnd() ){
    
    // fragment search mass:
    MS2Fragment* frag = &(P->second);
    double searchMZ = frag->getFragmentMz();
    
    MS2Fragment* matchedFrag = this->findMS2Fragment( searchMZ );    
    
    // in case where the fragment has been detected in a previous MS/MS
    // merge thes in one MS2 fragment
    if( matchedFrag != NULL ){
      this->mergeMS2Fragments( matchedFrag, frag );
    }
    else{
      this->addMS2Fragment( frag );
    }
    
    P++;
  }
}

//////////////////////////////////////////////////
// add a MS2 fragment:
void ClusteredMS2ConsensusSpectrum::addMS2ConsensusSpectrum( MS2ConsensusSpectrum* in ){
  
  // extract directly the MS/MS clustered spectrum
  extractFragmentsFromSpectra( in );
  
  // store the scan numbers:
  MS2Scans.push_back( in->getApexScan() );
  
  
}


//////////////////////////////////////////////////////
// merge a MS2 fragment into the target MS2 fragment:
void ClusteredMS2ConsensusSpectrum::mergeMS2Fragments(MS2Fragment* target, MS2Fragment* toMerge ){
  
  // sum up intensities:
  target->setFragmentPeakArea( target->getFragmentPeakArea() + toMerge->getFragmentPeakArea());
  
  // average m/z:
  target->setFragmentMz( ( target->getFragmentMz() + toMerge->getFragmentMz() ) / 2.0);
  
  // average m/z:
  target->setPrecursorMZ( ( target->getPrecursorMZ() + toMerge->getPrecursorMZ() ) / 2.0);
  
  
}

}
