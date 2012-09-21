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

#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>

namespace OpenMS
{

// outlier selection attribute:
// 1 -> retention time
// 2 -> preciursor MZ
// 
int MS2Fragment::OutlierAttribute = 1;

////////////////////////////////////////////////
// constructor for the object MS2Fragment:
MS2Fragment::MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ, double iIntensityArea,
                         int iScanStart, int iScanEnd, double iTrStart, double iTrEnd){

  precursorMZ = iPrecursorMZ;
  precursorCHRG = iPrecursorCHRG;
  TR = iTR;
  scan = iScan;
  z = iZ;
  fragmentMZ = iFragmentMZ;
  intensityArea = iIntensityArea;  
  scanStart = iScanStart;
  scanEnd = iScanEnd;
  trStart = iTrStart;
  trEnd = iTrEnd;
  
}



////////////////////////////////////////////////
// constructor for the object MS2Fragment:
MS2Fragment::MS2Fragment(double iPrecursorMZ, int iPrecursorCHRG, double iTR, int iScan, int iZ, double iFragmentMZ, double iIntensityArea ){
  
  precursorMZ = iPrecursorMZ;
  precursorCHRG = iPrecursorCHRG;
  TR = iTR;
  scan = iScan;
  z = iZ;
  fragmentMZ = iFragmentMZ;
  intensityArea = iIntensityArea;  
  scanStart = -1;
  scanEnd = -1;
  trStart = -1;
  trEnd = -1;
  
}


//////////////////////////////////////////////////
// class desctructor of MS2Fragment
MS2Fragment::~MS2Fragment(){
}

//////////////////////////////////////////////////
// class copy constructor of MS2Fragment
MS2Fragment::MS2Fragment(const MS2Fragment& tmp){
  precursorMZ = tmp.precursorMZ;
  precursorCHRG = tmp.precursorCHRG;
  TR = tmp.TR;
  scan = tmp.scan;
  z = tmp.z;
  fragmentMZ = tmp.fragmentMZ;
  intensityArea = tmp.intensityArea;  
  scanStart = tmp.scanStart;
  scanEnd = tmp.scanEnd;
  trStart = tmp.trStart;
  trEnd = tmp.trEnd;  
}

//////////////////////////////////////////////////
// class copy constructor of MS2Fragment
MS2Fragment::MS2Fragment(const MS2Fragment* tmp){
  precursorMZ = tmp->precursorMZ;
  precursorCHRG = tmp->precursorCHRG;
  TR = tmp->TR;
  scan = tmp->scan;
  z = tmp->z;
  fragmentMZ = tmp->fragmentMZ;
  intensityArea = tmp->intensityArea;  
  scanStart = tmp->scanStart;
  scanEnd = tmp->scanEnd;
  trStart = tmp->trStart;
  trEnd = tmp->trEnd;
  
}


//////////////////////////////////////////////////
// copy constructor:
MS2Fragment& MS2Fragment::operator=(const MS2Fragment& tmp){
  precursorMZ = tmp.precursorMZ;
  precursorCHRG = tmp.precursorCHRG;
  TR = tmp.TR;
  scan = tmp.scan;
  z = tmp.z;
  fragmentMZ = tmp.fragmentMZ;
  intensityArea = tmp.intensityArea;  
  scanStart = tmp.scanStart;
  scanEnd = tmp.scanEnd;
  trStart = tmp.trStart;
  trEnd = tmp.trEnd;
  return *this;
}

/////////////////////////////////////////////
// show info of the MS2 fragment
void MS2Fragment::show_info(){
  
  // print AMRT tag info:
  printf("\tm/z=%0.2f|precursor=%0.4f|TR=%0.2f:", 
         getFragmentMz(), getPrecursorMZ(), getTR());
  
  // scan/tr range:
  printf("[%d-%d],[%0.2f-%0.2f],",scanStart,scanEnd,trStart,trEnd);
  
  // area/intensity:
  printf("A=%0.1f",getFragmentPeakArea());
  
  printf("\n");
  
}

    
////////////////////////////////////////////
// get the attribute of the fragment
// according to which outliers are removed
double MS2Fragment::getOutlierDetectionAttribute(){
  
  switch( MS2Fragment::OutlierAttribute ){
    // retention time:
    case 1:
      return getTR();
      
      // precursor mass
    case 2:
      return getPrecursorMZ();
  
  
    default:
      break;
  }
  
  // otherwise use retention time:
  return getTR();
}

}
