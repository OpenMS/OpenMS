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
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
// 


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/featureLCprofile.h>

namespace OpenMS
{

using namespace std;

////////////////////////////////////////////////
// constructor for the object featureLCprofile:
featureLCprofile::featureLCprofile(){
}

////////////////////////////////////////////////
// constructor for the object featureLCprofile:
featureLCprofile::featureLCprofile(double apex_MZ, double apex_TR, double apex_Intensity, int apex_scan, int charge_state, double peak_area){
  
  // set the apex:
  apexMS1Signal.mass = apex_MZ;
  apexMS1Signal.TR = apex_TR;
  apexMS1Signal.intensity = apex_Intensity;
  apexMS1Signal.scan = apex_scan;
  apexMS1Signal.charge = charge_state;
  
  // set the computed peak area:
  LCelutionArea = peak_area;
  
}


////////////////////////////////////////////////
// constructor for the object featureLCprofile:
featureLCprofile::featureLCprofile(double apex_MZ, double apex_TR, int charge_state, double peak_area){
  
  // set the apex:
  apexMS1Signal.mass = apex_MZ;
  apexMS1Signal.TR = apex_TR;
  apexMS1Signal.intensity = -1;
  apexMS1Signal.scan = -1;
  apexMS1Signal.charge = charge_state;
  
  // set the computed peak area:
  LCelutionArea = peak_area;
}

//////////////////////////////////////////////////
// class desctructor of featureLCprofile
featureLCprofile::~featureLCprofile(){
  LCelutionSignals.clear();
  if( !outsideLCelutionSignals.empty() ){
    outsideLCelutionSignals.clear();
  }
}

//////////////////////////////////////////////////
// class copy constructor of featureLCprofile
featureLCprofile::featureLCprofile(const featureLCprofile& tmp){
  LCelutionSignals = tmp.LCelutionSignals;
  outsideLCelutionSignals = tmp.outsideLCelutionSignals;
  apexMS1Signal = tmp.apexMS1Signal;
  LCelutionArea = tmp.LCelutionArea;
}

//////////////////////////////////////////////////
// class copy constructor of featureLCprofile
featureLCprofile::featureLCprofile(const featureLCprofile* tmp){
  LCelutionSignals = tmp->LCelutionSignals;
  outsideLCelutionSignals = tmp->outsideLCelutionSignals;
  apexMS1Signal = tmp->apexMS1Signal;
  LCelutionArea = tmp->LCelutionArea;
}

//////////////////////////////////////////////////
// copy constructor:
featureLCprofile& featureLCprofile::operator=(const featureLCprofile& tmp){
  LCelutionSignals = tmp.LCelutionSignals;
  outsideLCelutionSignals = tmp.outsideLCelutionSignals;
  apexMS1Signal = tmp.apexMS1Signal;
  LCelutionArea = tmp.LCelutionArea;
  return *this;
}
    

/////////////////////////////////////////////////
// add / get signals:
void featureLCprofile::addMS1elutionSignal(  double mass, double intensity, int scan, int charge, double TR){
  MS1Signal tmp;
  tmp.mass = mass;
  tmp.intensity = intensity;
  tmp.scan = scan;
  tmp.charge = charge;
  tmp.TR = TR;
  LCelutionSignals.insert( std::make_pair( scan, tmp ) );

}

/////////////////////////////////////////////////
// add / get signals:
void featureLCprofile::addMS1elutionSignal( MS1Signal* in){
  LCelutionSignals.insert( std::make_pair( in->scan, *in ) );
}

/////////////////////////////////////////////////
// add / get signals:
void featureLCprofile::addOutsideMS1elutionSignal(  double mass, double intensity, int scan, int charge, double TR){
  MS1Signal tmp;
  tmp.mass = mass;
  tmp.intensity = intensity;
  tmp.scan = scan;
  tmp.charge = charge;
  tmp.TR = TR;
  outsideLCelutionSignals.insert( std::make_pair( scan, tmp ) );
  
}



/////////////////////////////////////////////////
// change all elution times by a factor:
void featureLCprofile::changeElutionTimesByFactor(double factor){
  apexMS1Signal.TR += factor;
  map<int, MS1Signal>::iterator P = getLCelutionSignalsStart();
  while( P != getLCelutionSignalsEnd() ){
  
    P->second.TR += factor;
    P++; 
  }

}
}
