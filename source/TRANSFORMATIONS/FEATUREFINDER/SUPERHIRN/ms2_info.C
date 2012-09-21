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

#include <string>
#include <map>
#include <list>
#include <vector>
#include <algorithm>
#include <string.h>
#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms2_info.h>

namespace OpenMS
{

using namespace std;

// static values:
const double ms2_info::_MONO_H = 1.00728;
const double ms2_info::_MONO_O = 15.99491;

//  monoisotopic masses of all amino acids:
const double ms2_info::mono_mass[26]={ 71.03711, 0.0, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406, 0.0, 128.09496, 113.08406, 131.04049, 114.04293, 0.0, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768, 0.0, 99.06841, 186.07931, 0.0, 163.06333, 0.0};
//  one letter code of all amino acids:
const char ms2_info:: AA[20]={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};

double ms2_info::MS2_TR_TOL;
bool ms2_info::THEO_MATCH_MODUS;
double ms2_info::MS2_MZ_PPM_TOLERANCE;




////////////////////////////////////////////////
// constructor for the object ms2_info:
ms2_info::ms2_info(){
  PEP_PROB = 0;
  DELTA_CN = 0;
  XCORR = 0;
  THEO_MZ = 0;
  MONO_MZ = 0;
  NEUTRAL_MR = 0;
  CHRG = 0;
  SCAN_START = 0;
  SCAN_END = 0;
  ID = -1;
  TR = -1.0;
}

////////////////////////////////////////////////
// constructor for the object ms2_info:
ms2_info::ms2_info(string IN_AC, string IN_SQ, float IN_PEP){
  PEP_PROB = IN_PEP;
  DELTA_CN = 0;
  XCORR = 0;
  THEO_MZ = 0;
  MONO_MZ = 0;
  NEUTRAL_MR = 0;
  CHRG = 0;
  ID = -1;
  TR = -1.0;
  SQ = IN_SQ;
  set_AC( IN_AC );  
  set_THEO_MASS_from_SQ();
  set_FULL_SQ();

}

////////////////////////////////////////////////
// constructor for the object ms2_info:
ms2_info::ms2_info(string IN_AC, string IN_SQ, int IN_CHRG, float IN_PEP){
  PEP_PROB = IN_PEP;
  DELTA_CN = 0;
  XCORR = 0;
  THEO_MZ = 0;
  MONO_MZ = 0;
  NEUTRAL_MR = 0;
  ID = -1;
  TR = -1.0;
  SQ = IN_SQ;
  set_AC( IN_AC );  
  CHRG = IN_CHRG;
  set_THEO_MASS_from_SQ();
  set_FULL_SQ();
}

////////////////////////////////////////////////
// constructor for the object ms2_info:
ms2_info::ms2_info(string IN_AC, string IN_SQ, float IN_PEP, int IN_CHRG, int IN_SCAN){
  PEP_PROB = IN_PEP;
  THEO_MZ = 0;
  MONO_MZ = 0;
  DELTA_CN = 0;
  XCORR = 0;
  NEUTRAL_MR = 0;
  ID = -1;
  TR = -1.0;
  SQ = IN_SQ;
  set_AC( IN_AC );  
  CHRG = IN_CHRG;
  SCAN_START = IN_SCAN;
  SCAN_END = IN_SCAN;
  set_THEO_MASS_from_SQ();
  set_FULL_SQ();

}


////////////////////////////////////////////////
// constructor for the object ms2_info:
ms2_info::ms2_info(int IN_ID){  
  ID = IN_ID;
  PEP_PROB = 0;
  THEO_MZ = 0;
  DELTA_CN = 0;
  XCORR = 0;
  MONO_MZ = 0;
  NEUTRAL_MR = 0;
  CHRG = 0;
  SCAN_START = 0;
  SCAN_END = 0;
  TR = -1.0;
}

//////////////////////////////////////////////////
// class desctructor of ms2_info
ms2_info::~ms2_info(){
  MOD_LIST.clear();
  FULL_SQ.clear();
  SQ.clear();
  AC.clear();
  TR = -1.0;
}

//////////////////////////////////////////////////
// class copy constructor of ms2_info
ms2_info::ms2_info(const ms2_info& tmp){
  ID = tmp.ID;
  PEP_PROB = tmp.PEP_PROB;
  DELTA_CN = tmp.DELTA_CN;
  XCORR = tmp.XCORR;
  THEO_MZ = tmp.THEO_MZ;
  MONO_MZ = tmp.MONO_MZ;
  NEUTRAL_MR = tmp.NEUTRAL_MR;
  CHRG = tmp.CHRG;
  SCAN_START = tmp.SCAN_START;
  SCAN_END = tmp.SCAN_END;
  TR = tmp.TR;
  AC = tmp.AC;
  SQ = tmp.SQ;
  PREV_AA = tmp.PREV_AA;
  FULL_SQ = tmp.FULL_SQ;
  MOD_LIST = tmp.MOD_LIST;
  ORIGINAL_INTERACT_FILE = tmp.ORIGINAL_INTERACT_FILE;
  MS2_TYPE_TAG = tmp.MS2_TYPE_TAG;
}

//////////////////////////////////////////////////
// class copy constructor of ms2_info
ms2_info::ms2_info(const ms2_info* tmp){
  ID = tmp->ID;
  PEP_PROB = tmp->PEP_PROB;
  DELTA_CN = tmp->DELTA_CN;
  XCORR = tmp->XCORR;
  THEO_MZ = tmp->THEO_MZ;
  MONO_MZ = tmp->MONO_MZ;
  NEUTRAL_MR = tmp->NEUTRAL_MR;
  CHRG = tmp->CHRG;
  SCAN_START = tmp->SCAN_START;
  SCAN_END = tmp->SCAN_END;
  TR = tmp->TR;
  AC = tmp->AC;
  SQ = tmp->SQ;
  PREV_AA = tmp->PREV_AA;
  FULL_SQ = tmp->FULL_SQ;
  MOD_LIST = tmp->MOD_LIST;
  ORIGINAL_INTERACT_FILE = tmp->ORIGINAL_INTERACT_FILE;
  MS2_TYPE_TAG = tmp->MS2_TYPE_TAG;
}

//////////////////////////////////////////////////
// copy constructor:
ms2_info& ms2_info::operator=(const ms2_info& tmp){
  ID = tmp.ID;
  PEP_PROB = tmp.PEP_PROB;
  DELTA_CN = tmp.DELTA_CN;
  XCORR = tmp.XCORR;
  THEO_MZ = tmp.THEO_MZ;
  MONO_MZ = tmp.MONO_MZ;
  NEUTRAL_MR = tmp.NEUTRAL_MR;
  CHRG = tmp.CHRG;
  SCAN_START = tmp.SCAN_START;
  SCAN_END = tmp.SCAN_END;
  AC = tmp.AC;
  SQ = tmp.SQ;
  TR = tmp.TR;
  PREV_AA = tmp.PREV_AA;
  FULL_SQ = tmp.FULL_SQ;
  MOD_LIST = tmp.MOD_LIST;
  ORIGINAL_INTERACT_FILE = tmp.ORIGINAL_INTERACT_FILE;
  MS2_TYPE_TAG = tmp.MS2_TYPE_TAG;
  return *this;
}

//////////////////////////////////////////////////
// equal operator:
bool ms2_info::operator==(const ms2_info& tmp){

  if(SQ == tmp.SQ)
    return true;
  else
    return false;
}


//////////////////////////////////////////////////
// calculates the theoretical mass from a sequence:
void ms2_info::set_THEO_MASS_from_SQ(){
  
  THEO_MZ = 0;
  int nb = 0;
  double TMP = 0;
    
  for(unsigned int POS = 0; POS < SQ.size(); POS++){
      
    // check for modification:
    map<int,double>::iterator P = MOD_LIST.find(POS);
    if(P == MOD_LIST.end()){
      
      /////////////////////////////////////
      // compute THEO M/Z
      switch( SQ[POS] ){
        
        case 'X':
          nb = (int)'L' - (int)'A';
          TMP += mono_mass[nb];
          break;
          
        default:
          nb = (int)SQ[POS] - (int)'A';
          
          if( ( nb >= 0 ) && ( nb < 26 ) ){
            TMP += mono_mass[nb];
          }
      }
      ////////////////////////////////
    }
    else{
      TMP += (*P).second;
    }
  }
  
  if( TMP > 0.0 ){
    // add the H20 (+18) and
    // add the charge state:
    TMP += (2.0 * _MONO_H + _MONO_O);
    TMP += double(CHRG) * _MONO_H; 
    TMP /= double(CHRG);
    THEO_MZ = TMP;
  }
  
}

//////////////////////////////////////////////////
// get the theoretical mass from a sequence:
double ms2_info::get_MONO_AA_MASS( int POS ){
  int nb = 0;
  switch( SQ[POS] ){
    case 'X':
      nb = (int)'L' - (int)'A';
      break;
    default:
      nb = (int)SQ[POS] - (int)'A';
  }
  
  return mono_mass[nb];
}

//////////////////////////////////////////////////
// compute the precursor 
void ms2_info::set_MONO_MZ(double IN){
  MONO_MZ = IN;
  NEUTRAL_MR = IN;
  NEUTRAL_MR *= double(CHRG);
  NEUTRAL_MR -= double(CHRG) * _MONO_H; 
}

//////////////////////////////////////////////////
// compute the precursor 
void ms2_info::set_NEUTRAL_MR(double IN){
  NEUTRAL_MR = IN;
  MONO_MZ = IN;
  MONO_MZ += double(CHRG) * _MONO_H; 
  MONO_MZ /= double(CHRG);
}

//////////////////////////////////////////////////
// add modification
void ms2_info::add_modification(int POS, double delta){
  
  map<int, double>::iterator M = MOD_LIST.find( POS );
  if( M != MOD_LIST.end() ){
    MOD_LIST.erase( M );
  }
  
  // add modification
  MOD_LIST.insert(pair<int,double>(POS,delta));
  // recompute the theoretical mass:
  set_THEO_MASS_from_SQ();
  // set the modified SQ:
  set_FULL_SQ();

}


/////////////////////////////////////////////////
// show info:
void ms2_info::show_info(){
  printf("\t\tMS2 ID: prec. m/z=%0.5f,theo. m/z=%0.5f,AC=%s,SQ=%s,P=%0.2f,scan=%d,tr=%0.2f,z=%d\n", MONO_MZ, THEO_MZ, get_AC().c_str(), get_TOTAL_SQ().c_str(), PEP_PROB, SCAN_START, TR, CHRG); 
}




/////////////////////////////////////////////////
// sets modificatied SQ:
void ms2_info::set_FULL_SQ(){
  
  FULL_SQ.clear();

  int pos = 0;
  for( unsigned int i=0; i<SQ.size(); i++){ 
    // insert normal AA letter:
    FULL_SQ += SQ[i];
    pos++;
    // check for modifications:
    map< int, double>::iterator F = find_Modification( i );
    if( F != get_Modification_list_end() ){
      char buffer[20];
      sprintf(buffer, "[%0.4f]",(*F).second);
      FULL_SQ += buffer;
      pos += strlen(buffer);
    }
  }  
}



//////////////////////////////////////////////////
// check if AC is in here:
bool ms2_info::find_AC( string IN ){
  
  vector<string>::iterator F = find( AC.begin(), AC.end(), IN );
  if( AC.end() == F ){
    return false;
  }
  else{
    return true;
  }
  
}



//////////////////////////////////////////////////
// check whethere proteotypic peptide:
bool ms2_info::get_PROTEO_TYPE(){
  
  if( AC.size() > 1 ){
    return false;
  }
  else{
    return true;
  }  
}


///////////////////////////////////////////////////
// add an AC to the ms2 scan:
void ms2_info::set_AC(string IN){
  
  vector<string>::iterator F = find( AC.begin(), AC.end(), IN);
  if( F == AC.end() ){
    AC.push_back( IN );
  }
}


///////////////////////////////////////////////////
// check if this AC or not:
bool ms2_info::compare_AC( string IN ){
  
  vector<string>::iterator F = find( AC.begin(), AC.end(), IN);
  if( F != AC.end() ){
    return true;
  }
  else{
    return false;
  }
}


///////////////////////////////////////////////////
// search a pattern in the  AC list:
bool ms2_info::search_AC_pattern(string IN){

  vector<string>::iterator F = AC.begin();
  while( F != AC.end() ){
    if( (*F).find( IN ) != string::npos){
      return true;
    }
    F++;
  }

  return false;
}


//////////////////////////////////////////////////
// check the tryptic state:
// 2: full tryptic
// 1: semi tryptic
// 0: non tryptic
int ms2_info::get_TRYPTIC_STATE(){
  
  int status = 0;
  // check C-terminus:
  if( ( SQ[ SQ.size() - 1] == 'R' ) || (  SQ[ SQ.size() - 1] == 'K' ) ){
    status++;
  }
  
  // check N-terminus:
  if( ( PREV_AA == "R" ) || (  PREV_AA == "K" ) ){
    status++;
  }
  return status;
}
}
