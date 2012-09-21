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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/consensIsotopePattern.h>

#include <map>
#include <vector>
#include <math.h>

namespace OpenMS
{

using namespace std;

double consensIsotopePattern::FT_MZ_TOLERANCE;


////////////////////////////////////////////////
// constructor for the object consensIsotopePattern:
consensIsotopePattern::consensIsotopePattern(){
}

//////////////////////////////////////////////////
// class desctructor of consensIsotopePattern
consensIsotopePattern::~consensIsotopePattern(){
  
  isotopesTrace.clear(); 
  mzIsotopesStDev.clear(); 
  intensIsotopesStDev.clear(); 
  rawIsotopes.clear();
  
}

//////////////////////////////////////////////////
// class copy constructor of consensIsotopePattern
consensIsotopePattern::consensIsotopePattern(const consensIsotopePattern& tmp){
  isotopesTrace = tmp.isotopesTrace; 
  mzIsotopesStDev = tmp.mzIsotopesStDev; 
  intensIsotopesStDev = tmp.intensIsotopesStDev; 
  rawIsotopes = tmp.rawIsotopes;
}

//////////////////////////////////////////////////
// class copy constructor of consensIsotopePattern
consensIsotopePattern::consensIsotopePattern(const consensIsotopePattern* tmp){
  isotopesTrace = tmp->isotopesTrace; 
  mzIsotopesStDev = tmp->mzIsotopesStDev; 
  intensIsotopesStDev = tmp->intensIsotopesStDev; 
  rawIsotopes = tmp->rawIsotopes;
}

//////////////////////////////////////////////////
// copy constructor:
consensIsotopePattern& consensIsotopePattern::operator=(const consensIsotopePattern& tmp){
  isotopesTrace = tmp.isotopesTrace; 
  mzIsotopesStDev = tmp.mzIsotopesStDev; 
  intensIsotopesStDev = tmp.intensIsotopesStDev; 
  rawIsotopes = tmp.rawIsotopes;
  return *this;
}
    
// copied from simple_math
bool simple_math_compareMassValuesAtPPMLevel3( double mzA, double mzB, double PPM_TOLERANCE ){
  
  // take the average mass:
  double avMass = (mzA + mzB) / 2.0;
  
  // define the parts per million:
  double ppmValue = avMass / 1000000.00;
  double ppmDeltaTol = ppmValue * PPM_TOLERANCE;
  
  double deltaMass = fabs( mzA - mzB);
  if( deltaMass > ppmDeltaTol ){
    return false;
  }
  
  return true;
}
  
/////////////////////////////////////////////////
// order an isotope trace in the correct cluster:
void consensIsotopePattern::addIsotopeTrace( double mz, double intens){
  
  map<double, pair< vector<double>, vector<double> > >::iterator F = rawIsotopes.lower_bound( mz );
  bool match = false;
  if( F != rawIsotopes.end() ){
  
    // compute teh delta:
    if( simple_math_compareMassValuesAtPPMLevel3( mz, (*F).first, consensIsotopePattern::FT_MZ_TOLERANCE) ){
      (*F).second.first.push_back(mz);
      (*F).second.second.push_back(mz);
      match = true;
    }
    else if( F != rawIsotopes.begin() ){
      F--;
      if( simple_math_compareMassValuesAtPPMLevel3( mz, (*F).first, consensIsotopePattern::FT_MZ_TOLERANCE ) ){
        (*F).second.first.push_back(mz);
        (*F).second.second.push_back(mz);
        match = true;
      }
    
    }
  
  }

  if( !match ){
    vector< double > mzTmp;
    mzTmp.push_back( mz );
    vector< double > intensTmp;
    intensTmp.push_back( intens );
    rawIsotopes.insert( make_pair( mz, make_pair( mzTmp, intensTmp ) ) );
  }
  
  
}


/////////////////////////////////////////////////
// construc the consenus pattern:
void consensIsotopePattern::constructConsusPattern( ){
  
  
  map<double, pair< vector<double>, vector<double> > >::iterator I = rawIsotopes.begin();
  while( I != rawIsotopes.end() ){
    // condens a isotope peak trace:
    condensIsotopePattern( &(*I).second );    
    I++;
  }
  
}

// copied from simple_math
pair<double, double> simple_math_AVERAGE_and_STDEV(vector<double>* IN){
  
  double AVERAGE = 0;
  double STDEV = 0;
  
  if( IN->empty() ){
    return make_pair(AVERAGE, STDEV);  
  }
  
  if( IN->size() > 1 ){
    vector<double>::iterator START = IN->begin();  
    while( START != IN->end() ){
      AVERAGE += (*START);
      START++;
    }
    AVERAGE /= double( IN->size() );
    
    START = IN->begin();  
    while( START != IN->end() ){
      STDEV += ( ( AVERAGE - (*START) )*( AVERAGE - (*START) ) );
      START++;
    }
    STDEV /= double( IN->size() );
    STDEV = sqrt( STDEV );
    return make_pair(AVERAGE, STDEV);  
  }
  else{
    return make_pair( (*IN->begin()), 0.0);
  }
}

//////////////////////////////////////////////////
// condens the pattern, make averge peaks from the traces:
void consensIsotopePattern::condensIsotopePattern( pair< vector<double>, vector<double> >* in){

  // mz 
  pair<double, double> mz = simple_math_AVERAGE_and_STDEV( &(in->first) );
  // intens:
  pair<double, double> intens = simple_math_AVERAGE_and_STDEV( &(in->second) );
  
  isotopesTrace.insert( make_pair( mz.first, intens.first) ); 
  mzIsotopesStDev.push_back( mz.second ); 
  intensIsotopesStDev.push_back( intens.second ); 
  

}

}
