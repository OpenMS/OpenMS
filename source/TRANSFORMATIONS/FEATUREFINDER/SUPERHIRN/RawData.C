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
/*
 *  RawData.cpp
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>

#include <iostream>
#include <iomanip>

namespace OpenMS
{

using namespace std;

// Constructor & destructor
///////////////////////////////////////////////////////////////////////////////


// More general constructor
RawData::RawData(
                 vector<double>& pMassValues, // mass sample values 
                 vector<double>& pIntensValues // intensity sample values 
                 )
{
  fProfileMasses = pMassValues;
  fProfileIntens = pIntensValues;
  LOW_INTENSITY_MS_SIGNAL_THRESHOLD = 1.0; 
}

// Destructor
RawData::~RawData()
{
}

// Operators

// Writes data to out stream using the << operator
ostream& operator<<(
                    ostream& pOut, // output stream 
                    RawData& pRawData) // 
{
  vector<double> m,h;
  vector<double>::iterator mi,hi;
  
  pRawData.get(m,h);
  for (mi=m.begin(),hi=h.begin();mi!=m.end();++mi,++hi) {
    pOut << fixed << setprecision(4) << *mi << " " << fixed << setprecision(2) << *hi << endl;
	}
  
  return pOut;
}


// Public methods

// Retrieve raw data as mass and intensity vectors
void RawData::get(
                  vector<double> &pProfileMasses, // Mass sample values in profile mode
                  vector<double> &pProfileIntens  // Intensity sample values in profile mode
                  ){
  
  pProfileMasses = fProfileMasses;
  pProfileIntens = fProfileIntens;

}

// Set raw data as mass and intensity vectors
void RawData::set(	vector<double> &pProfileMasses, // Mass sample values in profile mode
                   vector<double> &pProfileIntens  // Intensity sample values in profile mode
                   ){
  
	fProfileMasses = pProfileMasses;
	fProfileIntens = pProfileIntens;

}
}
