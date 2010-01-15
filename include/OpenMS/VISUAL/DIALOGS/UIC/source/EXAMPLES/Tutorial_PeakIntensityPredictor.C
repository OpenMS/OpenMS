// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Scherbart $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/PIP/PeakIntensityPredictor.h>

using namespace OpenMS;
using namespace std;

Int main() 
{
  //Create a vector for the predicted values that is large enough to hold them all
  vector<AASequence> peptides;
  peptides.push_back(AASequence("IVGLMPHPEHAVEK"));
  peptides.push_back(AASequence("LADNISNAMQGISEATEPR"));
  peptides.push_back(AASequence("ELDHSDTIEVIVNPEDIDYDAASEQAR"));
  peptides.push_back(AASequence("AVDTVR"));
  peptides.push_back(AASequence("AAWQVK"));
  peptides.push_back(AASequence("FLGTQGR"));
  peptides.push_back(AASequence("NYPSDWSDVDTK"));
  peptides.push_back(AASequence("GSPSFGPESISTETWSAEPYGR"));
  peptides.push_back(AASequence("TELGFDPEAHFAIDDEVIAHTR"));
  
  //Create new predictor model with vector of AASequences
  PeakIntensityPredictor model;
  
  //Perform prediction with LLM model
  vector< DoubleReal> predicted = model.predict(peptides);
  
  //for each element in peptides print sequence as well as corresponding predicted peak intensity value.
  for(Size i=0; i<peptides.size(); i++)
  {
    cout << "Intensity of " << peptides[i] << " is " << predicted[i] << endl;
  }

  return 0;
} //end of main

