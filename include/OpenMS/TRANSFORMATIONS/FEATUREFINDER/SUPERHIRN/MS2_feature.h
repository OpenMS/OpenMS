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

#ifndef MS2_FEATURE_H
#define MS2_FEATURE_H

#include "MS2Fragment.h"
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ClusteredMS2ConsensusSpectrum.h>

namespace OpenMS
{

class OPENMS_DLLAPI MS2_feature : public ClusteredMS2ConsensusSpectrum {

  
  using ClusteredMS2ConsensusSpectrum::operator=;
    
  ////////////////////////////////////////////////
  // declaration of the private members:
  
private:
  
  int ID;
  
  ////////////////////////////////////////////////
  // declaration of the public members:
  
public:
  
  
  // class destructor
  ~MS2_feature();
  // class constructor
  MS2_feature();
  MS2_feature(MS2Fragment*);
  MS2_feature(double iPrecursorMZ, double iTR, int iChrg, int iApexScan);
  
  // class copy constructor
  MS2_feature(const MS2_feature&);
  // class copy constructor
  MS2_feature(const MS2_feature*);
  
  
  //////////////////////////////////////////////////
  // overload operators:
  //MS2_feature& operator=(const MS2_feature&);
  bool operator==(const MS2_feature&);
  MS2_feature& operator<=(const MS2_feature&);
  MS2_feature& operator>=(const MS2_feature&);
  MS2_feature& operator<(const MS2_feature&);
  MS2_feature& operator>(const MS2_feature&);
  
  
  // show info 
  void show_info();

  
  ///////////////////////////////
  // start here all the get / set
  // function to access the
  // variables of the class

  void setID( int in ){ID = in;};
  int getID( ){return ID;};

  
};

} // ns

#endif

    
