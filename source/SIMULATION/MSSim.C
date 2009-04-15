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
// $Maintainer: Stephan Aiche$
// $Authors: Stephan Aiche Chris Bielow$
// --------------------------------------------------------------------------

#include<OpenMS/SIMULATION/MSSim.h>

namespace OpenMS {

  MSSim::MSSim()
    : DefaultParamHandler("MSSim")
  {}

  MSSim::MSSim(const MSSim& source)
    : DefaultParamHandler(source)
  {}

  MSSim& MSSim::operator = (const MSSim& source)
  {
    return *this;
  }
  
  MSSim::~MSSim()
  {}
  
  void MSSim::simulate()
  {
    /*
      General progress should be 
        1. Digest Proteins
        2. Predict retention times
        3. add Post Translational modifications 
        4. predict detectibility 
        5. simulate ionization
        6. simulate the (lc)ms signal
        7. select features for MS2
        8. generate MS2 signals for selected features
     */
  
  }
}
