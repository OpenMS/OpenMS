// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/MarkerMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ComplementMarker.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>

using namespace std;

namespace OpenMS
{

  /**
  IsotopeMarker, ComplementMarker, NeutralLossMarker are used per default<br>
  */
  MarkerMower::MarkerMower()
    : PreprocessingFunctor()
  {
		check_defaults_ = false;
		setName(MarkerMower::getProductName());
    //todo remove
    //insertmarker(new IsotopeMarker());
    //insertmarker(new ComplementMarker());
    //insertmarker(new NeutralLossMarker());
		defaultsToParam_();
  }

  MarkerMower::MarkerMower(const MarkerMower& source)
    : PreprocessingFunctor(source)
  {
  }

  MarkerMower::~MarkerMower()
  {
  }

  MarkerMower& MarkerMower::operator=(const MarkerMower& source)
  {
		if (this != &source)
		{
    	PreprocessingFunctor::operator=(source);
		}
    return *this;
  }

  void MarkerMower::filterPeakSpectrum(PeakSpectrum& spectrum)
  {
    filterSpectrum(spectrum);
  }

  void MarkerMower::filterPeakMap(PeakMap& exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

	
  ///@todo violates DefaultParamHandler interface (Andreas)
  void MarkerMower::insertmarker(PeakMarker* pm)
  {
    markers_.push_back(pm);
  }

}
