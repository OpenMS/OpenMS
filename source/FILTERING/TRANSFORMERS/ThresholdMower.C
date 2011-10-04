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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>

using namespace std;
namespace OpenMS
{
	
	ThresholdMower::ThresholdMower()
		: DefaultParamHandler("ThresholdMower")
	{
    defaults_.setValue("threshold", 0.05, "Intensity threshold, peaks below this threshold are discarded");
		defaultsToParam_();
	}
	
	ThresholdMower::~ThresholdMower()
  {
  }
	
  ThresholdMower::ThresholdMower(const ThresholdMower& source)
    : DefaultParamHandler(source)
  {
  }

  ThresholdMower& ThresholdMower::operator=(const ThresholdMower& source)
  {
		if (this != &source)
		{
    	DefaultParamHandler::operator=(source);
		}
    return *this;
  }

  void ThresholdMower::filterPeakSpectrum(PeakSpectrum& spectrum)
  {
    filterSpectrum(spectrum);
  }

  void ThresholdMower::filterPeakMap(PeakMap& exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
