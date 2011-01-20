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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>

using namespace std;
namespace OpenMS
{

  ThresholdMower::ThresholdMower()
    : PreprocessingFunctor()
  {
		setName(ThresholdMower::getProductName());
    defaults_.setValue("threshold", 0.05, "Intensity threshold, peaks which are below this threshold are thrown away");
		defaultsToParam_();
  }

  ThresholdMower::ThresholdMower(const ThresholdMower& source)
    : PreprocessingFunctor(source)
  {
  }
  
  ThresholdMower::~ThresholdMower()
  {
  }

  ThresholdMower& ThresholdMower::operator = (const ThresholdMower& source)
  {
		if (this != &source)
		{
    	PreprocessingFunctor::operator=(source);
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
