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
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>

using namespace std;
namespace OpenMS
{
  BernNorm::BernNorm()
    : PreprocessingFunctor()
  {
		setName(BernNorm::getProductName());
    // values from the paper
    // they should be good for GoodDiff and Complements
    // IsotopeDiffs needs lower peaks
    defaults_.setValue("C1", 28.0, "C1 value of the normalization.", StringList::create("advanced"));
    defaults_.setValue("C2", 400.0, "C2 value of the normalization.", StringList::create("advanced"));
    defaults_.setValue("threshold", 0.1, "Threshold of the Bern et al. normalization."); // i.e. what is a significant peak
		defaultsToParam_();
  }

  BernNorm::BernNorm(const BernNorm& source)
    : PreprocessingFunctor(source)
  {
  }

  BernNorm::~BernNorm()
  {
  }

  BernNorm& BernNorm::operator = (const BernNorm& source)
  {
		if (this != &source)
		{
    	PreprocessingFunctor::operator = (source);
		}
    return *this;
  }

	void BernNorm::filterPeakSpectrum(PeakSpectrum& spectrum)
	{
		filterSpectrum(spectrum);
	}

	void BernNorm::filterPeakMap(PeakMap& exp)
	{
		for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
		{
			filterSpectrum(*it);
		}
	}
}
