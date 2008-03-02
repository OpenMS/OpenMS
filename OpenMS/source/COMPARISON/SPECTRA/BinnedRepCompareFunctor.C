// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepCompareFunctor.h>

#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSpectrumContrastAngle.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSharedPeakCount.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSumAgreeingIntensities.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepMutualInformation.h>

#include <OpenMS/CONCEPT/Factory.h>

#include <cmath>
#include <sstream>

using namespace std;

namespace OpenMS
{
  BinnedRepCompareFunctor::BinnedRepCompareFunctor()
		: FactoryProduct(BinnedRepCompareFunctor::getProductName())
  {
  }
  
  BinnedRepCompareFunctor::BinnedRepCompareFunctor(const BinnedRepCompareFunctor& source)
		:	FactoryProduct(source)
  {
  }
 
	BinnedRepCompareFunctor::~BinnedRepCompareFunctor()
	{
	}
 
  BinnedRepCompareFunctor& BinnedRepCompareFunctor::operator = (const BinnedRepCompareFunctor& source)
  {
		if (this != &source)
		{
    	FactoryProduct::operator=(source);
		}
    return *this;
  }
 
	void BinnedRepCompareFunctor::registerChildren()
	{
    Factory<BinnedRepCompareFunctor>::registerProduct(BinnedRepSpectrumContrastAngle::getProductName(), &BinnedRepSpectrumContrastAngle::create);
    Factory<BinnedRepCompareFunctor>::registerProduct(BinnedRepMutualInformation::getProductName(), &BinnedRepMutualInformation::create);
    Factory<BinnedRepCompareFunctor>::registerProduct(BinnedRepSumAgreeingIntensities::getProductName(), &BinnedRepSumAgreeingIntensities::create);
    Factory<BinnedRepCompareFunctor>::registerProduct(BinnedRepSharedPeakCount::getProductName(), &BinnedRepSharedPeakCount::create);
	}

}
