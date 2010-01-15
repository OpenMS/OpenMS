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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSumAgreeingIntensities.h>
#include <OpenMS/CONCEPT/Factory.h>



using namespace std;

namespace OpenMS
{
  BinnedSpectrumCompareFunctor::BinnedSpectrumCompareFunctor()
		: DefaultParamHandler(BinnedSpectrumCompareFunctor::getProductName())
  {
  }
  
  BinnedSpectrumCompareFunctor::BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor& source)
		:	DefaultParamHandler(source)
  {
  }
 
	BinnedSpectrumCompareFunctor::~BinnedSpectrumCompareFunctor()
	{
	}
 
  BinnedSpectrumCompareFunctor& BinnedSpectrumCompareFunctor::operator=(const BinnedSpectrumCompareFunctor& source)
  {
		if (this != &source)
		{
    	DefaultParamHandler::operator=(source);
		}
    return *this;
  }
 
  void BinnedSpectrumCompareFunctor::registerChildren()
  {
  	Factory<BinnedSpectrumCompareFunctor>::registerProduct(BinnedSharedPeakCount::getProductName(), &BinnedSharedPeakCount::create);
  	Factory<BinnedSpectrumCompareFunctor>::registerProduct(BinnedSpectralContrastAngle::getProductName(), &BinnedSpectralContrastAngle::create);
  	Factory<BinnedSpectrumCompareFunctor>::registerProduct(BinnedSumAgreeingIntensities::getProductName(), &BinnedSumAgreeingIntensities::create);
  }
  
  BinnedSpectrumCompareFunctor::IncompatibleBinning::IncompatibleBinning(const char* file, int line, const char* function, const char* message) throw()
          : BaseException(file, line, function, "BinnedSpectrumCompareFunctor::IncompatibleBinning", message)
  {
  } 
  BinnedSpectrumCompareFunctor::IncompatibleBinning::~IncompatibleBinning() throw()
  {	
  }

}
