// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

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
		//setName(BinnedRepCompareFunctor::getName())
    defaults_.setValue("filterwindow", 2.3);
		defaultsToParam_();
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
	
  /**
  does a check if comparison makes sense at all<br>
  counter-examples for useful comparisons would be spectra that have a mass 
  difference bigger than the <i>filterwindow</i>, spectra with different parent 
  charge states or bin-representation spectra with different binning parameters  
  */
	/*
	double BinnedRepCompareFunctor::filter(const ClusterSpectrum& a, const ClusterSpectrum& b) const
  {
    double filterwindow = (double)param_.getValue("filterwindow");
    double factor = 1;
    if (a.getParentionCharge() != b.getParentionCharge())
    {
      factor = 0;
    }
    if (usebins_)
    {
      stringstream ss;
      ss << a.getBinSize() << ":" << a.getBinSpread() << " != " << b.getBinSize() << ":" << b.getBinSpread() << " ids: " << a.id() << " " << b.id();
      // its not very informative to compare spectra with different bin sizes and spread
      if (fabs(a.getBinSize() - b.getBinSize()) > 1e-8 || a.getBinSpread() != b.getBinSpread())
      {
        throw ClusterSpectrum::WrongRepresentation(__FILE__, __LINE__, __PRETTY_FUNCTION__, ss.str().c_str());
      }
    }
    if (fabs(a.getParentMass() - b.getParentMass()) > filterwindow) return 0;
    return factor;
  }
	*/

}
