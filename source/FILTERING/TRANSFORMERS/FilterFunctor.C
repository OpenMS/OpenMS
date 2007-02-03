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
#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityDistBins.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/KellerQuality.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ParentFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakDensityFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakDiffBins.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakPosBins.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TradSeqQuality.h>

#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{
  FilterFunctor::FilterFunctor()
    : FactoryProduct(FilterFunctor::getProductName())
  {	
  }
  
  FilterFunctor::FilterFunctor(const FilterFunctor& source)
    : FactoryProduct(source)
  {
  }
  
  FilterFunctor& FilterFunctor::operator = (const FilterFunctor& source)
  {
		if (this != &source)
		{
    	FactoryProduct::operator=(source);
		}
    return *this;
  }
	
  FilterFunctor::~FilterFunctor()
  {
  }

	void FilterFunctor::registerChildren()
	{
		Factory<FilterFunctor>::instance()->registerProduct(ComplementFilter::getProductName(), &ComplementFilter::create);
    Factory<FilterFunctor>::instance()->registerProduct(GoodDiffFilter::getProductName(), &GoodDiffFilter::create);
    Factory<FilterFunctor>::instance()->registerProduct(IntensityBalanceFilter::getProductName(), &IntensityBalanceFilter::create);
    Factory<FilterFunctor>::instance()->registerProduct(IntensityDistBins::getProductName(), &IntensityDistBins::create);
    Factory<FilterFunctor>::instance()->registerProduct(NeutralLossDiffFilter::getProductName(), &NeutralLossDiffFilter::create);
    Factory<FilterFunctor>::instance()->registerProduct(IsotopeDiffFilter::getProductName(), &IsotopeDiffFilter::create);
    Factory<FilterFunctor>::instance()->registerProduct(KellerQuality::getProductName(), &KellerQuality::create);
    Factory<FilterFunctor>::instance()->registerProduct(ParentFilter::getProductName(), &ParentFilter::create);
    Factory<FilterFunctor>::instance()->registerProduct(TICFilter::getProductName(), &TICFilter::create);
    Factory<FilterFunctor>::instance()->registerProduct(PeakDensityFilter::getProductName(), &PeakDensityFilter::create);
    Factory<FilterFunctor>::instance()->registerProduct(PeakDiffBins::getProductName(), &PeakDiffBins::create);
    Factory<FilterFunctor>::instance()->registerProduct(PeakPosBins::getProductName(), &PeakPosBins::create);
    Factory<FilterFunctor>::instance()->registerProduct(TradSeqQuality::getProductName(), &TradSeqQuality::create);
	}
}
