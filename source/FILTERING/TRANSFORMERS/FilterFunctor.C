// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>

#include <OpenMS/CONCEPT/Factory.h>

namespace OpenMS
{
  FilterFunctor::FilterFunctor()
    : DefaultParamHandler("FilterFunctor")
  {	
  }
  
  FilterFunctor::FilterFunctor(const FilterFunctor& source)
    : DefaultParamHandler(source)
  {
  }
  
  FilterFunctor& FilterFunctor::operator = (const FilterFunctor& source)
  {
		if (this != &source)
		{
    	DefaultParamHandler::operator=(source);
		}
    return *this;
  }
	
  FilterFunctor::~FilterFunctor()
  {
  }

	void FilterFunctor::registerChildren()
	{
		Factory<FilterFunctor>::registerProduct(ComplementFilter::getProductName(), &ComplementFilter::create);
    Factory<FilterFunctor>::registerProduct(GoodDiffFilter::getProductName(), &GoodDiffFilter::create);
    Factory<FilterFunctor>::registerProduct(IntensityBalanceFilter::getProductName(), &IntensityBalanceFilter::create);
    Factory<FilterFunctor>::registerProduct(NeutralLossDiffFilter::getProductName(), &NeutralLossDiffFilter::create);
    Factory<FilterFunctor>::registerProduct(IsotopeDiffFilter::getProductName(), &IsotopeDiffFilter::create);
    Factory<FilterFunctor>::registerProduct(TICFilter::getProductName(), &TICFilter::create);
	}
}
