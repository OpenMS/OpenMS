// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Correlation.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <iostream>

namespace OpenMS
{

	Correlation::Correlation():	BaseQuality()
	{
		name_ = Correlation::getName();
	}

	Correlation::~Correlation(){}

	double Correlation::evaluate(const IndexSet& set, const BaseModel<2>& model)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		typedef BaseModel<2>::IntensityType Intensity;
		Intensity cross_product_sum = 0;
		Intensity data_square_sum   = 0;
		Intensity model_square_sum  = 0;

		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			const DRawDataPoint<2>& peak = traits_->getPeak(*it);
			Intensity model_it = model.getIntensity(peak.getPosition());
			Intensity data_it  = peak.getIntensity();
			cross_product_sum += model_it * data_it;
			data_square_sum   += data_it  * data_it;
			model_square_sum  += model_it * model_it;
		}
		if ( ! data_square_sum || ! model_square_sum ) return 0;
		return (cross_product_sum * cross_product_sum) / (data_square_sum * model_square_sum);
	}
	
	double Correlation::evaluate(const IndexSet& set, const BaseModel<1>& model, UnsignedInt dim)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		typedef BaseModel<2>::IntensityType Intensity;
		Intensity cross_product_sum = 0;
		Intensity data_square_sum   = 0;
		Intensity model_square_sum  = 0;

		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			const DRawDataPoint<2>& peak = traits_->getPeak(*it);
			Intensity model_it = model.getIntensity(peak.getPosition()[dim]);
			Intensity data_it  = peak.getIntensity();
			cross_product_sum += model_it * data_it;
			data_square_sum   += data_it  * data_it;
			model_square_sum  += model_it * model_it;
		}
		if ( ! data_square_sum || ! model_square_sum ) return 0;
		return (cross_product_sum * cross_product_sum) / (data_square_sum * model_square_sum);
	}

}
