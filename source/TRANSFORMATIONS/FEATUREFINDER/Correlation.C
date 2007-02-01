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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Correlation.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

namespace OpenMS
{

	Correlation::Correlation()
		:	BaseQuality()
	{
		setName(getProductName());
		check_defaults_ = false;
	}

	Correlation::~Correlation()
	{
	}

	double Correlation::evaluate(const IndexSet& set, const BaseModel<2>& model)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		typedef BaseModel<2>::IntensityType Intensity;
		Intensity cross_product_sum = 0;
		Intensity data_square_sum   = 0;
		Intensity model_square_sum  = 0;
		
		Intensity model_avg = 0.0;
		Intensity data_avg    = 0.0;
		
		std::vector<Intensity> model_intensities(set.size());
		std::vector<Intensity> data_intensities(set.size());
		
		std::vector<Intensity>::iterator model_iter = model_intensities.begin();
		std::vector<Intensity>::iterator data_iter   = data_intensities.begin();

		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			*data_iter    = traits_->getPeakIntensity(*it);
			*model_iter = model.getIntensity(traits_->getPeakPos(*it));
			
			model_avg += *model_iter;
			data_avg    += *data_iter;
			
			++model_iter;
			++data_iter;			
		}
		
		// compute average intensities for data and model
		model_avg /= set.size();
		data_avg /= set.size();
		
		model_iter = model_intensities.begin();
		data_iter   = data_intensities.begin();
				
		for ( ;model_iter != model_intensities.end(); ++model_iter)
		{
			cross_product_sum += ( *model_iter - model_avg) * ( *data_iter - data_avg);
			data_square_sum    += ( *data_iter - data_avg)  * ( *data_iter - data_avg);
			model_square_sum  += ( *model_iter - model_avg)  * ( *model_iter - model_avg);			
		
			 ++data_iter;
		}
		
		if ( ! data_square_sum || ! model_square_sum ) return 0;
		double corr = cross_product_sum / sqrt(data_square_sum * model_square_sum);
				
		UnsignedInt df = set.size()-2;
		double t_stat = sqrt(df) * (corr / sqrt(1 - corr*corr)); 
		
		// t_stat follows t-distributoin with n-2 degrees of freedom
		pval_ = (1 - gsl_cdf_tdist_P(t_stat, df ));	
			
		return (fabs(corr));
	}
	
	double Correlation::evaluate(const IndexSet& set, const BaseModel<1>& model, UnsignedInt dim)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		typedef BaseModel<2>::IntensityType Intensity;
		Intensity cross_product_sum = 0;
		Intensity data_square_sum   = 0;
		Intensity model_square_sum  = 0;
		
		Intensity model_avg = 0.0;
		Intensity data_avg    = 0.0;
		
		std::vector<Intensity> model_intensities(set.size());
		std::vector<Intensity> data_intensities(set.size());
		
		std::vector<Intensity>::iterator model_iter = model_intensities.begin();
		std::vector<Intensity>::iterator data_iter   = data_intensities.begin();

		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			*data_iter    = traits_->getPeakIntensity(*it);
			*model_iter = model.getIntensity(traits_->getPeakPos(*it)[dim] );
			
			model_avg += *model_iter;
			data_avg    += *data_iter;
			
			++model_iter;
			++data_iter;			
		}
		
		// compute average intensities for data and model
		model_avg /= set.size();
		data_avg /= set.size();
		
		model_iter = model_intensities.begin();
		data_iter   = data_intensities.begin();
		
		for ( ;model_iter != model_intensities.end(); ++model_iter)
		{
			cross_product_sum += ( *model_iter - model_avg) * ( *data_iter - data_avg);
			data_square_sum    += ( *data_iter - data_avg)  * ( *data_iter - data_avg);
			model_square_sum  += ( *model_iter - model_avg)  * ( *model_iter - model_avg);			
		
			 ++data_iter;
		}
		
		if ( ! data_square_sum || ! model_square_sum ) return 0;
		
		double corr = cross_product_sum / sqrt(data_square_sum * model_square_sum);
		
 		return fabs(corr);
	}

}
