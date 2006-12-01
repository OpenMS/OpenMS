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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/RankCorrelation.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <iostream>

namespace OpenMS
{

	RankCorrelation::RankCorrelation():	BaseQuality()
	{
		name_ = RankCorrelation::getName();
	}

	RankCorrelation::~RankCorrelation(){}

	double RankCorrelation::evaluate(const IndexSet& set, const BaseModel<2>& model)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				
		// store and sort intensities of model and data
		std::vector<IntensityType> data_intensities;
		std::vector<IntensityType> model_intensities;
		
		std::vector<unsigned int> ranks_data;
		std::vector<unsigned int> ranks_model;
		
		unsigned int rank_count = 0;
		
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			const DRawDataPoint<2>& peak = traits_->getPeak(*it);
			model_intensities.push_back( model.getIntensity(peak.getPosition()) );
			data_intensities.push_back( peak.getIntensity() );
			
			ranks_data.push_back(rank_count);
			ranks_model.push_back(rank_count);
			++rank_count;
		}
				
		// compute ranks of data
		RankCorrelation::RankComp data_comp = RankCorrelation::RankComp::RankComp(data_intensities);
		std::sort(ranks_data.begin(),ranks_data.end(),data_comp);
		
		// compute ranks of model
		RankCorrelation::RankComp model_comp = RankCorrelation::RankComp::RankComp(model_intensities);
		std::sort(ranks_model.begin(),ranks_model.end(),model_comp);
		
		int mu = (data_intensities.size() + 1) / 2; // average of ranks
		
		IntensityType sum_model_data   = 0;
				
		IntensityType sqsum_data   = 0;
		IntensityType sqsum_model = 0;
		
		for (unsigned int i=0; i<ranks_data.size();++i)
		{
			sum_model_data  += (ranks_data[i] - mu) *(ranks_model[i] - mu);
			
			sqsum_data   += (ranks_data[i] - mu) * (ranks_data[i] - mu);
			sqsum_model += (ranks_model[i] - mu) * (ranks_model[i] - mu);
		}
		
		// check for division by zero
		if ( ! sqsum_data || ! sqsum_model ) return 0;		
		
		double corr = sum_model_data / (  sqrt(sqsum_data) * sqrt(sqsum_model) ); 
		
		UnsignedInt df = set.size()-1;
		double t_stat = sqrt(df) * corr; 
		
		// t_stat follows t-distributuin with n-2 degrees of freedom
		pval_ = (1 - gsl_cdf_ugaussian_P(t_stat));	
				
// 		std::cout << "RankCorrelation: " << corr << std::endl;		
// 		std::cout << "RankCorrelation: t(1-a/2,n-1) = " << gsl_ran_gaussian_pdf(0.975, df ) << std::endl;
// 		std::cout << "RankCorrelation: t_stat = " << fabs(t_stat) << std::endl;
// 		std::cout << "RankCorrelation: P(t_stat) = " << pval_ << std::endl;		
		
		return ( fabs(corr));
	}
	
	double RankCorrelation::evaluate(const IndexSet& set, const BaseModel<1>& model, UnsignedInt dim)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				
		// store and sort intensities of model and data
		std::vector<IntensityType> data_intensities;
		std::vector<IntensityType> model_intensities;
		
		std::vector<unsigned int> ranks_data;
		std::vector<unsigned int> ranks_model;
		
		unsigned int rank_count = 0;
		
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			const CoordinateType coord = traits_->getPeak(*it).getPosition()[dim];;
			model_intensities.push_back( model.getIntensity( coord ) );
			data_intensities.push_back( traits_->getPeakIntensity(*it) );
			ranks_data.push_back(rank_count);
			ranks_model.push_back(rank_count);
			++rank_count;
		}
				
		// compute ranks of data
		RankCorrelation::RankComp data_comp = RankCorrelation::RankComp::RankComp(data_intensities);
		std::sort(ranks_data.begin(),ranks_data.end(),data_comp);
		
		// compute ranks of model
		RankCorrelation::RankComp model_comp = RankCorrelation::RankComp::RankComp(model_intensities);
		std::sort(ranks_model.begin(),ranks_model.end(),model_comp);
		
		int mu = (data_intensities.size() + 1) / 2; // average of ranks
		
		IntensityType sum_model_data   = 0;
				
		IntensityType sqsum_data   = 0;
		IntensityType sqsum_model = 0;
		
		for (unsigned int i=0; i<ranks_data.size();++i)
		{
			sum_model_data  += (ranks_data[i] - mu) *(ranks_model[i] - mu);
			
			sqsum_data   += (ranks_data[i] - mu) * (ranks_data[i] - mu);
			sqsum_model += (ranks_model[i] - mu) * (ranks_model[i] - mu);
		}
		
		// check for division by zero
		if ( ! sqsum_data || ! sqsum_model ) return 0;		
		return (sum_model_data / (sqsum_data * sqsum_model) );
		}
}
