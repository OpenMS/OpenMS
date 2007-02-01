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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/RankCorrelation.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <iostream>

namespace OpenMS
{

	RankCorrelation::RankCorrelation()
		:	BaseQuality()
	{
		setName(getProductName());
		check_defaults_ = false;
	}

	RankCorrelation::~RankCorrelation()
	{
	}

	double RankCorrelation::evaluate(const IndexSet& set, const BaseModel<2>& model)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				
		// store and sort intensities of model and data
		std::vector<IntensityType> ranks_data;
		std::vector<IntensityType> ranks_model;
			
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			ranks_model.push_back( model.getIntensity( traits_->getPeakPos(*it) ) );
			ranks_data.push_back( traits_->getPeakIntensity(*it));
		}
		
// 		for (UnsignedInt i=0; i< ranks_data.size(); ++i)
// 		{
// 			std::cout << "Intensity data: " << ranks_data[i] << std::endl;
// 		}
// 		for (UnsignedInt i=0; i< ranks_model.size(); ++i)
// 		{
// 			std::cout << "Intensity model: " << ranks_model[i] << std::endl;
// 		}		
				
		// compute ranks of data
		std::sort(ranks_data.begin(),ranks_data.end());
		std::sort(ranks_model.begin(),ranks_model.end());
		
		computeRank_(ranks_data);
		computeRank_(ranks_model);
		
		int mu = (ranks_data.size() + 1) / 2; // mean of ranks

// 		for (UnsignedInt i=0; i< ranks_data.size(); ++i)
// 		{
// 			std::cout << "Rank data: " << ranks_data[i] << std::endl;
// 		}
// 		for (UnsignedInt i=0; i< ranks_model.size(); ++i)
// 		{
// 			std::cout << "Rank model: " << ranks_model[i] << std::endl;
// 		}		

		IntensityType sum_model_data   = 0;
				
		IntensityType sqsum_data   = 0;
		IntensityType sqsum_model = 0;
		
		for (unsigned int i=0; i<ranks_data.size();++i)
		{
			sum_model_data  += (ranks_data[i] - mu) *(ranks_model[i] - mu);
			sqsum_data   += (ranks_data[i] - mu) * (ranks_data[i] - mu);
			sqsum_model += (ranks_model[i] - mu) * (ranks_model[i] - mu);
		}
//		std::cout << "sqsum_data : " << sqsum_data << std::endl;
//		std::cout << "sqsum_model : " << sqsum_model << std::endl;
		
		// check for division by zero
		if ( ! sqsum_data || ! sqsum_model ) return 0;		
		
		double corr = sum_model_data / (  sqrt(sqsum_data) * sqrt(sqsum_model) ); 
		
		UnsignedInt df = set.size()-1;
		double t_stat = sqrt(df) * corr; 
		
		// t_stat follows Normal Gaussian distribution 
		pval_ = (1 - gsl_cdf_ugaussian_P(t_stat));	
				
		return ( fabs(corr));
	}
	
	double RankCorrelation::evaluate(const IndexSet& set, const BaseModel<1>& model, UnsignedInt dim)
	{
		if (!traits_) throw Exception::NullPointer(__FILE__, __LINE__, __PRETTY_FUNCTION__);
				
		// store and sort intensities of model and data
		std::vector<IntensityType> ranks_data;
		std::vector<IntensityType> ranks_model;
		
// 		std::cout << "Size of index set: " << set.size() << std::endl;
				
		for (IndexSet::const_iterator it=set.begin(); it!=set.end(); ++it)
		{
			const CoordinateType coord = traits_->getPeakPos(*it)[dim];
			ranks_data.push_back(traits_->getPeakIntensity(*it));
			ranks_model.push_back( model.getIntensity( coord ) );
		}
		
// 		std::cout << "Size of ranks set: " << ranks_data.size() << std::endl;
// 		for (UnsignedInt i=0; i< ranks_data.size(); ++i)
// 		{
// 			std::cout << "Intensity data: " << ranks_data[i] << std::endl;
// 		}
// 		for (UnsignedInt i=0; i< ranks_model.size(); ++i)
// 		{
// 			std::cout << "Intensity model: " << ranks_model[i] << std::endl;
// 		}		
				
		// compute ranks of data and model
		std::sort(ranks_data.begin(),ranks_data.end() );
		std::sort(ranks_model.begin(),ranks_model.end() );
		
		computeRank_(ranks_data);
		computeRank_(ranks_model);
		
		int mu = (ranks_data.size() + 1) / 2; // mean of ranks
		
	/*	for (UnsignedInt i=0; i< ranks_data.size(); ++i)
		{
			std::cout << "Rank data: " << ranks_data[i] << std::endl;
		}
		for (UnsignedInt i=0; i< ranks_model.size(); ++i)
		{
			std::cout << "Rank model: " << ranks_model[i] << std::endl;
		}		*/		
		
		IntensityType sum_model_data   = 0;
				
		IntensityType sqsum_data   = 0;
		IntensityType sqsum_model = 0;
		
		for (unsigned int i=0; i<ranks_data.size();++i)
		{
			sum_model_data  += (ranks_data[i] - mu) *(ranks_model[i] - mu);
			
			sqsum_data   += (ranks_data[i] - mu) * (ranks_data[i] - mu);
			sqsum_model += (ranks_model[i] - mu) * (ranks_model[i] - mu);
		}
		
// 		std::cout << "sqsum_data : " << sqsum_data << std::endl;
// 		std::cout << "sqsum_model : " << sqsum_model << std::endl;
				
			// check for division by zero
		if ( ! sqsum_data || ! sqsum_model ) return 0;		
		
		double corr = sum_model_data / (  sqrt(sqsum_data) * sqrt(sqsum_model) ); 
			
		return fabs(corr);
		}
}
