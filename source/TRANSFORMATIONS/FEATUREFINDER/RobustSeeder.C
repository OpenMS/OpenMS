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


#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/RobustSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <iostream>

namespace OpenMS
{

	RobustSeeder::RobustSeeder(): 
		BaseSeeder(), is_initialised_(false), 
		noise_threshold_(0), nr_seeds_(1)
	{
		name_ = RobustSeeder::getName();
		defaults_.setValue("intensity_perc",0.03f);
		defaults_.setValue("min_intensity",-1.0f);
		param_ = defaults_;
	}

	RobustSeeder::~RobustSeeder(){}

  Index RobustSeeder::nextSeed() throw (NoSuccessor)
	{
		if (!is_initialised_) 
		{
			int endIndex = traits_->getNumberOfPeaks();
			for (int i=0;i<endIndex;i++) indizes_.push_back(i);
						
			sort_(); // sort index vector by intensity of peaks
			
			current_peak_ = indizes_.begin();
			is_initialised_ = true;
		}
						
		// while the current peak is either already used or in a feature
		// jump to next peak...
		while (current_peak_ != indizes_.end() 
		       && traits_->getPeakFlag(*current_peak_) != FeaFiTraits::UNUSED) 
		{
			current_peak_++;
		}
		
		if (current_peak_ == indizes_.end()) 
		{
			throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, *current_peak_);
		}		
		std::cout << " Processing seed " << nr_seeds_
		              << " with intensity " << traits_->getPeakIntensity(*current_peak_) << std::endl;
												
		// if the intensity of the next seed is below this threshold,
		// seeding stops and we throw an execption.
		if (traits_->getPeakIntensity(*current_peak_) < noise_threshold_) 
		{
			throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, *current_peak_); 	
		}
		
		// set flag
		traits_->getPeakFlag(*current_peak_) = FeaFiTraits::SEED;
		
		nr_seeds_++;
		return *current_peak_++;
	}
	
	void RobustSeeder::sort_() 
	{
		
		RobustSeeder::IntensityLess comp = RobustSeeder::IntensityLess::IntensityLess(traits_);
		sort(indizes_.begin(),indizes_.end(),comp);
		
		// we want to retrieve the peak with the highest
		// intensity first, therefore we reverse the order of peaks
		reverse(indizes_.begin(),indizes_.end());
			
		// compute noise threshold i.e. the minimum intensity of a seed
		// we estimate this threshold by computing the average of the
		// 300 data points with the largest intensity.
		unsigned int count = 0;

		noise_threshold_ = param_.getValue("min_intensity");
		if (noise_threshold_ < 0.0)
		{
			float int_perc = param_.getValue("intensity_perc");;
			std::vector<UnsignedInt>::const_iterator citer = indizes_.begin();
			FeaFiTraits::IntensityType intens_sum = 0;
			while (citer != indizes_.end() && count < 300)
			{
				intens_sum += traits_->getPeakIntensity(*current_peak_);
				citer++;
			}
			intens_sum /= count;
			noise_threshold_ = int_perc * intens_sum;
		}
		
		//std::cout << " Setting noise threshold to value " << noise_threshold_ << std::endl;
	}

}
