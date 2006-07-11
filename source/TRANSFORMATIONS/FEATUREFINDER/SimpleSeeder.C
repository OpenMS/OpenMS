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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <iostream>

namespace OpenMS
{

	SimpleSeeder::SimpleSeeder(): 
		BaseSeeder(), is_initialised_(false),
		nr_seeds_(1), noise_threshold_(0)
	{
		name_ = SimpleSeeder::getName();
		defaults_.setValue("intensity_perc",0.03f);
		defaults_.setValue("min_intensity",-1.0f);
		param_ = defaults_;
	}

	SimpleSeeder::~SimpleSeeder(){}

  Index SimpleSeeder::nextSeed() throw (NoSuccessor)
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
		
		
		std::cout	<< "Processing seed " << nr_seeds_	<< " ("
							<< traits_->getPeakRt(*current_peak_) << ","
							<< traits_->getPeakMz(*current_peak_)
							<< ") with intensity " << traits_->getPeakIntensity(*current_peak_) << std::endl;
				
		nr_seeds_++;
		
		// we set the intensity threshold to a fixed percentage of the 5th largest intensity
		if (nr_seeds_ == 6) 
		{
			noise_threshold_ = param_.getValue("min_intensity");
			
			if (noise_threshold_ < 0.0)
			{	
				float int_perc = param_.getValue("intensity_perc");;
				noise_threshold_ = int_perc * traits_->getPeakIntensity(*current_peak_);
				//std::cout << " Setting noise threshold to relative value..." << std::endl;
			}
		
		}	
		std::cout << "SimpleSeeder: Intensity threshold for seeds is " << noise_threshold_ << std::endl;
		// if the intensity of the next seed is below this threshold,
		// seeding stops and we throw an execption.
		if (traits_->getPeakIntensity(*current_peak_) < noise_threshold_) 
		{
			std::cout << "Intensity below threshold: " << noise_threshold_ << std::endl;
			throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, *current_peak_); 	
		}
		
		// set flag
		traits_->getPeakFlag(*current_peak_) = FeaFiTraits::SEED;
				
		return *current_peak_++;
	}
	
	void SimpleSeeder::sort_() 
	{
		
		SimpleSeeder::IntensityLess comp = SimpleSeeder::IntensityLess::IntensityLess(traits_);
		sort(indizes_.begin(),indizes_.end(),comp);
		
		// we want to retrieve the peak with the highest
		// intensity first, therefore we reverse the order of peaks
		reverse(indizes_.begin(),indizes_.end());
			
	}

}
