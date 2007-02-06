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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleSeeder.h>

#include <iostream>

using namespace std;

namespace OpenMS
{

	SimpleSeeder::SimpleSeeder()
		: BaseSeeder(), 
			is_initialized_(false),
			nr_seeds_(1)
	{
		setName(getProductName());
		
		defaults_.setValue("intensity_perc",0.03f);
		defaults_.setValue("min_intensity",0.0f);
		
		defaultsToParam_();
	}

	SimpleSeeder::~SimpleSeeder()
	{
	}

  SimpleSeeder::SimpleSeeder(const SimpleSeeder& rhs)
    : BaseSeeder(rhs),
    	is_initialized_(false)
  {
  }
  
  SimpleSeeder& SimpleSeeder::operator= (const SimpleSeeder& rhs)
  {
    if (&rhs == this) return *this;
    
    BaseSeeder::operator=(rhs);
    is_initialized_ = false;
    
    return *this;
  }

  FeaFiModule::IndexSet SimpleSeeder::nextSeed() throw (NoSuccessor)
	{
		if (!is_initialized_) 
		{
			// determine mininum intensity for last seed
			IntensityType noise_threshold  = param_.getValue("min_intensity");
			if (noise_threshold == 0.0)
			{
				IntensityType int_perc = param_.getValue("intensity_perc");;
				noise_threshold = int_perc * traits_->getData().getMaxInt();			
			}
			
			//reserve space for a quater of the peaks
			indizes_.reserve((std::vector<IDX>::size_type)round(traits_->getData().getSize() / 4.0));
			//fill indices for peaks above noise threshold
			IDX tmp = make_pair(0,0);
			while (tmp.first < traits_->getData().size())
			{
				tmp.second = 0;
				while (tmp.second < traits_->getData()[tmp.first].size())
				{
					if (traits_->getPeakIntensity(tmp)>noise_threshold)
					{
						indizes_.push_back(tmp);
					}
					++tmp.second;
				}
				++tmp.first;
			}
#ifdef DEBUG_FEATUREFINDER
 		std::cout	<< "Number of peaks above threshold (" << noise_threshold	<< "): " << indizes_.size() << endl;
#endif
			
			// sort index vector by intensity of peaks (highest first)
			sort(indizes_.rbegin(),indizes_.rend(),SimpleSeeder::IntensityLess::IntensityLess(traits_));
		
			current_peak_ = indizes_.begin();
			is_initialized_ = true;
		}
		
		// while the current peak is either already used or in a feature
		// jump to next peak...
		while (current_peak_ != indizes_.end() && traits_->getPeakFlag(*current_peak_) != FeaFiTraits::UNUSED) 
		{
			current_peak_++;
		}

		if (current_peak_ == indizes_.end()) 
		{
			throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, *current_peak_);
		}
		
		nr_seeds_++;
		
		// set flag
		traits_->getPeakFlag(*current_peak_) = FeaFiTraits::SEED;
		
		IndexSet result;
		result.insert( *current_peak_++ );
				
		return result;
	}
}
