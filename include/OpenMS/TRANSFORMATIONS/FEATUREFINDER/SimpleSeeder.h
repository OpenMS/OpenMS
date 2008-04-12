// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marcel Grunert, Clemens Groepl$
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>

#include <algorithm>
#include <vector>
#include <iostream>

namespace OpenMS
{
	/** 
		@brief Simple seeding class that uses the strongest peak as next seed.

		This class simply sorts the peaks according to intensity and proposes
		the highest peak, which is not yet included in a feature, as next seed.

		@ref SimpleSeeder_Parameters are explained on a separate page.

		@ingroup FeatureFinder
	*/
	template<class PeakType, class FeatureType> class SimpleSeeder :
		public FeaFiModule<PeakType, FeatureType>,
		public FeatureFinderDefs
	{
		public:
			typedef FeaFiModule<PeakType,FeatureType> Base;

			/// Constructor
			SimpleSeeder(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff) :
				Base(map,features,ff),
				initialized_(false)
			{
				this->setName("SimpleSeeder");				
				this->defaults_.setValue("intensity_perc",10.0, "Minimum percentage of the intensity of the largest peak that a seed has to have (used only if min_intensity is set to 0).", false);
        this->defaults_.setMinFloat("intensity_perc",0.0);
        this->defaults_.setMaxFloat("intensity_perc",100.0);
        this->defaults_.setValue("min_intensity",0.0, "Absolute value for the minimum intensity of a seed. If set to 0, a fixed percentage of the intensity of the largest peak is taken (see intensity_perc).", false);
        this->defaults_.setMinFloat("min_intensity",0.0);
			
				this->defaultsToParam_();
			}
			/// destructor 
			virtual ~SimpleSeeder()
			{
			}

			/// return the next seed
			IndexPair nextSeed() throw (NoSuccessor)
			{	
				if (!initialized_)
				{
					initialize_();
				}

				// while the current peak is either already used or in a feature jump to next peak...
				while (current_peak_ != indices_.end() && this->ff_->getPeakFlag(*current_peak_) == USED) 
				{
					++current_peak_;
				}

				if (current_peak_ == indices_.end()) 
				{
					throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, *current_peak_);
				}

				this->ff_->setProgress(current_peak_-indices_.begin());

				// set flag
				this->ff_->getPeakFlag(*current_peak_) = USED;

				return *(current_peak_++);
			} // nextSeed

		protected:
	
			void initialize_()
			{
				// determine mininum intensity for last seed
				typename FeatureType::IntensityType noise_threshold  = this->param_.getValue("min_intensity");
				if (noise_threshold == 0.0)
				{
					noise_threshold =
						typename FeatureType::IntensityType(this->param_.getValue("intensity_perc"))
						* (*this->map_).getMaxInt()
						/ 100.0;
				}
#ifdef DEBUG_FEATUREFINDER
				std::cout << "Threshold: " << noise_threshold << std::endl;			
#endif
			
				// fill indices_ for peaks above noise threshold
				IndexPair tmp = std::make_pair(0,0);
				while (tmp.first < (*this->map_).size())
				{
					tmp.second = 0;
					while (tmp.second < (*this->map_)[tmp.first].size())
					{
						if (this->getPeakIntensity(tmp)>noise_threshold)
						{
							indices_.push_back(tmp);
						}
						++tmp.second;
					}
					++tmp.first;
				}

#ifdef DEBUG_FEATUREFINDER
				std::cout	<< "Number of peaks above threshold (" << noise_threshold	<< "): " << indices_.size() << std::endl;
#endif

				// sort index vector by intensity of peaks (highest first)
				sort(indices_.begin(),indices_.end(),
						reverseComparator(Internal::IntensityLess<Base>(*this))
						);

				// progress logger
				this->ff_->startProgress(0, indices_.size() , "FeatureFinder");

				current_peak_ = indices_.begin();

				initialized_ = true;
				return;
			}

			/// contains the indizes 
			std::vector<IndexPair> indices_;

			/// Points to the next peak in the peak vector 
			std::vector<IndexPair>::const_iterator current_peak_;

			/// Flag that indicates of the indices are initialized
			bool initialized_;
			
		private:
			/// Not implemented
			SimpleSeeder();
			/// Not implemented
			SimpleSeeder& operator=(const SimpleSeeder&);
			/// Not implemented
			SimpleSeeder(const SimpleSeeder&);		

	}; // class SimpleSeeder

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H

