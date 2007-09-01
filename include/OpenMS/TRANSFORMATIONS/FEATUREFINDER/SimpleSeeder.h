// -*- mode: C++; tab-width: 2; -*-
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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>

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
	template<class PeakType, class FeatureType>
  class SimpleSeeder
		: public FeaFiModule<PeakType, FeatureType>,
			public FeatureFinderDefs
  {
		public:
			typedef FeaFiModule<PeakType,FeatureType> Base;

			///Functor that allows to compare the indizes of two peaks by their intensity.
			template<class FeaFiModuleType>
			class IntensityLess 
			{			
				public:
					///Constructor that takes a FeaFiModule reference
					IntensityLess(const FeaFiModuleType& module)
						: module_(module)
					{
					}
					
					inline bool operator() (const typename FeatureFinderDefs::IDX& x, const typename FeatureFinderDefs::IDX& y)
					{
						return module_.getPeakIntensity(x) < module_.getPeakIntensity(y);
					}
					
				protected:
					const FeaFiModuleType& module_;
			};
			
			/// Default constructor
			SimpleSeeder(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff)
		: Base(map,features,ff), 
				is_initialized_(false),
				nr_seeds_(1)
			{
				this->setName("SimpleSeeder");
				
				this->defaults_.setValue("min_intensity",0.03f,"Absolute value for the minimum intensity of a seed. If set to 0, a fixed percentage of the intensity of the largest peak is taken (see intensity_perc).");
				this->defaults_.setValue("intensity_perc",0.0f,"Minimum percentage of the intensity of the largest peak that a seed has to have (used only if min_nitensity is set to 0).");
				
				this->defaultsToParam_();
			}

			/// destructor 
			virtual ~SimpleSeeder()
			{
			}
		
			/// return next seed 
			ChargedIndexSet nextSeed() throw (NoSuccessor)
			{
				if (!is_initialized_) 
				{
					// determine mininum intensity for last seed
					Real noise_threshold  = this->param_.getValue("min_intensity");
					if (noise_threshold == 0.0)
					{
						Real int_perc = this->param_.getValue("intensity_perc");;
						noise_threshold = int_perc * (*this->map_).getMaxInt();			
					}
					
					//reserve space for a quarter of the peaks
					indizes_.reserve((std::vector<IDX>::size_type)round((*this->map_).getSize() / 4.0));
					//fill indices for peaks above noise threshold
					IDX tmp = std::make_pair(0,0);
					while (tmp.first < (*this->map_).size())
					{
						tmp.second = 0;
						while (tmp.second < (*this->map_)[tmp.first].size())
						{
							if (this->getPeakIntensity(tmp)>noise_threshold)
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
					sort(indizes_.rbegin(),indizes_.rend(),typename IntensityLess<SimpleSeeder>::IntensityLess(*this));

					// progress logger
					this->ff_->startProgress(1, indizes_.size() , "FeatureFinder");
							
					current_peak_ = indizes_.begin();
					is_initialized_ = true;
				}
				
				// while the current peak is either already used or in a feature
				// jump to next peak...
				while (current_peak_ != indizes_.end() && this->ff_->getPeakFlag(*current_peak_) == USED) 
				{
					++current_peak_;
				}

				if (current_peak_ == indizes_.end()) 
				{
					throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, *current_peak_);
				}
				
				nr_seeds_++;
				this->ff_->setProgress(nr_seeds_);
				
				// set flag
				this->ff_->getPeakFlag(*current_peak_) = USED;
				
				ChargedIndexSet result;
				result.insert( *current_peak_++ );
						
				return result;
			}

		protected:
			/// contains the indizes 
			std::vector<IDX> indizes_;

			/// Indicates whether the vector of indizes is sorted 
			bool is_initialized_;
			
			/// Points to the next peak in the peak vector 
			std::vector<IDX>::const_iterator current_peak_;
			
			/// counts the number of seeds that we returned so far
			UInt nr_seeds_;

		private:
			/// Not implemented
			SimpleSeeder();
			/// Not implemented
			SimpleSeeder& operator=(const SimpleSeeder&);
			/// Not implemented
			SimpleSeeder(const SimpleSeeder&);
  };
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H

