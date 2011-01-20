// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLESEEDER_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>

#include <algorithm>
#include <vector>
#include <iostream>

namespace OpenMS
{
	/** 
		@brief Simple seeding class that uses the strongest peak as next seed.

		This class simply sorts the peaks according to intensity and proposes
		the highest peak, which is not yet included in a feature, as next seed.

		@htmlinclude OpenMS_SimpleSeeder.parameters

		@ingroup FeatureFinder
	*/
	template<class PeakType, class FeatureType> class SimpleSeeder :
		public FeaFiModule<PeakType, FeatureType>,
		public FeatureFinderDefs
	{
		public:
			typedef FeaFiModule<PeakType,FeatureType> Base;
			typedef MSExperiment<PeakType> MapType;
      
			/// Constructor
			SimpleSeeder(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff) :
				Base(map,features,ff),
				initialized_(false)
			{
				this->setName("SimpleSeeder");				

        this->defaults_.setValue("min_intensity",0.0, "Absolute value for the minimum intensity required for a seed.");
        this->defaults_.setMinFloat("min_intensity",0.0);
        this->defaults_.setValue("signal_to_noise", 10.0, "Minimal required SignalToNoise (S/N) ratio for a seed.");
        this->defaults_.setMinFloat("signal_to_noise",0.0);
        
        //this->subsections_.push_back("SignalToNoiseEstimationParameter");
				SignalToNoiseEstimatorMedian < typename MapType::SpectrumType > sne; // make sure this is the same as in pick()!
				this->defaults_.insert ("SignalToNoiseEstimationParameter:", sne.getDefaults());
        
        this->defaultsToParam_();
			}
			/// destructor 
			virtual ~SimpleSeeder()
			{
			}

			/// return the next seed
			IndexPair nextSeed()
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
					// if no seed was found:
					if (indices_.size()==0) throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, IndexPair());
					else throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, *(current_peak_-1));
				}

				this->ff_->setProgress(current_peak_-indices_.begin());

				// set flag
				this->ff_->getPeakFlag(*current_peak_) = USED;

				return *(current_peak_++);
			} // nextSeed

		protected:
      
  		void initialize_()
			{
				// determine mininum intensity and signal-to-noise parameter for last seed
				typename FeatureType::IntensityType noise_threshold  = this->param_.getValue("min_intensity");
 			  typename FeatureType::IntensityType sn  = this->param_.getValue("signal_to_noise");
	     
#ifdef DEBUG_FEATUREFINDER
				std::cout << "Intensity threshold: " << noise_threshold << std::endl;	
				std::cout << "S/N: " << sn << std::endl;
#endif

				// fill indices_ for peaks above noise threshold and S/N
				IndexPair tmp = std::make_pair(0,0);
				if (sn == 0)
				{
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
				}
				else
				{
          SignalToNoiseEstimatorMedian < typename MapType::SpectrumType > estimator;
          Param param(this->param_.copy("SignalToNoiseEstimationParameter:",true));
          estimator.setParameters(param);
          
          for (typename MapType::ConstIterator it = (*this->map_).begin(); it != (*this->map_).end(); ++it)
          {
            estimator.init(it->begin(),it->end()); 
            tmp.second = 0;   
            for (typename MapType::SpectrumType::ConstIterator spec = it->begin(); spec != it->end(); ++spec)
            {
              if (estimator.getSignalToNoise(spec) > sn && this->getPeakIntensity(tmp) > noise_threshold) 
              {
                indices_.push_back(tmp);
              }
              ++tmp.second;
            }
            ++tmp.first;
          }
 			  }

#ifdef DEBUG_FEATUREFINDER
				std::cout	<< "Number of peaks above threshold (" << noise_threshold	<< ") and S/N (" << sn << "): " << indices_.size() << std::endl;
#endif

				// sort index vector by intensity of peaks (highest first)
				sort(indices_.begin(),indices_.end(),
						 reverseComparator(Internal::IntensityLess<Base>(*this))
						);

				// progress logger
				this->ff_->startProgress(0, indices_.size() , "FeatureFinder");

				current_peak_ = indices_.begin();

				initialized_ = true;
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

