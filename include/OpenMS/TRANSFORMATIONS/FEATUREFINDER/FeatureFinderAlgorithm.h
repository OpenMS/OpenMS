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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

	// forward declaration
	class FeatureFinder;
	
	/// Summary of fitting results
	struct OPENMS_DLLAPI Summary
	{
		std::map<String,UInt> exception; //count exceptions
		UInt no_exceptions;
		std::map<String,UInt> mz_model; //count used mz models
		std::map<float,UInt> mz_stdev; //count used mz standard deviations
		std::vector<UInt> charge; //count used charges
		DoubleReal corr_mean, corr_max, corr_min; 	//boxplot for correlation
		
		/// Initial values
		Summary() :
			no_exceptions(0),
			corr_mean(0),
			corr_max(0),
			corr_min(1)
		{}
	
	};

	/**
		@brief Abstract base class for FeatureFinder algorithms

	*/
	template<class PeakType, class FeatureType> class FeatureFinderAlgorithm
		: public DefaultParamHandler
	{
		public:
			/// Input map type
			typedef MSExperiment<PeakType> MapType;
			/// Coordinate/Position type of peaks
			typedef typename MapType::CoordinateType CoordinateType;
			/// Intensity type of peaks
			typedef typename MapType::IntensityType IntensityType;
			/// Output feature type
			typedef FeatureMap<FeatureType> FeatureMapType;

			/// default constructor
			FeatureFinderAlgorithm() :
				DefaultParamHandler("FeatureFinderAlgorithm"),
				map_(0),
				features_(0),
				ff_(0)
			{
			}

			/// destructor
			virtual ~FeatureFinderAlgorithm()
			{
			}

			/// register all derived classes here (see FeatureFinderAlgorithm_impl.h)
			static void registerChildren();

			/// Main method that implements the actual algorithm
			virtual void run()=0;

			/**
				@brief Returns the default parameters. Reimplment
				
				Reimplment if you derive a class and have to incoopreate sub-algorithm default parameters.
			*/
			virtual Param getDefaultParameters() const
			{
				return this->defaults_;
			}

			/// Sets a reference to the calling FeatureFinder
			void setData(const MapType& map, FeatureMapType& features, FeatureFinder& ff)
			{
				map_ = &map;
				features_ = &features;
				ff_ = &ff;
			}
			
			/**
				@brief Sets a reference to the calling FeatureFinder
				
				@exception Exception::IllegalArgument is thrown if the algorithm does not support user-specified seed lists
			*/
			virtual void setSeeds(const FeatureMapType& seeds)
			{
				if (seeds.size()!=0)
				{
					throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"The used feature detection algorithm does not support user-specified seed lists!");
				}
			}
			
		protected:

			/// Input data pointer
			const MapType* map_;

			/// Output data pointer
			FeatureMapType* features_;

			/// Pointer to the calling FeatureFinder that is used to access the feature flags
			FeatureFinder* ff_;
		
		private:

			/// Not implemented
			FeatureFinderAlgorithm& operator=(const FeatureFinderAlgorithm&);

			/// Not implemented
			FeatureFinderAlgorithm(const FeatureFinderAlgorithm&);

	};
}

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H
