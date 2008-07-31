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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHM_H

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

namespace OpenMS
{

	// forward declaration
	class FeatureFinder;
	
	/**
		@brief Abstract base class for FeatureFinder algorithms

	*/
	template<class PeakType, class FeatureType> class FeatureFinderAlgorithm : 
		public FactoryProduct
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
				FactoryProduct("FeatureFinderAlgorithm"),
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
