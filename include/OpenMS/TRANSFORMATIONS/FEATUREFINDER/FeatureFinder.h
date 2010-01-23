// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderDefs.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

	/**@brief The main feature finder class.

		- Stores the flags for (indices of) data points ("used", "unused")
		- The algorithm itself is a factory product (derived from FeatureFinderAlgorithm)
		- The main method is run(), which is a template so that we can deal with different types of input and output
		- The run() method takes five arguments: algorithm_name, input_map, output, parameters, seeds
		.	

		@ingroup FeatureFinder
	*/
	class OPENMS_DLLAPI FeatureFinder :
		public ProgressLogger,
	 	public FeatureFinderDefs
	{

		public:
			/// Default constructor.
			FeatureFinder();

			/// Destructor
			virtual ~FeatureFinder();

			/**
				@brief Executes the FeatureFinder using the given algorithm

				There are several constraints for the @p input_map.  They are tested before
				the algorithm starts.  It must only contain MS 1 level scans and you
				have to call updateRanges() before passing it to this method.

				@param algorithm_name Name of the feature finding algorithm to use
				@param input_map Input peak map
				@param features Output feature map
				@param param Algorithm parameters
				@param seeds List of seeds to use

				Implemented in FeatureFinder_impl.h
			*/
			template<class PeakType, class FeatureType>
			void run(const String& algorithm_name, MSExperiment<PeakType> const & input_map, FeatureMap<FeatureType> & features, const Param& param, const FeatureMap<FeatureType>& seeds);

			/// Returns a non-mutable reference to a peak flag
			const Flag& getPeakFlag(const IndexPair& index) const
			{
				return flags_[index.first][index.second];
			}

			/// Returns mutable reference to a peak flag
			Flag& getPeakFlag(const IndexPair& index) 
			{ 
				return flags_[index.first][index.second];
			}

			/// Returns the default parameters for the algorithm with name @p algorithm_name
			Param getParameters(const String& algorithm_name) const;

		protected:

			/// Container for flags attached to input data
			std::vector< std::vector<Flag> > flags_;

	}; // class FeatureFinder

} // namespace OpenMS

#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDER_H
