// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMUNLABELED_H
#define OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMUNLABELED_H

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

namespace OpenMS
{
	/**
		@brief A map feature grouping algorithm for unlabeled data.
		
		It takes many maps and searches for corresponding features.
		The corresponding features must be aligned, but may have small position deviations.

	  @htmlinclude OpenMS_FeatureGroupingAlgorithmUnlabeled.parameters
	  
		@ingroup FeatureGrouping
	*/
	class OPENMS_DLLAPI FeatureGroupingAlgorithmUnlabeled
	 : public FeatureGroupingAlgorithm
	{
		public:
			/// Default constructor
			FeatureGroupingAlgorithmUnlabeled();

			/// Destructor
			virtual ~FeatureGroupingAlgorithmUnlabeled();
			
			/**
				@brief Applies the algorithm
				
				@exception IllegalArgument is thrown if less than two input maps are given.
			*/
			virtual void group(const std::vector< FeatureMap<> >& maps, ConsensusMap& out);

			/// Creates a new instance of this class (for Factory)
			static FeatureGroupingAlgorithm* create()
			{
				return new FeatureGroupingAlgorithmUnlabeled();
			}
			
			/// Returns the product name (for the Factory)
			static String getProductName()
			{
				return "unlabeled";
			}
			
		private:

			/// Copy constructor intentionally not implemented -> private
			FeatureGroupingAlgorithmUnlabeled(const FeatureGroupingAlgorithmUnlabeled& );
			/// Assignment operator intentionally not implemented -> private
			FeatureGroupingAlgorithmUnlabeled& operator=(const FeatureGroupingAlgorithmUnlabeled& );
			
	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMUNLABELED_H
