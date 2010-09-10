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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMQT_H
#define OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMQT_H

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

namespace OpenMS
{
	/**
		@brief A map feature grouping algorithm for unlabeled data.

		It takes many maps and searches for corresponding features.
		The corresponding features must be aligned, but may have small position deviations.
		The input maps can be either of type ConsensusMap or FeatureMap.

	  @htmlinclude OpenMS_FeatureGroupingAlgorithmQT.parameters

		@ingroup FeatureGrouping
	*/
	class OPENMS_DLLAPI FeatureGroupingAlgorithmQT
	 : public FeatureGroupingAlgorithm
	{
		public:
			/// Default constructor
			FeatureGroupingAlgorithmQT();

			/// Destructor
			virtual ~FeatureGroupingAlgorithmQT();

			/**
				@brief Applies the algorithm

				@exception IllegalArgument is thrown if less than two input maps are given.
			*/
			virtual void group(const std::vector< FeatureMap<> >& maps, ConsensusMap& out);

			virtual void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out);

			/// Creates a new instance of this class (for Factory)
			static FeatureGroupingAlgorithm* create()
			{
				return new FeatureGroupingAlgorithmQT();
			}

			/// Returns the product name (for the Factory)
			static String getProductName()
			{
				return "unlabeled_qt";
			}

		private:

			/// Copy constructor intentionally not implemented -> private
			FeatureGroupingAlgorithmQT(const FeatureGroupingAlgorithmQT& );
			/// Assignment operator intentionally not implemented -> private
			FeatureGroupingAlgorithmQT& operator=(const FeatureGroupingAlgorithmQT& );

	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMQT_H
