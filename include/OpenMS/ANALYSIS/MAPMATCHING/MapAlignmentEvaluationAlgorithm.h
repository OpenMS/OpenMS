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
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTEVALUATIONALGORITHM_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTEVALUATIONALGORITHM_H

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

	/**
		@brief Base class for all Caap evaluation algorithms
		
		These algorithms evaluates alignment results against a ground truth.
	*/
	class OPENMS_DLLAPI MapAlignmentEvaluationAlgorithm
	{

		protected:
			typedef ConsensusFeature::HandleSetType::const_iterator HandleIterator;

		public:

			/// Default constructor
			MapAlignmentEvaluationAlgorithm();

			/// Destructor
			virtual ~MapAlignmentEvaluationAlgorithm();

			
			///Applies the algorithm. The input consensus map is compared to the ground truth.
			virtual void evaluate(const ConsensusMap & conensus_map_in, const ConsensusMap & consensus_map_gt, const DoubleReal& rt_dev, const DoubleReal& mz_dev, const Peak2D::IntensityType& int_dev, DoubleReal & out)=0;

			///Decides if two features are the same, based on maximum allowed deviations for retention time, m/z and intensity.
			bool isSameHandle(const FeatureHandle & lhs, const FeatureHandle & rhs, const DoubleReal& rt_dev, const DoubleReal& mz_dev, const Peak2D::IntensityType& int_dev);

			/// Register all derived classes in this method
			static void registerChildren();
		
		private:
			///Copy constructor is not implemented -> private
			MapAlignmentEvaluationAlgorithm(const MapAlignmentEvaluationAlgorithm& );
			///Assignment operator is not implemented -> private
			MapAlignmentEvaluationAlgorithm& operator=(const MapAlignmentEvaluationAlgorithm& );
			
	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTEVALUATIONALGORITHM_H
