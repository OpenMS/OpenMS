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
// $Maintainer: Katharina Albers, Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTEVALUATIONALGORITHM_H
#define OPENMS_ANALYSIS_MAPMATCHING_MAPALIGNMENTEVALUATIONALGORITHM_H

#include <OpenMS/CONCEPT/FactoryProduct.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

	/**
		@brief Base class for all Caap evaluation algorithms
		
		These algorithms evaluates alignment results against a ground truth.
	*/
	class OPENMS_DLLAPI MapAlignmentEvaluationAlgorithm
	 : public FactoryProduct
	{

		protected:
			typedef ConsensusFeature::HandleSetType::const_iterator HandleIterator; //geht nicht private! fehler:/home/bude/albers/RAID/cmakeOpenMS/include/OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h:46: error: 'typedef struct std::_Rb_tree_const_iterator<OpenMS::FeatureHandle> OpenMS::MapAlignmentEvaluationAlgorithm::HandleIterator' is private

		public:

			/// Default constructor
			MapAlignmentEvaluationAlgorithm();

			/// Destructor
			virtual ~MapAlignmentEvaluationAlgorithm();

			
			///Applies the algorithm. The input consensus map is compared to the ground truth.
			virtual void evaluate(const ConsensusMap & conensus_map_in, const ConsensusMap & consensus_map_gt, DoubleReal & out)=0;

			bool isSameHandle(const FeatureHandle & lhs, const FeatureHandle & rhs);

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
