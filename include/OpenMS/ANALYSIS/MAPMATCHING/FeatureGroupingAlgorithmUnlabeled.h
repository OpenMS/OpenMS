// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMUNLABELED_H
#define OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMUNLABELED_H

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

namespace OpenMS
{
	/**
		@brief A map feature grouping algorithm for unlabeled data.
		
		It takes many maps and searches for corresponding features.
	*/
	class FeatureGroupingAlgorithmUnlabeled
	 : public FeatureGroupingAlgorithm
	{
		public:
			/// Default constructor
			FeatureGroupingAlgorithmUnlabeled();

			/// Destructor
			virtual ~FeatureGroupingAlgorithmUnlabeled();
			
			//Docu in base class
			virtual void group(const std::vector< FeatureMap<> >& maps, ConsensusMap<>& out);

			///Creates a new instance of this class (for Factory)
			static FeatureGroupingAlgorithm* create()
			{
				return new FeatureGroupingAlgorithmUnlabeled();
			}
			
			///Returns the product name (for the Factory)
			static String getProductName()
			{
				return "unlabeled";
			}
			
		private:

			///Copy constructor is not implemented -> private
			FeatureGroupingAlgorithmUnlabeled(const FeatureGroupingAlgorithmUnlabeled& );
			///Assignment operator is not implemented -> private
			FeatureGroupingAlgorithmUnlabeled& operator=(const FeatureGroupingAlgorithmUnlabeled& );
			
	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMUNLABELED_H
