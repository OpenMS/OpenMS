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
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer, Oliver Kohlbacher $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHM_H
#define OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHM_H

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

	/**
	 * @brief Algorithms of ConsensusMapNormalizer
	 *
	 */
	class OPENMS_DLLAPI ConsensusMapNormalizerAlgorithm
	{
	private:
    /// copy constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithm(const ConsensusMapNormalizerAlgorithm &copyin);

    /// assignment operator is not implemented -> private
    ConsensusMapNormalizerAlgorithm& operator = (const ConsensusMapNormalizerAlgorithm &rhs);

  public:
		/// default constructor is not implemented -> private
		ConsensusMapNormalizerAlgorithm();

		/// destructor is not implemented -> private
		virtual ~ConsensusMapNormalizerAlgorithm();

    /**
		 * @brief determines the ratio of all maps to the map with the most features
		 * @param map ConsensusMap
		 * @param ratio_threshold threshold for the ratio
		 */
		static std::vector<double> computeCorrelation(const ConsensusMap& map, const double& ratio_threshold);
	
		/**
		 * @brief applies the given ratio to the maps of the consensusMap
		 * @param map ConsensusMap
		 * @param ratios ratios for the normalization
		 */
		static void normalizeMaps(ConsensusMap& map, const std::vector<double>& ratios);
	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHM_H
