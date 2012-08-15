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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHMMEDIAN_H
#define OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHMMEDIAN_H

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

	/**
	 * @brief Algorithms of ConsensusMapNormalizer
	 *
	 */
  class OPENMS_DLLAPI ConsensusMapNormalizerAlgorithmMedian
	{
	private:
    /// copy constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmMedian(const ConsensusMapNormalizerAlgorithmMedian &copyin);

    /// assignment operator is not implemented -> private
    ConsensusMapNormalizerAlgorithmMedian& operator = (const ConsensusMapNormalizerAlgorithmMedian &rhs);

  public:
		/// default constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmMedian();

		/// destructor is not implemented -> private
    virtual ~ConsensusMapNormalizerAlgorithmMedian();

		/**
     * @brief normalizes the maps of the consensusMap
		 * @param map ConsensusMap
		 */
    static void normalizeMaps(ConsensusMap& map);

    /**
     * @brief computes the normalization factors for all maps
     * @param map ConsensusMap
     */
    static std::vector<double> computeNormalizationFactors(const ConsensusMap& map);
	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHMMEDIAN_H
