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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHMQUANTILE_H
#define OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHMQUANTILE_H

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

	/**
	 * @brief Algorithms of ConsensusMapNormalizer
	 *
	 */
  class OPENMS_DLLAPI ConsensusMapNormalizerAlgorithmQuantile
	{
	private:
    /// copy constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmQuantile(const ConsensusMapNormalizerAlgorithmQuantile &copyin);

    /// assignment operator is not implemented -> private
    ConsensusMapNormalizerAlgorithmQuantile& operator = (const ConsensusMapNormalizerAlgorithmQuantile &rhs);

  public:
    /// default constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmQuantile();

    /// destructor is not implemented -> private
    virtual ~ConsensusMapNormalizerAlgorithmQuantile();

    /**
     * @brief normalizes the maps of the consensusMap
     * @param map ConsensusMap
     */
    static void normalizeMaps(ConsensusMap& map);

    /**
     * @brief resamples data_in and writes the results to data_out
     * @param data_in the data to be resampled
     * @param data_out the results are written to this vector
     * @param n_resampling_points the number of points to resample from data_in
     */
    static void resample(const std::vector<double>& data_in, std::vector<double>& data_out, UInt n_resampling_points);

    /**
     * @brief extracts the intensities of the features of the different maps
     * @param map ConsensusMap
     * @param out_intensities resulting data, contains the feature intensities for each map of the consensus map
     */
    static void extractIntensityVectors(const ConsensusMap& map, std::vector<std::vector<double> >& out_intensities);

    /**
     * @brief writes the intensity values in feature_ints to the corresponding features in map
     * @param feature_ints contains the new feature intensities for each map of the consensus map
     * @param map ConsensusMap the map to be updated
     */
    static void setNormalizedIntensityValues(const std::vector<std::vector<double> >& feature_ints, ConsensusMap& map);
	};

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHMQUANTILE_H
