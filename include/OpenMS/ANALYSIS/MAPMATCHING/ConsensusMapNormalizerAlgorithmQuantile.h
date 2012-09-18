// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
