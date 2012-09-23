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
// $Maintainer: Lars Nilse $
// $Authors: Hendrik Brauer, Oliver Kohlbacher, Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHMTHRESHOLD_H
#define OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHMTHRESHOLD_H

#include <OpenMS/KERNEL/ConsensusMap.h>

namespace OpenMS
{

  /**
   * @brief Algorithms of ConsensusMapNormalizer
   *
   */
  class OPENMS_DLLAPI ConsensusMapNormalizerAlgorithmThreshold
  {
private:
    /// copy constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmThreshold(const ConsensusMapNormalizerAlgorithmThreshold & copyin);

    /// assignment operator is not implemented -> private
    ConsensusMapNormalizerAlgorithmThreshold & operator=(const ConsensusMapNormalizerAlgorithmThreshold & rhs);

public:
    /// default constructor is not implemented -> private
    ConsensusMapNormalizerAlgorithmThreshold();

    /// destructor is not implemented -> private
    virtual ~ConsensusMapNormalizerAlgorithmThreshold();

    /**
         * @brief determines the ratio of all maps to the map with the most features
         * @param map ConsensusMap
         * @param ratio_threshold threshold for the ratio
         */
    static std::vector<double> computeCorrelation(const ConsensusMap & map, const double & ratio_threshold);

    /**
     * @brief applies the given ratio to the maps of the consensusMap
     * @param map ConsensusMap
     * @param ratios ratios for the normalization
     */
    static void normalizeMaps(ConsensusMap & map, const std::vector<double> & ratios);
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_CONSENSUSMAPNORMALIZERALGORITHMTHRESHOLD_H
