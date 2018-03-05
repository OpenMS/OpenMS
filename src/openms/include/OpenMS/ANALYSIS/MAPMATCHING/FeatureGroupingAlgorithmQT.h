// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMQT_H
#define OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMQT_H

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>

namespace OpenMS
{
  /**
      @brief A feature grouping algorithm for unlabeled data.

      The algorithm takes a number of feature or consensus maps and searches for corresponding (consensus) features across different maps. The maps have to be aligned (i.e. retention time distortions corrected, e.g. using one of the @ref MapAlignmentAlgorithm "MapAlignmentAlgorithms"), but small deviations are tolerated.

      This particular algorithm accumulates the features from all input maps, then applies a variant of QT clustering to find groups of corresponding features. For more details, see QTClusterFinder.

    @htmlinclude OpenMS_FeatureGroupingAlgorithmQT.parameters

      @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI FeatureGroupingAlgorithmQT :
    public FeatureGroupingAlgorithm
  {
public:
    /// Default constructor
    FeatureGroupingAlgorithmQT();

    /// Destructor
    ~FeatureGroupingAlgorithmQT() override;

    /**
        @brief Applies the algorithm to feature maps

        @pre The data ranges of the input maps have to be up-to-date (use FeatureMap::updateRanges).

        @exception IllegalArgument is thrown if less than two input maps are given.
    */
    void group(const std::vector<FeatureMap>& maps, ConsensusMap& out) override;

    /**
        @brief Applies the algorithm to consensus maps

         @pre The data ranges of the input maps have to be up-to-date (use ConsensusMap::updateRanges).

        @exception IllegalArgument is thrown if less than two input maps are given.
    */
    void group(const std::vector<ConsensusMap>& maps, 
                       ConsensusMap& out) override;

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
    FeatureGroupingAlgorithmQT(const FeatureGroupingAlgorithmQT&);

    /// Assignment operator intentionally not implemented -> private
    FeatureGroupingAlgorithmQT& operator=(const FeatureGroupingAlgorithmQT&);

    /**
        @brief Applies the algorithm to feature or consensus maps

        @pre The data ranges of the input maps have to be up-to-date (use MapType::updateRanges).

        @exception IllegalArgument is thrown if less than two input maps are given.
    */
    template <typename MapType>
    void group_(const std::vector<MapType>& maps, ConsensusMap& out);
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_MAPMATCHING_FEATUREGROUPINGALGORITHMQT_H
