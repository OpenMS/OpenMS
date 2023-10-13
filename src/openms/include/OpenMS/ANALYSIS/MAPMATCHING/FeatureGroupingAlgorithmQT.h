// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

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

