// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <set>

namespace OpenMS
{
  class BaseFeature;
  class AASequence;

  /**
   * @brief Representation of a feature in a hash grid.
   *
   * A GridFeature can be stored in a HashGrid and points to a BaseFeature (Feature or ConsensusFeature). Used for QT feature grouping (see QTClusterFinder).
   */
  class OPENMS_DLLAPI GridFeature
  {
private:
    /// Reference to the contained feature
    const BaseFeature& feature_;

    /// Index of the feature map or consensus map
    Size map_index_;

    /// Index of the feature in the map
    Size feature_index_;

    /// Set of peptide sequences annotated to the feature
    std::set<AASequence> annotations_;

public:
    /**
     * @brief Detailed constructor
     * @param feature Reference to the contained feature
     * @param map_index Index of the feature map or consensus map
     * @param feature_index Index of the feature in the map
     */
    GridFeature(const BaseFeature& feature, Size map_index, Size feature_index);

    /// Returns the feature
    const BaseFeature& getFeature() const;

    /// Destructor
    virtual ~GridFeature();

    /// Returns the map index
    Size getMapIndex() const;

    /// Returns the feature index
    Size getFeatureIndex() const;

    /// Returns the ID of the GridFeature (same as the feature index)
    Int getID() const;

    /// Returns the set of peptide sequences annotated to the cluster center
    const std::set<AASequence>& getAnnotations() const;

    /// Returns the feature RT
    double getRT() const;

    /// Returns the feature m/z
    double getMZ() const;
  };
}

