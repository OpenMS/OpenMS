// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  /// Enum to specify the overlap detection mode
  enum class FeatureOverlapMode
  {
    CONVEX_HULL,     ///< Use convex hull bounding boxes (default)
    TRACE_LEVEL,     ///< Check overlap at trace level
    CENTROID_BASED   ///< Check overlap based on centroid distances
  };
  
  /// Structure to hold centroid-based overlap tolerances
  struct CentroidTolerances
  {
    double rt_tolerance = 5.0;        ///< Maximum RT difference in seconds
    double mz_tolerance = 0.05;       ///< Maximum m/z difference in Da
    bool require_same_charge = true;  ///< Whether to require identical charge states
  };

  class OPENMS_DLLAPI FeatureOverlapFilter
  {
    public:
    
    /// Enum to specify the overlap detection mode (alias for backward compatibility)
    using OverlapMode = FeatureOverlapMode;
    
    /*
        @brief Filter overlapping features using a spatial datastructure (quadtree). 
               Retains only the best feature in each cluster of overlapping features.

        @param FeatureComparator must implement the concept of a less comparator.
               If several features overlap, the feature that evaluates as "smallest" is considered the best (according to the passed comparator) and is kept.
               The other overlapping features are removed and FeatureOverlapCallback evaluated on them.
               Default: overall feature quality.

        @param FeatureOverlapCallback(best_in_cluster, f) is called if a feature f overlaps with a feature best_in_cluster.
               FeatureOverlapCallback provides a customization point to e.g.:
              - transfer information from the soon-to-be-removed feature f over to the best_in_cluster feature
              - gather overlap statistics
              - help in debugging
              - etc.
              in form of a callable.
              If the FeatureOverlapCallback returns false, the overlapping feature will be treated as not overlapping with best_in_cluster (and not removed).
              Default: function that just returns true.

        @ingroup Datareduction
    */
    static void filter(FeatureMap& fmap, 
      std::function<bool(const Feature&, const Feature&)> FeatureComparator = [](const Feature& left, const Feature& right){ return left.getOverallQuality() > right.getOverallQuality(); },
      std::function<bool(Feature&, Feature&)> FeatureOverlapCallback = [](Feature&, Feature&){ return true; },
      bool check_overlap_at_trace_level = true);
      
    /*
        @brief Filter overlapping features with configurable overlap detection mode.
               Extended version that allows choosing between different overlap detection strategies.

        @param fmap The feature map to filter
        @param FeatureComparator Comparator to determine the best feature in overlapping clusters
        @param FeatureOverlapCallback Callback function called when features overlap
        @param mode The overlap detection mode to use
        @param tolerances Tolerances for centroid-based overlap detection (only used when mode == CENTROID_BASED)

        @ingroup Datareduction
    */
    static void filter(FeatureMap& fmap,
      std::function<bool(const Feature&, const Feature&)> FeatureComparator,
      std::function<bool(Feature&, Feature&)> FeatureOverlapCallback,
      FeatureOverlapMode mode,
      const CentroidTolerances& tolerances = CentroidTolerances());
  };

}


