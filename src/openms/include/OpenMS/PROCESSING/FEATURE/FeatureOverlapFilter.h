// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{
  /// Enum to specify the overlap detection mode
  enum class FeatureOverlapMode
  {
    CONVEX_HULL,     ///< Use convex hull bounding boxes (default)
    TRACE_LEVEL,     ///< Check overlap at trace level
    CENTROID_BASED   ///< Check overlap based on centroid distances
  };

  /// Enum to specify how intensities are combined when merging features
  enum class MergeIntensityMode
  {
    SUM,  ///< Sum intensities of merged features (default)
    MAX   ///< Keep maximum intensity
  };

  /// Structure to hold centroid-based overlap tolerances
  struct CentroidTolerances
  {
    double rt_tolerance = 5.0;        ///< Maximum RT difference in seconds
    double mz_tolerance = 0.05;       ///< Maximum m/z difference in Da
    bool require_same_charge = true;  ///< Whether to require identical charge states
    bool require_same_im = false;     ///< Whether to require identical FAIMS CV (or both missing)
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

    /**
        @brief Create a callback function for merging overlapping features.

        Creates a callback suitable for use with filter() that merges feature information
        when features overlap. The callback:
        - Combines intensities according to intensity_mode (SUM or MAX)
        - Optionally stores merge tracking information as meta values:
          - "merged_centroid_rts": vector of RT positions from all merged features
          - "merged_centroid_mzs": vector of m/z positions from all merged features
          - "merged_centroid_IMs": vector of FAIMS CV values (only if FAIMS_CV present on features)
          - "FAIMS_merge_count": count of merged FAIMS CV values

        This callback is designed to work with features that may or may not have FAIMS CV
        annotations. Features without FAIMS_CV will simply not contribute to merged_centroid_IMs.

        @param intensity_mode How to combine intensities: SUM adds all intensities,
                              MAX keeps the highest intensity (default: SUM)
        @param write_meta_values If true, write merge tracking meta values to the surviving
                                 feature for debugging/analysis (default: true)
        @return A callback function suitable for use with filter()
    */
    static std::function<bool(Feature&, Feature&)> createFAIMSMergeCallback(
      MergeIntensityMode intensity_mode = MergeIntensityMode::SUM,
      bool write_meta_values = true);

    /**
        @brief Merge overlapping features based on centroid distances.

        Identifies features whose centroids are within the specified RT and m/z tolerances
        and merges them into a single representative feature. This is useful for:
        - Combining features detected at the same location across different FAIMS CV values
        - Deduplicating features from different detection methods
        - Consolidating split features

        The feature with the highest intensity is kept as the representative. Depending on
        intensity_mode, intensities are either summed or the maximum is kept.

        When write_meta_values is true, the following meta values are stored on merged features:
        - "merged_centroid_rts": RT positions of all features that were merged
        - "merged_centroid_mzs": m/z positions of all features that were merged
        - "merged_centroid_IMs": FAIMS CV values (if present on the original features)
        - "FAIMS_merge_count": number of FAIMS CV values that were merged

        @param feature_map The feature map to process (modified in place)
        @param max_rt_diff Maximum RT difference in seconds for considering features as
                          overlapping (default: 5.0)
        @param max_mz_diff Maximum m/z difference in Da for considering features as
                          overlapping (default: 0.05)
        @param require_same_charge If true, only merge features with identical charge states.
                                   If false, features with different charges can be merged
                                   (default: true)
        @param require_same_im If true, only merge features with identical FAIMS CV values.
                               Features without FAIMS_CV are treated as a separate group and
                               can only merge with other features lacking FAIMS_CV (default: false)
        @param intensity_mode How to combine intensities: SUM adds all intensities,
                              MAX keeps the highest intensity (default: SUM)
        @param write_meta_values If true, write merge tracking meta values to merged features
                                 (default: true)
    */
    static void mergeOverlappingFeatures(FeatureMap& feature_map,
                                         double max_rt_diff = 5.0,
                                         double max_mz_diff = 0.05,
                                         bool require_same_charge = true,
                                         bool require_same_im = false,
                                         MergeIntensityMode intensity_mode = MergeIntensityMode::SUM,
                                         bool write_meta_values = true);
  };

}


