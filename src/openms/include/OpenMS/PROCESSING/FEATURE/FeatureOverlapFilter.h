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
  class OPENMS_DLLAPI FeatureOverlapFilter
  {
    public:   
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
  };

}


