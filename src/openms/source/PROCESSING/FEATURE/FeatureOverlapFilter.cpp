// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h>

#include <Quadtree.h>
#include <Box.h>
#include <unordered_set>

#include <cmath>

namespace OpenMS
{
  /// Boundaries for a mass trace in a feature
  struct MassTraceBounds
  {
    Size sub_index;
    double rt_min, rt_max, mz_min, mz_max;
  };

  /// Boundaries for all mass traces per feature
  using FeatureBoundsMap = std::map<UInt64, std::vector<MassTraceBounds>>;

  /// Get bounding boxes for all mass traces in all features of a feature map
  FeatureBoundsMap getFeatureBounds(const FeatureMap& features)
  {
    FeatureBoundsMap feature_bounds;
    for (const auto& feat : features)
    {
      for (Size i = 0; i < feat.getSubordinates().size(); ++i)
      {
        MassTraceBounds mtb;
        mtb.sub_index = i;
        const ConvexHull2D::PointArrayType& points =
        feat.getConvexHulls()[i].getHullPoints();
        mtb.mz_min = points.front().getY();
        mtb.mz_max = points.back().getY();
        const Feature& sub = feat.getSubordinates()[i];
        // convex hulls should be written out by "MRMFeatureFinderScoring" (see
        // parameter "write_convex_hull"):
        if (sub.getConvexHulls().empty())
        {
          String error = "convex hulls for mass traces missing";
          throw Exception::MissingInformation(__FILE__, __LINE__,
                                            OPENMS_PRETTY_FUNCTION, error);
        }
        const ConvexHull2D& hull = sub.getConvexHulls()[0];
        // find beginning of mass trace (non-zero intensity):
        if (hull.getHullPoints().empty())
        {
          continue;
        }
        double rt_min = hull.getHullPoints().back().getX();
        for (auto p_it = hull.getHullPoints().begin(); p_it != hull.getHullPoints().end(); ++p_it)
        {
          if (p_it->getY() > 0)
          {
            rt_min = p_it->getX();
            break;
          }
        }
        // find end of mass trace (non-zero intensity):
        double rt_max = hull.getHullPoints().front().getX();
        for (auto p_it =
             hull.getHullPoints().rbegin(); p_it !=
             hull.getHullPoints().rend(); ++p_it)
        {
          if (p_it->getX() < rt_min)
          {
            break;
          }
          if (p_it->getY() > 0)
          {
            rt_max = p_it->getX();
            break;
          }
        }
        if (rt_min > rt_max)
        {
          continue; // no peak -> skip
        }
        mtb.rt_min = rt_min;
        mtb.rt_max = rt_max;
        feature_bounds[feat.getUniqueId()].push_back(std::move(mtb));
      }
    }
    return feature_bounds;
  }

  /// Check if two sets of mass trace boundaries overlap
  bool hasOverlappingBounds(const std::vector<MassTraceBounds>& mtb1, const std::vector<MassTraceBounds>& mtb2)
  {
    for (const MassTraceBounds& mt1 : mtb1)
    {
      for (const MassTraceBounds& mt2 : mtb2)
      {
        if (!((mt1.rt_max < mt2.rt_min) ||
              (mt1.rt_min > mt2.rt_max) ||
              (mt1.mz_max < mt2.mz_min) ||
              (mt1.mz_min > mt2.mz_max)))
        {
          return true;
        }
      }
    }
    return false;
  }

  bool tracesOverlap(const Feature& a, const Feature& b, const FeatureBoundsMap& feature_bounds)
  {
    auto fbm_it1 = feature_bounds.find(a.getUniqueId());
    auto fbm_it2 = feature_bounds.find(b.getUniqueId());
    return hasOverlappingBounds(fbm_it1->second, fbm_it2->second);
  }

  void FeatureOverlapFilter::filter(FeatureMap& fmap, 
    std::function<bool(const Feature&, const Feature&)> FeatureComparator, 
    std::function<bool(Feature&, Feature&)> FeatureOverlapCallback,
    bool check_overlap_at_trace_level)
  {
    // Delegate to the new overload with appropriate mode
    FeatureOverlapMode mode = check_overlap_at_trace_level ? FeatureOverlapMode::TRACE_LEVEL : FeatureOverlapMode::CONVEX_HULL;
    CentroidTolerances tolerances; // Use default values
    filter(fmap, FeatureComparator, FeatureOverlapCallback, mode, tolerances);
  }
  
  void FeatureOverlapFilter::filter(FeatureMap& fmap,
    std::function<bool(const Feature&, const Feature&)> FeatureComparator,
    std::function<bool(Feature&, Feature&)> FeatureOverlapCallback,
    FeatureOverlapMode mode,
    const CentroidTolerances& tolerances)
  {
    fmap.updateRanges();
    // Sort all features according to the comparator. After the sort, the "smallest" == best feature will be the first entry we will start processing with...
    std::stable_sort(fmap.begin(), fmap.end(), FeatureComparator);

    // Define getBox function based on mode
    std::function<quadtree::Box<float>(const Feature*)> getBox;
    
    if (mode == FeatureOverlapMode::CENTROID_BASED)
    {
      // For centroid-based mode, create tolerance boxes around centroids
      getBox = [&tolerances](const Feature* f)
      {
        float rt = f->getRT();
        float mz = f->getMZ();
        return quadtree::Box<float>(
          mz - tolerances.mz_tolerance, 
          rt - tolerances.rt_tolerance, 
          2 * tolerances.mz_tolerance, 
          2 * tolerances.rt_tolerance
        );
      };
    }
    else
    {
      // For convex hull/trace modes, use full convex hull bounding boxes
      getBox = [](const Feature* f)
      {
        const auto& bb = f->getConvexHull().getBoundingBox();
        return quadtree::Box<float>(bb.minY(), bb.minX(), bb.maxY()-bb.minY(), bb.maxX()-bb.minX());
      };
    }

    float minMZ = fmap.getMinMZ();
    float maxMZ = fmap.getMaxMZ();
    float minRT = fmap.getMinRT();
    float maxRT = fmap.getMaxRT();

    // Expand boundaries for centroid mode to accommodate tolerance boxes
    if (mode == FeatureOverlapMode::CENTROID_BASED)
    {
      minMZ -= tolerances.mz_tolerance;
      maxMZ += tolerances.mz_tolerance;
      minRT -= tolerances.rt_tolerance;
      maxRT += tolerances.rt_tolerance;
    }

    // Build quadtree with all features
    quadtree::Box<float> fullExp(minMZ-1, minRT-1, maxMZ-minMZ+2, maxRT-minRT+2);
    auto quadtree = quadtree::Quadtree<Feature*, decltype(getBox)>(fullExp, getBox);
    for (auto& f : fmap)
    {
      quadtree.add(&f);
    }        

    // If we check for overlapping traces we need a faster lookup structure
    FeatureBoundsMap fbm;
    if (mode == FeatureOverlapMode::TRACE_LEVEL)
    {
      fbm = getFeatureBounds(fmap);
    }

    std::unordered_set<Size> removed_uids;
    for (auto& f : fmap)
    {
      if (removed_uids.count(f.getUniqueId()) == 0)
      {
        for (auto& overlap : quadtree.query(getBox(&f)))
        {
          if ((overlap != &f))
          {
            bool is_true_overlap = true;
            
            if (mode == FeatureOverlapMode::CENTROID_BASED)
            {
              // Check charge requirement
              if (tolerances.require_same_charge && f.getCharge() != overlap->getCharge())
              {
                is_true_overlap = false;
              }
              else
              {
                // Check exact centroid distances within tolerance
                double rt_diff = std::abs(f.getRT() - overlap->getRT());
                double mz_diff = std::abs(f.getMZ() - overlap->getMZ());
                is_true_overlap = (rt_diff <= tolerances.rt_tolerance && mz_diff <= tolerances.mz_tolerance);
              }
            }
            else if (mode == FeatureOverlapMode::TRACE_LEVEL)
            {            
              is_true_overlap = tracesOverlap(f, *overlap, fbm);
            }
            // For CONVEX_HULL mode, is_true_overlap remains true (quadtree query already handles overlap)

            if (is_true_overlap)
            {
              // callback allows to e.g., transfer information from the to-be-removed feature to the representative feature
              // if the callback returns false, overlap will not be removed (at least not because of an overlap with f)
              if (FeatureOverlapCallback(f, *overlap)) 
              {
                removed_uids.insert(overlap->getUniqueId());
              }                            
            }
          }
        }
      }
    }

    const auto filtered = [&removed_uids](const Feature& f)
    {
      return removed_uids.count(f.getUniqueId()) == 1;
    };
    fmap.erase(std::remove_if(fmap.begin(), fmap.end(), filtered), fmap.end());
  }

}
