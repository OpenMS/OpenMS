// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h>
#include <OpenMS/CONCEPT/Constants.h>

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
              // Check FAIMS CV requirement
              else if (tolerances.require_same_im)
              {
                bool f_has_im = f.metaValueExists(Constants::UserParam::FAIMS_CV);
                bool overlap_has_im = overlap->metaValueExists(Constants::UserParam::FAIMS_CV);

                if (f_has_im != overlap_has_im)
                {
                  // One has FAIMS CV, the other doesn't - not same group
                  is_true_overlap = false;
                }
                else if (f_has_im && overlap_has_im)
                {
                  // Both have FAIMS CV - must match
                  double f_cv = f.getMetaValue(Constants::UserParam::FAIMS_CV);
                  double overlap_cv = overlap->getMetaValue(Constants::UserParam::FAIMS_CV);
                  if (f_cv != overlap_cv)
                  {
                    is_true_overlap = false;
                  }
                }
                // else: both don't have FAIMS CV - same group, continue to distance check

                if (is_true_overlap)
                {
                  // Check exact centroid distances within tolerance
                  double rt_diff = std::abs(f.getRT() - overlap->getRT());
                  double mz_diff = std::abs(f.getMZ() - overlap->getMZ());
                  is_true_overlap = (rt_diff <= tolerances.rt_tolerance && mz_diff <= tolerances.mz_tolerance);
                }
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

  std::function<bool(Feature&, Feature&)> FeatureOverlapFilter::createFAIMSMergeCallback(
    MergeIntensityMode intensity_mode,
    bool write_meta_values)
  {
    return [intensity_mode, write_meta_values](Feature& best_in_cluster, Feature& f) -> bool
    {
      double best_intensity = best_in_cluster.getIntensity();
      double f_intensity = f.getIntensity();

      if (write_meta_values)
      {
        // Collect centroid RT positions
        std::vector<double> merged_rts;
        if (best_in_cluster.metaValueExists("merged_centroid_rts"))
        {
          merged_rts = best_in_cluster.getMetaValue("merged_centroid_rts");
        }
        else
        {
          merged_rts.push_back(best_in_cluster.getRT());
        }
        merged_rts.push_back(f.getRT());
        best_in_cluster.setMetaValue("merged_centroid_rts", merged_rts);

        // Collect centroid m/z positions
        std::vector<double> merged_mzs;
        if (best_in_cluster.metaValueExists("merged_centroid_mzs"))
        {
          merged_mzs = best_in_cluster.getMetaValue("merged_centroid_mzs");
        }
        else
        {
          merged_mzs.push_back(best_in_cluster.getMZ());
        }
        merged_mzs.push_back(f.getMZ());
        best_in_cluster.setMetaValue("merged_centroid_mzs", merged_mzs);

        // Collect FAIMS CV values (only if present on features)
        std::vector<double> merged_ims;
        if (best_in_cluster.metaValueExists("merged_centroid_IMs"))
        {
          merged_ims = best_in_cluster.getMetaValue("merged_centroid_IMs");
        }
        else if (best_in_cluster.metaValueExists(Constants::UserParam::FAIMS_CV))
        {
          merged_ims.push_back(best_in_cluster.getMetaValue(Constants::UserParam::FAIMS_CV));
          best_in_cluster.removeMetaValue(Constants::UserParam::FAIMS_CV);
        }

        if (f.metaValueExists(Constants::UserParam::FAIMS_CV))
        {
          merged_ims.push_back(f.getMetaValue(Constants::UserParam::FAIMS_CV));
        }

        if (!merged_ims.empty())
        {
          best_in_cluster.setMetaValue("merged_centroid_IMs", merged_ims);
          best_in_cluster.setMetaValue("FAIMS_merge_count", static_cast<int>(merged_ims.size()));
        }
      }

      // Combine intensities according to mode
      double new_intensity = best_intensity + f_intensity; // default: SUM
      if (intensity_mode == MergeIntensityMode::MAX)
      {
        new_intensity = std::max(best_intensity, f_intensity);
      }
      best_in_cluster.setIntensity(new_intensity);

      return true; // Remove the overlapping feature
    };
  }

  void FeatureOverlapFilter::mergeOverlappingFeatures(FeatureMap& feature_map,
                                                      double max_rt_diff,
                                                      double max_mz_diff,
                                                      bool require_same_charge,
                                                      bool require_same_im,
                                                      MergeIntensityMode intensity_mode,
                                                      bool write_meta_values)
  {
    CentroidTolerances tolerances;
    tolerances.rt_tolerance = max_rt_diff;
    tolerances.mz_tolerance = max_mz_diff;
    tolerances.require_same_charge = require_same_charge;
    tolerances.require_same_im = require_same_im;

    // Use intensity-based comparator (higher intensity = better = "smaller" in sort order)
    auto intensity_comparator = [](const Feature& left, const Feature& right)
    {
      return left.getIntensity() > right.getIntensity();
    };

    filter(feature_map,
           intensity_comparator,
           createFAIMSMergeCallback(intensity_mode, write_meta_values),
           FeatureOverlapMode::CENTROID_BASED,
           tolerances);
  }

  void FeatureOverlapFilter::mergeFAIMSFeatures(FeatureMap& feature_map,
                                                double max_rt_diff,
                                                double max_mz_diff)
  {
    // Check if any features have FAIMS_CV - if not, nothing to do
    bool has_faims_features = false;
    for (const auto& f : feature_map)
    {
      if (f.metaValueExists(Constants::UserParam::FAIMS_CV))
      {
        has_faims_features = true;
        break;
      }
    }

    if (!has_faims_features)
    {
      return; // No FAIMS features, nothing to merge
    }

    // Separate features into FAIMS and non-FAIMS groups
    FeatureMap faims_features;
    FeatureMap non_faims_features;

    for (auto& f : feature_map)
    {
      if (f.metaValueExists(Constants::UserParam::FAIMS_CV))
      {
        faims_features.push_back(std::move(f));
      }
      else
      {
        non_faims_features.push_back(std::move(f));
      }
    }

    // Only merge the FAIMS features if we have more than one
    if (faims_features.size() > 1)
    {
      CentroidTolerances tolerances;
      tolerances.rt_tolerance = max_rt_diff;
      tolerances.mz_tolerance = max_mz_diff;
      tolerances.require_same_charge = true;
      tolerances.require_same_im = false; // We handle IM check in callback

      // Custom callback that only merges features with DIFFERENT FAIMS CV values
      auto merge_callback = [](Feature& best_in_cluster, Feature& f) -> bool
      {
        // Only merge if FAIMS CVs are DIFFERENT
        // (same CV features should not be merged - they are different analytes)
        double best_cv = best_in_cluster.getMetaValue(Constants::UserParam::FAIMS_CV);
        double f_cv = f.getMetaValue(Constants::UserParam::FAIMS_CV);

        if (best_cv == f_cv)
        {
          return false; // Don't merge features with same CV
        }

        // Merge features with different CVs - sum intensities
        double best_intensity = best_in_cluster.getIntensity();
        double f_intensity = f.getIntensity();

        // Collect centroid RT positions
        std::vector<double> merged_rts;
        if (best_in_cluster.metaValueExists("merged_centroid_rts"))
        {
          merged_rts = best_in_cluster.getMetaValue("merged_centroid_rts");
        }
        else
        {
          merged_rts.push_back(best_in_cluster.getRT());
        }
        merged_rts.push_back(f.getRT());
        best_in_cluster.setMetaValue("merged_centroid_rts", merged_rts);

        // Collect centroid m/z positions
        std::vector<double> merged_mzs;
        if (best_in_cluster.metaValueExists("merged_centroid_mzs"))
        {
          merged_mzs = best_in_cluster.getMetaValue("merged_centroid_mzs");
        }
        else
        {
          merged_mzs.push_back(best_in_cluster.getMZ());
        }
        merged_mzs.push_back(f.getMZ());
        best_in_cluster.setMetaValue("merged_centroid_mzs", merged_mzs);

        // Collect FAIMS CV values
        std::vector<double> merged_ims;
        if (best_in_cluster.metaValueExists("merged_centroid_IMs"))
        {
          merged_ims = best_in_cluster.getMetaValue("merged_centroid_IMs");
        }
        else
        {
          merged_ims.push_back(best_cv);
          best_in_cluster.removeMetaValue(Constants::UserParam::FAIMS_CV);
        }
        merged_ims.push_back(f_cv);
        best_in_cluster.setMetaValue("merged_centroid_IMs", merged_ims);
        best_in_cluster.setMetaValue("FAIMS_merge_count", static_cast<int>(merged_ims.size()));

        // Sum intensities
        best_in_cluster.setIntensity(best_intensity + f_intensity);

        return true; // Remove the merged feature
      };

      // Use intensity-based comparator
      auto intensity_comparator = [](const Feature& left, const Feature& right)
      {
        return left.getIntensity() > right.getIntensity();
      };

      filter(faims_features,
             intensity_comparator,
             merge_callback,
             FeatureOverlapMode::CENTROID_BASED,
             tolerances);
    }

    // Combine back: merged FAIMS features + untouched non-FAIMS features
    feature_map.clear();
    for (auto& f : faims_features)
    {
      feature_map.push_back(std::move(f));
    }
    for (auto& f : non_faims_features)
    {
      feature_map.push_back(std::move(f));
    }
  }

}
