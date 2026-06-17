// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>
#include <OpenMS/DATASTRUCTURES/KDTree.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureNode.h>

namespace OpenMS
{

/**
  @brief Pools features from multiple maps into a single 2D kd-tree on (RT, m/z) for fast range queries.

  The class stores raw pointers to @ref BaseFeature instances supplied by the caller (no copy is made),
  along with the map of origin for each feature. A balanced kd-tree over the (RT, m/z) plane is built so
  that subsequent neighbourhood searches across maps run in @f$\mathcal{O}(\log n + k)@f$ instead of a
  linear scan over all features.

  Typical use:
    - KD-tree based map alignment (see @ref MapAlignmentAlgorithmKD), where pairwise RT correspondences
      between maps are collected from features that fall within an RT/m/z tolerance window.
    - KD-tree based feature grouping (see @ref FeatureGroupingAlgorithmKD), where consensus features
      are built by clustering compatible features from different maps.

  RT values are stored separately from the wrapped features so that an alignment transformation can
  later be applied via applyTransformations() without modifying the original features. m/z, intensity,
  and charge are always read through the feature pointer and reflect the current state of the
  underlying feature.

  Ownership: the caller retains ownership of the wrapped @ref BaseFeature objects. The kd-tree only
  stores pointers, so the wrapped features must outlive this object.

  @note This class derives from @ref DefaultParamHandler but disables the parameter check
        (@c check_defaults_ is set to false), so subclasses or callers may pass arbitrary @ref Param
        entries through the constructor without triggering a check against a fixed default schema.

  @ingroup FeatureGrouping
*/
class OPENMS_DLLAPI KDTreeFeatureMaps : public DefaultParamHandler
{

public:

  /// 2D kd-tree on (RT, m/z) holding @ref KDTreeFeatureNode entries
  typedef KDTree::KDTree<2,KDTreeFeatureNode> FeatureKDTree;

  /// Default constructor; creates an empty container with no features and zero maps
  KDTreeFeatureMaps() :
    DefaultParamHandler("KDTreeFeatureMaps")
  {
    check_defaults_ = false;
  }

  /**
    @brief Construct from a list of feature maps and a parameter set.

    Equivalent to default-constructing, calling setParameters(@p param), and addMaps(@p maps).

    @param[in] maps  Container of feature-map-like objects (each iterable over @ref BaseFeature
                     subclasses). All features from every map are inserted into the kd-tree.
    @param[in] param Parameter set forwarded to setParameters().
  */
  template <typename MapType>
  KDTreeFeatureMaps(const std::vector<MapType>& maps, const Param& param) :
    DefaultParamHandler("KDTreeFeatureMaps")
  {
    check_defaults_ = false;
    setParameters(param);
    addMaps(maps);
  }

  /// Destructor (does not delete the wrapped features)
  ~KDTreeFeatureMaps() override
  {
  }

  /**
    @brief Insert all features from @p maps and balance the kd-tree.

    Stores @p maps.size() as the number of input maps. Each feature is added via addFeature()
    using its map index, then optimizeTree() is called to balance the tree.

    @param[in] maps Container of feature-map-like objects (each iterable over @ref BaseFeature
                    subclasses). Pointers to the contained features are stored, so the maps must
                    outlive this object.
  */
  template <typename MapType>
  void addMaps(const std::vector<MapType>& maps)
  {
    num_maps_ = maps.size();

    for (Size i = 0; i < num_maps_; ++i)
    {
      const MapType& m = maps[i];
      for (typename MapType::const_iterator it = m.begin(); it != m.end(); ++it)
      {
        addFeature(i, &(*it));
      }
    }
    optimizeTree();
  }

  /**
    @brief Insert a single feature into the kd-tree.

    The feature's RT is cached internally (so it can be re-mapped by applyTransformations() later);
    m/z, intensity, and charge are always read through @p feature on demand.

    After many incremental inserts call optimizeTree() to re-balance the kd-tree.

    @param[in] mt_map_index Index of the source map this feature belongs to.
    @param[in] feature      Pointer to the feature to add (must outlive this object).
  */
  void addFeature(Size mt_map_index, const BaseFeature* feature);

  /// Pointer to the feature with global index @p i
  const BaseFeature* feature(Size i) const;

  /// (Possibly transformed) retention time of feature @p i
  double rt(Size i) const;

  /// m/z of feature @p i (read through the wrapped @ref BaseFeature)
  double mz(Size i) const;

  /// Intensity of feature @p i (read through the wrapped @ref BaseFeature)
  float intensity(Size i) const;

  /// Charge of feature @p i (read through the wrapped @ref BaseFeature)
  Int charge(Size i) const;

  /// Index of the source map that feature @p i was added with
  Size mapIndex(Size i) const;

  /// Number of features stored (across all maps)
  Size size() const;

  /// Number of points currently in the kd-tree
  Size treeSize() const;

  /// Number of source maps recorded in the most recent addMaps() call
  Size numMaps() const;

  /// Discard all stored features and clear the kd-tree (does not delete the wrapped features)
  void clear();

  /// Re-balance the kd-tree after incremental inserts
  void optimizeTree();

  /**
    @brief Collect indices of all features compatible with the feature at @p index.

    Performs a queryRegion() with an RT/m/z tolerance window centred on the feature at @p index.
    By default, features from the same source map as @p index are excluded (typical setting for
    cross-map matching). If @p max_pairwise_log_fc is non-negative, each candidate is additionally
    filtered by |log10(I_candidate / I_index)| <= @p max_pairwise_log_fc; candidates whose log
    fold change is NaN or +inf (e.g. due to zero/negative intensities) are silently dropped.

    Results are appended to @p result_indices (the existing contents are preserved).

    @param[in]  index                          Global index of the query feature.
    @param[in,out] result_indices              Indices of neighbouring features are appended here.
    @param[in]  rt_tol                         Tolerance on retention time (absolute, same unit as RT).
    @param[in]  mz_tol                         Tolerance on m/z (interpreted as ppm if @p mz_ppm is true, otherwise as Th).
    @param[in]  mz_ppm                         If true, @p mz_tol is treated as parts-per-million.
    @param[in]  include_features_from_same_map If true, features from the same source map as @p index are kept too.
    @param[in]  max_pairwise_log_fc            If >= 0, filter on absolute log10 intensity ratio; if < 0 (default), no intensity filter is applied.
  */
  void getNeighborhood(Size index, std::vector<Size>& result_indices, double rt_tol, double mz_tol, bool mz_ppm, bool include_features_from_same_map = false, double max_pairwise_log_fc = -1.0) const;

  /**
    @brief Range query on the kd-tree: collect indices of features within the given (RT, m/z) rectangle.

    Replaces (does not append to) @p result_indices.

    @param[in]  rt_low             Lower bound on RT.
    @param[in]  rt_high            Upper bound on RT.
    @param[in]  mz_low             Lower bound on m/z.
    @param[in]  mz_high            Upper bound on m/z.
    @param[out] result_indices     Indices of features inside the rectangle (cleared first, then filled).
    @param[in]  ignored_map_index  If different from std::numeric_limits<Size>::max() (the default),
                                   features whose map index equals @p ignored_map_index are excluded
                                   from the result.
  */
  void queryRegion(double rt_low, double rt_high, double mz_low, double mz_high, std::vector<Size>& result_indices, Size ignored_map_index = std::numeric_limits<Size>::max()) const;

  /**
    @brief Apply a per-map RT transformation to the cached retention times.

    For each feature @c i, replaces the cached RT by @c trafos[mapIndex(i)]->evaluate(feature(i)->getRT()).
    The underlying @ref BaseFeature objects are not modified.

    @param[in] trafos One transformation per source map. Size must equal numMaps().
  */
  void applyTransformations(const std::vector<TransformationModelLowess*>& trafos);

protected:

  void updateMembers_() override;

  /// Feature data (pointers owned by the caller)
  std::vector<const BaseFeature*> features_;

  /// Map of origin for each feature (parallel to @c features_)
  std::vector<Size> map_index_;

  /// (Potentially transformed) retention times, parallel to @c features_
  std::vector<double> rt_;

  /// Number of source maps (set by addMaps())
  Size num_maps_;

  /// 2D kd-tree on features from all input maps (queried in (RT, m/z))
  FeatureKDTree kd_tree_;

};
}

