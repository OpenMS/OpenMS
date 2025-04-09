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

/// Stores a set of features, together with a 2D tree for fast search
class OPENMS_DLLAPI KDTreeFeatureMaps : public DefaultParamHandler
{

public:

  /// 2D tree on features
  typedef KDTree::KDTree<2,KDTreeFeatureNode> FeatureKDTree;

  /// Default constructor
  KDTreeFeatureMaps() :
    DefaultParamHandler("KDTreeFeatureMaps")
  {
    check_defaults_ = false;
  }

  /// Constructor
  template <typename MapType>
  KDTreeFeatureMaps(const std::vector<MapType>& maps, const Param& param) :
    DefaultParamHandler("KDTreeFeatureMaps")
  {
    check_defaults_ = false;
    setParameters(param);
    addMaps(maps);
  }

  /// Destructor
  ~KDTreeFeatureMaps() override
  {
  }

  /// Add @p maps and balance kd-tree
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

  /// Add feature
  void addFeature(Size mt_map_index, const BaseFeature* feature);

  /// Return pointer to feature i
  const BaseFeature* feature(Size i) const;

  /// RT
  double rt(Size i) const;

  /// m/z
  double mz(Size i) const;

  /// Intensity
  float intensity(Size i) const;

  /// Charge
  Int charge(Size i) const;

  /// Map index
  Size mapIndex(Size i) const;

  /// Number of features stored
  Size size() const;

  /// Number of points in the tree
  Size treeSize() const;

  /// Number of maps
  Size numMaps() const;

  /// Clear all data
  void clear();

  /// Optimize the kD tree
  void optimizeTree();

  /// Fill @p result with indices of all features compatible (wrt. RT, m/z, map index) to the feature with @p index
  void getNeighborhood(Size index, std::vector<Size>& result_indices, double rt_tol, double mz_tol, bool mz_ppm, bool include_features_from_same_map = false, double max_pairwise_log_fc = -1.0) const;

  /// Fill @p result with indices of all features within the specified boundaries
  void queryRegion(double rt_low, double rt_high, double mz_low, double mz_high, std::vector<Size>& result_indices, Size ignored_map_index = std::numeric_limits<Size>::max()) const;

  /// Apply RT transformations
  void applyTransformations(const std::vector<TransformationModelLowess*>& trafos);

protected:

  void updateMembers_() override;

  /// Feature data
  std::vector<const BaseFeature*> features_;

  /// Map indices
  std::vector<Size> map_index_;

  /// (Potentially transformed) retention times
  std::vector<double> rt_;

  /// Number of maps
  Size num_maps_;

  /// 2D tree on features from all input maps.
  FeatureKDTree kd_tree_;

};
}

