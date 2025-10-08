// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{

///
/**
    @brief Proxy for a (potential) cluster.

    Proxy for a (potential) cluster. Instead of storing the entire cluster,
    this stores only its size, average distance to center, and the index of
    the center point. Objects of this class are kept in a sorted binary search
    tree (aka std::set) and operator< is defined in such a way that the first
    element of the set is always a cluster proxy for a cluster of current
    maximum size and smallest intra-cluster distance. The actual cluster points
    are then retrieved again from the kd-tree and a consensus feature is added
    to the output consensus map.
*/
class OPENMS_DLLAPI ClusterProxyKD
{

public:

  /// Default constructor
  ClusterProxyKD() :
    size_(0), // => isValid() returns false
    avg_distance_(0),
    center_index_(0)
  {
  }

  /// Constructor
  ClusterProxyKD(Size size, double avg_distance, Size center_index) :
    size_(size),
    avg_distance_(avg_distance),
    center_index_(center_index)
  {
  }

  /// Copy constructor
  ClusterProxyKD(const ClusterProxyKD& rhs) :
    size_(rhs.size_),
    avg_distance_(rhs.avg_distance_),
    center_index_(rhs.center_index_)
  {
  }

  /// Destructor (non-virtual to save memory)
  ~ClusterProxyKD()
  {
  }

  /// Assignment operator
  ClusterProxyKD& operator=(const ClusterProxyKD& rhs)
  {
    size_ = rhs.size_;
    avg_distance_ = rhs.avg_distance_;
    center_index_ = rhs.center_index_;

    return *this;
  }

  /// Less-than operator for sorting / equality check in std::set. We use the ordering in std::set as a "priority queue", hence a < b means cluster a will be preferred over b.
  bool operator<(const ClusterProxyKD& rhs) const
  {
    if (size_ > rhs.size_) return true;
    if (size_ < rhs.size_) return false;

    if (avg_distance_ < rhs.avg_distance_) return true;
    if (avg_distance_ > rhs.avg_distance_) return false;

    // arbitrary, but required for finding unambiguous elements in std::set
    if (center_index_ > rhs.center_index_) return true;
    if (center_index_ < rhs.center_index_) return false;

    // they are equal
    return false;
  }

  /// Inequality operator
  bool operator!=(const ClusterProxyKD& rhs) const
  {
    return *this < rhs || rhs < *this;
  }

  /// Equality operator
  bool operator==(const ClusterProxyKD& rhs) const
  {
    return !(*this != rhs);
  }

  /// Cluster size
  Size getSize() const
  {
    return size_;
  }

  /// Valid?
  bool isValid() const
  {
    return size_;
  }

  /// Average distance to center
  double getAvgDistance() const
  {
    return avg_distance_;
  }

  /// Index of center point
  Size getCenterIndex() const
  {
    return center_index_;
  }

private:

  /// Cluster size
  Size size_;

  /// Average distance to center
  double avg_distance_;

  /// Index of center point
  Size center_index_;
};


  /**
      @brief A feature grouping algorithm for unlabeled data.

      The algorithm takes a number of feature or consensus maps and searches
      for corresponding (consensus) features across different maps.

      @htmlinclude OpenMS_FeatureGroupingAlgorithmKD.parameters

      @ingroup FeatureGrouping
  */
  class OPENMS_DLLAPI FeatureGroupingAlgorithmKD :
    public FeatureGroupingAlgorithm,
    public ProgressLogger
  {

public:

    /// Default constructor
    FeatureGroupingAlgorithmKD();

    /// Destructor
    ~FeatureGroupingAlgorithmKD() override;

    /**
        @brief Applies the algorithm to feature maps

        @exception IllegalArgument is thrown if less than two input maps are given.
    */
    void group(const std::vector<FeatureMap>& maps, ConsensusMap& out) override;

    /**
        @brief Applies the algorithm to consensus maps

        @exception IllegalArgument is thrown if less than two input maps are given.
    */
    void group(const std::vector<ConsensusMap>& maps,
                       ConsensusMap& out) override;

private:

    /// Copy constructor intentionally not implemented -> private
    FeatureGroupingAlgorithmKD(const FeatureGroupingAlgorithmKD&);

    /// Assignment operator intentionally not implemented -> private
    FeatureGroupingAlgorithmKD& operator=(const FeatureGroupingAlgorithmKD&);

    /**
        @brief Applies the algorithm to feature or consensus maps

        @exception IllegalArgument is thrown if less than two input maps are given.
    */
    template <typename MapType>
    void group_(const std::vector<MapType>& input_maps, ConsensusMap& out);

    /// Run the actual clustering algorithm
    void runClustering_(const KDTreeFeatureMaps& kd_data, ConsensusMap& out);

    /// Update maximum possible sizes of potential consensus features for indices specified in @p update_these
    void updateClusterProxies_(std::set<ClusterProxyKD>& potential_clusters, std::vector<ClusterProxyKD>& cluster_for_idx, const std::set<Size>& update_these, const std::vector<Int>& assigned, const KDTreeFeatureMaps& kd_data);

    /// Compute the current best cluster with center index @p i (mutates @p proxy and @p cf_indices)
    ClusterProxyKD computeBestClusterForCenter_(Size i, std::vector<Size>& cf_indices, const std::vector<Int>& assigned, const KDTreeFeatureMaps& kd_data) const;

    /// Construct consensus feature and add to out map
    void addConsensusFeature_(const std::vector<Size>& indices, const KDTreeFeatureMaps& kd_data, ConsensusMap& out) const;

    /// Current progress for logging
    SignedSize progress_;

    /// RT tolerance
    double rt_tol_secs_;

    /// m/z tolerance
    double mz_tol_;

    /// m/z unit ppm?
    bool mz_ppm_;

    /// Feature distance functor
    FeatureDistance feature_distance_;
  };

} // namespace OpenMS

