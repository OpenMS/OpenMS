// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <vector>
#include <set>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/ML/CLUSTERING/ClusterFunctor.h>

namespace OpenMS
{
  /**
      @brief SingleLinkage ClusterMethod

      The details of the method can be found in:
      SLINK: An optimally efficient algorithm for the single-link cluster method, The Computer Journal 1973 16(1):30-34; doi:10.1093/comjnl/16.1.30
      @see ClusterFunctor() base class.

      @ingroup SpectraClustering
  */
  class OPENMS_DLLAPI SingleLinkage :
    public ClusterFunctor, public ProgressLogger
  {
public:

    /// default constructor
    SingleLinkage();

    /// copy constructor
    SingleLinkage(const SingleLinkage & source);

    /// destructor
    ~SingleLinkage() override;

    /// assignment operator
    SingleLinkage & operator=(const SingleLinkage & source);

    /**
        @brief clusters the indices according to their respective element distances

    @param original_distance DistanceMatrix<float> containing the distances of the elements to be clustered
    @param cluster_tree vector< BinaryTreeNode >, represents the clustering, each node contains the next two clusters merged and their distance, strict order is kept: left_child < right_child
    @param threshold float value to meet Base class interface, will not be used because algorithm used is considerably fast and does not work by growing distances
    @throw ClusterFunctor::InsufficientInput thrown if input is <2
        The clustering method is single linkage, where the updated distances after merging two clusters are each the minimal distance between the elements of their clusters.
    @see ClusterFunctor , BinaryTreeNode
    */
    void operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold = 1) const override;

  };



}
