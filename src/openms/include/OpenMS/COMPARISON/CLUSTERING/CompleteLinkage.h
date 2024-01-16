// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <cmath>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>

namespace OpenMS
{
  /**
      @brief CompleteLinkage ClusterMethod

      The details of the method can be found in:
      Backhaus, Erichson, Plinke, Weiber Multivariate Analysemethoden, Springer 2000 and
      Ellen M. Voorhees: Implementing agglomerative hierarchic clustering algorithms for use in document retrieval. Inf. Process. Manage. 22(6): 465-476 (1986)
      @see ClusterFunctor

      @ingroup SpectraClustering
  */
  class OPENMS_DLLAPI CompleteLinkage :
    public ClusterFunctor, public ProgressLogger
  {
public:

    /// default constructor
    CompleteLinkage();

    /// copy constructor
    CompleteLinkage(const CompleteLinkage & source);

    /// destructor
    ~CompleteLinkage() override;

    /// assignment operator
    CompleteLinkage & operator=(const CompleteLinkage & source);

    /**
        @brief clusters the indices according to their respective element distances

    @param original_distance DistanceMatrix<float> containing the distances of the elements to be clustered, will be changed during clustering process, make sure to have a copy or be able to redo
    @param cluster_tree vector< BinaryTreeNode >, represents the clustering, each node contains the next merged clusters (not element indices) and their distance, strict order is kept: left_child < right_child
    @param threshold float value, the minimal distance from which on cluster merging is considered unrealistic. By default set to 1, i.e. complete clustering until only one cluster remains
    @throw ClusterFunctor::InsufficientInput thrown if input is <2
        The clustering method is complete linkage, where the updated distances after merging two clusters are each the maximal distance between the elements of their clusters. After @p threshold is exceeded, @p cluster_tree is filled with dummy clusteringsteps (children: (0,1), distance:-1) to the root.
    @see ClusterFunctor , BinaryTreeNode
    */
    void operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold = 1) const override;

    /// creates a new instance of a CompleteLinkage object
    static ClusterFunctor * create();

    /// get the identifier for this object
    static const String getProductName();

  };

}
