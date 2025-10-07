// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>

#include <vector>

namespace OpenMS
{
  class String;

  /**
    @brief Bundles analyzing tools for a clustering (given as sequence of BinaryTreeNode's)

    @ingroup SpectraClustering
  */
  class OPENMS_DLLAPI ClusterAnalyzer
  {
public:
    /// default constructor
    ClusterAnalyzer();

    /// copy constructor
    ClusterAnalyzer(const ClusterAnalyzer & source);

    /// destructor
    virtual ~ClusterAnalyzer();

    /**
        @brief Method to calculate the average silhouette widths for a clustering

        @param tree vector of BinaryTreeNode's representing the clustering
        @param original DistanceMatrix for all clustered elements started from
        @return a vector filled with the average silhouette widths for each cluster step

        The average silhouette width will be calculated for each clustering step beginning with the first step(n-1 cluster) ending with the last (1 cluster, average silhouette width is 0 by definition).
        @see BinaryTreeNode
    */
    std::vector<float> averageSilhouetteWidth(const std::vector<BinaryTreeNode> & tree, const DistanceMatrix<float> & original);

    /**
        @brief Method to calculate Dunns indices for a clustering

        @param tree vector of BinaryTreeNode's representing the clustering
        @param original DistanceMatrix for all clustered elements started from
        @param tree_from_singlelinkage true if tree was created by SingleLinkage, i.e. the distances are the minimal distances in increasing order and can be used to speed up the calculation
        @see BinaryTreeNode
    */
    std::vector<float> dunnIndices(const std::vector<BinaryTreeNode> & tree, const DistanceMatrix<float> & original, const bool tree_from_singlelinkage = false);

    /**
        @brief Method to calculate the cohesions of a certain partition

        @param clusters vector of vectors holding the clusters (with indices to the actual elements)
        @param original DistanceMatrix for all clustered elements started from
        @return a vector that holds the cohesions of each cluster given with @p clusters (order corresponds to @p clusters)
    */
    std::vector<float> cohesion(const std::vector<std::vector<Size> > & clusters, const DistanceMatrix<float> & original);

    /**
        @brief Method to calculate the average aberration from average population in partition resulting from a certain step in clustering

        @param cluster_quantity desired partition Size analogue to ClusterAnalyzer::cut
        @param tree vector of BinaryTreeNode's representing the clustering
        @throw invalid_parameter if desired clustering is invalid
        @return the average aberration from the average cluster population (number of elements/cluster_quantity) at cluster_quantity
        @see BinaryTreeNode
    */
    float averagePopulationAberration(Size cluster_quantity, std::vector<BinaryTreeNode> & tree);

/**
        @brief Method to calculate a partition resulting from a certain step in clustering given by the number of clusters

    If you want to fetch all clusters which were created with a threshold, you simply count the number of tree-nodes which are
    not -1, and subtract that from the number of leaves, to get the number of clusters formed
    , i.e. cluster_quantity = data.size() - real_leaf_count;

        @param cluster_quantity Size giving the number of clusters (i.e. starting elements - cluster_quantity = cluster step)
        @param tree vector of BinaryTreeNode's representing the clustering
        @param clusters vector of vectors holding the clusters (with indices to the actual elements)
        @throw invalid_parameter if desired clusterstep is invalid
        @see BinaryTreeNode

        after call of this method the argument clusters is filled corresponding to the given @p cluster_quantity with the indices of the elements clustered
    */
    void cut(const Size cluster_quantity, const std::vector<BinaryTreeNode> & tree, std::vector<std::vector<Size> > & clusters);

/**
        @brief Method to calculate subtrees from a given tree resulting from a certain step in clustering given by the number of clusters

        @param cluster_quantity Size giving the number of clusters (i.e. starting elements - cluster_quantity = cluster step)
        @param tree vector of BinaryTreeNode's representing the clustering
        @param subtrees vector of trees holding the trees, tree is composed of cut at given size
        @throw invalid_parameter if desired clusterstep is invalid
        @see BinaryTreeNode

        after call of this method the argument clusters is filled corresponding to the given @p cluster_quantity with the indices of the elements clustered
    */
    void cut(const Size cluster_quantity, const std::vector<BinaryTreeNode> & tree, std::vector<std::vector<BinaryTreeNode> > & subtrees);

/**
        @brief Returns the hierarchy described by a clustering tree as Newick-String

        @param tree vector of BinaryTreeNode's representing the clustering
        @param include_distance bool value indicating whether the distance shall be included to the string
        @see BinaryTreeNode

    */
    String newickTree(const std::vector<BinaryTreeNode> & tree, const bool include_distance = false);

private:
    /// assignment operator
    ClusterAnalyzer & operator=(const ClusterAnalyzer & source);

  };
  ///returns the value of (x.distance < y.distance) for use with sort
  bool compareBinaryTreeNode(const BinaryTreeNode & x, const BinaryTreeNode & y);

}
