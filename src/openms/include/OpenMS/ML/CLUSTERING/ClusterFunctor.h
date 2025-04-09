// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/ML/CLUSTERING/ClusterAnalyzer.h>

#include <vector>

namespace OpenMS
{

  /**
      @brief Base class for cluster functors

      Each cluster functor employs a different method for stepwise merging clusters up to a given threshold, starting from the most elementary partition of data. Elements are represented by indices of a given distance matrix, which also should represent the order of input.

      @ingroup SpectraClustering
  */
  class OPENMS_DLLAPI ClusterFunctor
  {

public:

    /**
        @brief Exception thrown if not enough data (<2) is used

            If the set of data to be clustered contains only one data point,
            clustering algorithms would fail for obvious reasons.
    */
    class OPENMS_DLLAPI InsufficientInput :
      public Exception::BaseException
    {
public:
      InsufficientInput(const char * file, int line, const char * function, const char * message = "not enough data points to cluster anything") throw();
      ~InsufficientInput() throw() override;
    };


    /// default constructor
    ClusterFunctor();

    /// copy constructor
    ClusterFunctor(const ClusterFunctor & source);

    /// destructor
    virtual ~ClusterFunctor();

    /// assignment operator
    ClusterFunctor & operator=(const ClusterFunctor & source);

    /**
        @brief abstract for clustering the indices according to their respective element distances

        @param original_distance DistanceMatrix<float> containing the distances of the elements to be clustered, will be changed during clustering process, make sure to have a copy or be able to redo
        @param cluster_tree vector< BinaryTreeNode >, represents the clustering, each node contains the next merged clusters (not element indices) and their distance, strict order is kept: left_child < right_child,
        @param threshold float value, the minimal distance from which on cluster merging is considered unrealistic. By default set to 1, i.e. complete clustering until only one cluster remains

        @p original_distance is considered mirrored at the main diagonal, so only entries up the main diagonal are used.
        The @p threshold can be taken from the maximal distance of two elements considered related and adapted in a way corresponding to the employed clustering method.
        The results are represented by @p cluster_tree, to get the actual clustering (with element indices) from a certain step of the clustering
        @see BinaryTreeNode , ClusterAnalyzer::cut
    */
    virtual void operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold = 1) const = 0;



  };

}
