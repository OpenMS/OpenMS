// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ML/CLUSTERING/CompleteLinkage.h>

#include "LinkageHelpers.h"

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

namespace OpenMS
{

  CompleteLinkage::CompleteLinkage() :
    ClusterFunctor(), ProgressLogger()
  {
  }

  CompleteLinkage::CompleteLinkage(const CompleteLinkage & source) :
    ClusterFunctor(source), ProgressLogger()
  {
  }

  CompleteLinkage::~CompleteLinkage() = default;

  CompleteLinkage & CompleteLinkage::operator=(const CompleteLinkage & source)
  {
    if (this != &source)
    {
      ClusterFunctor::operator=(source);
      ProgressLogger::operator=(source);
    }
    return *this;
  }

  void CompleteLinkage::operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold /*=1*/) const
  {
    // attention: clustering process is done by clustering the indices
    // pointing to elements in input vector and distances in input matrix

    // input MUST have >= 2 elements!
    if (original_distance.dimensionsize() < 2)
    {
      throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Distance matrix to start from only contains one element");
    }

    //complete linkage: new distance between clusters is the maximum distance between elements of each cluster
    //lance-williams update for d((i,j),k): 0.5* d(i,k) + 0.5* d(j,k) + 0.5* |d(i,k)-d(j,k)|
    Internal::LinkageHelpers::lanceWilliamsClustering(original_distance, cluster_tree, threshold, *this,
      [](const Size /*m_i*/, const Size /*m_j*/)
      {
        return [](const float dik, const float djk) -> float { return (0.5f * dik + 0.5f * djk + 0.5f * std::fabs(dik - djk)); };
      });
  }

}
