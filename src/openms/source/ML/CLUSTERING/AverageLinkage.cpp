// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ML/CLUSTERING/AverageLinkage.h>

#include "LinkageHelpers.h"

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

namespace OpenMS
{

  AverageLinkage::AverageLinkage() :
    ClusterFunctor(), ProgressLogger()
  {
  }

  AverageLinkage::AverageLinkage(const AverageLinkage & source)  = default;

  AverageLinkage::~AverageLinkage() = default;

  AverageLinkage & AverageLinkage::operator=(const AverageLinkage & source)
  {
    if (this != &source)
    {
      ClusterFunctor::operator=(source);
      ProgressLogger::operator=(source);
    }
    return *this;
  }

  void AverageLinkage::operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold /*=1*/) const
  {
    // input MUST have >= 2 elements!
    if (original_distance.dimensionsize() < 2)
    {
      throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Distance matrix to start from only contains one element");
    }

    //average linkage: new distance between clusters is the average distance between elements of each cluster
    //lance-williams update for d((i,j),k): (m_i/m_i+m_j)* d(i,k) + (m_j/m_i+m_j)* d(j,k) ; m_x is the number of elements in cluster x
    Internal::LinkageHelpers::lanceWilliamsClustering(original_distance, cluster_tree, threshold, *this,
      [](const Size m_i, const Size m_j)
      {
        //calculate parameter for lance-williams formula
        const float alpha_i = (float)(m_i / (float)(m_i + m_j));
        const float alpha_j = (float)(m_j / (float)(m_i + m_j));
        return [alpha_i, alpha_j](const float dik, const float djk) -> float { return (alpha_i * dik + alpha_j * djk); };
      });
  }

}
