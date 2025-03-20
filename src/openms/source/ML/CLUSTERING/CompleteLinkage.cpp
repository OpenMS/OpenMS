// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ML/CLUSTERING/CompleteLinkage.h>

#include <OpenMS/DATASTRUCTURES/String.h>

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

    std::vector<std::set<Size> > clusters(original_distance.dimensionsize());
    for (Size i = 0; i < original_distance.dimensionsize(); ++i)
    {
      clusters[i].insert(i);
    }

    cluster_tree.clear();
    cluster_tree.reserve(original_distance.dimensionsize() - 1);

    // Initial minimum-distance pair
    original_distance.updateMinElement();
    std::pair<Size, Size> min = original_distance.getMinElementCoordinates();

    Size overall_cluster_steps(original_distance.dimensionsize());
    startProgress(0, original_distance.dimensionsize(), "clustering data");

    while (original_distance(min.first, min.second) < threshold)
    {
      //grow the tree
      cluster_tree.emplace_back(*(clusters[min.second].begin()), *(clusters[min.first].begin()), original_distance(min.first, min.second));
      if (cluster_tree.back().left_child > cluster_tree.back().right_child)
      {
        std::swap(cluster_tree.back().left_child, cluster_tree.back().right_child);
      }

      if (original_distance.dimensionsize() > 2)
      {
        //pick minimum-distance pair i,j and merge them

        //pushback elements of second to first (and then erase second)
        clusters[min.second].insert(clusters[min.first].begin(), clusters[min.first].end());
        // erase first one
        clusters.erase(clusters.begin() + min.first);

        //update original_distance matrix
        //complete linkage: new distance between clusters is the minimum distance between elements of each cluster
        //lance-williams update for d((i,j),k): 0.5* d(i,k) + 0.5* d(j,k) + 0.5* |d(i,k)-d(j,k)|
        for (Size k = 0; k < min.second; ++k)
        {
          float dik = original_distance.getValue(min.first, k);
          float djk = original_distance.getValue(min.second, k);
          original_distance.setValueQuick(min.second, k, (0.5f * dik + 0.5f * djk + 0.5f * std::fabs(dik - djk)));
        }
        for (Size k = min.second + 1; k < original_distance.dimensionsize(); ++k)
        {
          float dik = original_distance.getValue(min.first, k);
          float djk = original_distance.getValue(min.second, k);
          original_distance.setValueQuick(k, min.second, (0.5f * dik + 0.5f * djk + 0.5f * std::fabs(dik - djk)));
        }

        //reduce
        original_distance.reduce(min.first);

        //update minimum-distance pair
        original_distance.updateMinElement();

        //get new min-pair
        min = original_distance.getMinElementCoordinates();
      }
      else
      {
        break;
      }
      setProgress(overall_cluster_steps - original_distance.dimensionsize());

      //repeat until only two cluster remains or threshold exceeded, last step skips matrix operations
    }
    //fill tree with dummy nodes
    Size sad(*clusters.front().begin());
    for (Size i = 1; i < clusters.size() && (cluster_tree.size() < cluster_tree.capacity()); ++i)
    {
      cluster_tree.emplace_back(sad, *clusters[i].begin(), -1.0);
    }
    //~ while(cluster_tree.size() < cluster_tree.capacity())
    //~ {
    //~ cluster_tree.push_back(BinaryTreeNode(0,1,-1.0));
    //~ }

    endProgress();
  }

}
