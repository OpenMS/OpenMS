// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ML/CLUSTERING/SingleLinkage.h>
#include <OpenMS/CONCEPT/LogStream.h>

namespace OpenMS
{

  SingleLinkage::SingleLinkage() :
    ClusterFunctor(), ProgressLogger()
  {
  }

  SingleLinkage::SingleLinkage(const SingleLinkage & source) :
    ClusterFunctor(source), ProgressLogger()
  {
  }

  SingleLinkage::~SingleLinkage() = default;

  SingleLinkage & SingleLinkage::operator=(const SingleLinkage & source)
  {
    if (this != &source)
    {
      ClusterFunctor::operator=(source);
      ProgressLogger::operator=(source);
    }
    return *this;
  }

  void SingleLinkage::operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold /*=1*/) const
  {
    // input MUST have >= 2 elements!
    if (original_distance.dimensionsize() < 2)
    {
      throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Distance matrix to start from only contains one element");
    }

    cluster_tree.clear();
    if (threshold < 1)
    {
      OPENMS_LOG_ERROR << "You tried to use Single Linkage clustering with a threshold. This is currently not supported!" << std::endl;
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    //SLINK
    std::vector<Size> pi;
    pi.reserve(original_distance.dimensionsize());
    std::vector<float> lambda;
    lambda.reserve(original_distance.dimensionsize());

    startProgress(0, original_distance.dimensionsize(), "clustering data");

    //initialize first pointer values
    pi.push_back(0);
    lambda.push_back(std::numeric_limits<float>::max());

    for (Size k = 1; k < original_distance.dimensionsize(); ++k)
    {
      std::vector<float> row_k;
      row_k.reserve(k);

      //initialize pointer values for element to cluster
      pi.push_back(k);
      lambda.push_back(std::numeric_limits<float>::max());

      // get the right distances
      for (Size i = 0; i < k; ++i)
      {
        row_k.push_back(original_distance.getValue(i, k));
      }

      //calculate pointer values for element k
      for (Size i = 0; i < k; ++i)
      {
        if (lambda[i] >= row_k[i])
        {
          row_k[pi[i]] = std::min(row_k[pi[i]], lambda[i]);
          lambda[i] = row_k[i];
          pi[i] = k;
        }
        else
        {
          row_k[pi[i]] = std::min(row_k[pi[i]], row_k[i]);
        }
      }

      //update clustering if necessary
      for (Size i = 0; i < k; ++i)
      {
        if (lambda[i] >= lambda[pi[i]])
        {
          pi[i] = k;
        }
      }
      setProgress(k);
    }

    for (Size i = 0; i < pi.size() - 1; ++i)
    {
      //strict order is always kept in algorithm: i < pi[i]
      cluster_tree.emplace_back(i, pi[i], lambda[i]);
      //~ std::cout << i << '\n' << pi[i] << '\n' << lambda[i] << std::endl;
    }

    //sort pre-tree
    std::sort(cluster_tree.begin(), cluster_tree.end(), compareBinaryTreeNode);

    // convert -pre-tree to correct format
    for (Size i = 0; i < cluster_tree.size(); ++i)
    {
      if (cluster_tree[i].right_child < cluster_tree[i].left_child)
      {
        std::swap(cluster_tree[i].left_child, cluster_tree[i].right_child);
      }
      for (Size k = i + 1; k < cluster_tree.size(); ++k)
      {
        if (cluster_tree[k].left_child == cluster_tree[i].right_child)
        {
          cluster_tree[k].left_child = cluster_tree[i].left_child;
        }
        else if (cluster_tree[k].left_child > cluster_tree[i].right_child)
        {
          --cluster_tree[k].left_child;
        }
        if (cluster_tree[k].right_child == cluster_tree[i].right_child)
        {
          cluster_tree[k].right_child = cluster_tree[i].left_child;
        }
        else if (cluster_tree[k].right_child > cluster_tree[i].right_child)
        {
          --cluster_tree[k].right_child;
        }
      }

    }
    //~ prepare to redo clustering to get all indices for binarytree in min index element representation
    std::vector<std::set<Size> > clusters(original_distance.dimensionsize());
    for (Size i = 0; i < original_distance.dimensionsize(); ++i)
    {
      clusters[i].insert(i);
    }
    for (Size cluster_step = 0; cluster_step < cluster_tree.size(); ++cluster_step)
    {
      Size new_left_child = *(clusters[cluster_tree[cluster_step].left_child].begin());
      Size new_right_child = *(clusters[cluster_tree[cluster_step].right_child].begin());
      clusters[cluster_tree[cluster_step].left_child].insert(clusters[cluster_tree[cluster_step].right_child].begin(), clusters[cluster_tree[cluster_step].right_child].end());
      clusters.erase(clusters.begin() + cluster_tree[cluster_step].right_child);
      std::swap(cluster_tree[cluster_step].left_child, new_left_child);
      std::swap(cluster_tree[cluster_step].right_child, new_right_child);
      if (cluster_tree[cluster_step].left_child > cluster_tree[cluster_step].right_child)
      {
        std::swap(cluster_tree[cluster_step].left_child, cluster_tree[cluster_step].right_child);
      }
    }

    endProgress();
  }

}
