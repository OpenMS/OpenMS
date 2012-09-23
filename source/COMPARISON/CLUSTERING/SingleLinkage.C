// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
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

  SingleLinkage::~SingleLinkage()
  {
  }

  SingleLinkage & SingleLinkage::operator=(const SingleLinkage & source)
  {
    if (this != &source)
    {
      ClusterFunctor::operator=(source);
      ProgressLogger::operator=(source);
    }
    return *this;
  }

  void SingleLinkage::operator()(DistanceMatrix<Real> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const Real threshold /*=1*/) const
  {
    // input MUST have >= 2 elements!
    if (original_distance.dimensionsize() < 2)
    {
      throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Distance matrix to start from only contains one element");
    }

    cluster_tree.clear();
    if (threshold < 1)
    {
      LOG_ERROR << "You tried to use Single Linkage clustering with a threshold. This is currently not supported!" << std::endl;
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }

    //SLINK
    std::vector<Size> pi;
    pi.reserve(original_distance.dimensionsize());
    std::vector<Real> lambda;
    lambda.reserve(original_distance.dimensionsize());

    startProgress(0, original_distance.dimensionsize(), "clustering data");

    //initialize first pointer values
    pi.push_back(0);
    lambda.push_back(std::numeric_limits<Real>::max());

    for (Size k = 1; k < original_distance.dimensionsize(); ++k)
    {
      std::vector<Real> row_k;
      row_k.reserve(k);

      //initialize pointer values for element to cluster
      pi.push_back(k);
      lambda.push_back(std::numeric_limits<Real>::max());

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

      //update clustering if neccessary
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
      cluster_tree.push_back(BinaryTreeNode(i, pi[i], lambda[i]));
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
