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

#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>

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

  CompleteLinkage::~CompleteLinkage()
  {
  }

  CompleteLinkage & CompleteLinkage::operator=(const CompleteLinkage & source)
  {
    if (this != &source)
    {
      ClusterFunctor::operator=(source);
      ProgressLogger::operator=(source);
    }
    return *this;
  }

  void CompleteLinkage::operator()(DistanceMatrix<Real> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const Real threshold /*=1*/) const
  {
    // attention: clustering process is done by clustering the indices
    // pointing to elements in inputvector and distances in inputmatrix

    // input MUST have >= 2 elements!
    if (original_distance.dimensionsize() < 2)
    {
      throw ClusterFunctor::InsufficientInput(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Distance matrix to start from only contains one element");
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
      cluster_tree.push_back(BinaryTreeNode(*(clusters[min.second].begin()), *(clusters[min.first].begin()), original_distance(min.first, min.second)));
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
        //complete linkage: new distcance between clusteres is the minimum distance between elements of each cluster
        //lance-williams update for d((i,j),k): 0.5* d(i,k) + 0.5* d(j,k) + 0.5* |d(i,k)-d(j,k)|
        for (Size k = 0; k < min.second; ++k)
        {
          Real dik = original_distance.getValue(min.first, k);
          Real djk = original_distance.getValue(min.second, k);
          original_distance.setValueQuick(min.second, k, (0.5f * dik + 0.5f * djk + 0.5f * std::fabs(dik - djk)));
        }
        for (Size k = min.second + 1; k < original_distance.dimensionsize(); ++k)
        {
          Real dik = original_distance.getValue(min.first, k);
          Real djk = original_distance.getValue(min.second, k);
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
      cluster_tree.push_back(BinaryTreeNode(sad, *clusters[i].begin(), -1.0));
    }
    //~ while(cluster_tree.size() < cluster_tree.capacity())
    //~ {
    //~ cluster_tree.push_back(BinaryTreeNode(0,1,-1.0));
    //~ }

    endProgress();
  }

}
