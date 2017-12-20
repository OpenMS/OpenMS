// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#ifndef OPENMS_COMPARISON_CLUSTERING_AVERAGELINKAGE_H
#define OPENMS_COMPARISON_CLUSTERING_AVERAGELINKAGE_H

#include <vector>
#include <set>

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>

namespace OpenMS
{
  /**
      @brief AverageLinkage ClusterMethod

      The details of the method can be found in:
      Backhaus, Erichson, Plinke, Weiber Multivariate Analysemethoden, Springer 2000 and
      Ellen M. Voorhees: Implementing agglomerative hierarchic clustering algorithms for use in document retrieval. Inf. Process. Manage. 22(6): 465-476 (1986)
      @see ClusterFunctor

      @ingroup SpectraClustering
  */
  class OPENMS_DLLAPI AverageLinkage :
    public ClusterFunctor, public ProgressLogger
  {
public:

    /// default constructor
    AverageLinkage();

    /// copy constructor
    AverageLinkage(const AverageLinkage & source);

    /// destructor
    ~AverageLinkage() override;

    /// assignment operator
    AverageLinkage & operator=(const AverageLinkage & source);


    /**
        @brief clusters the indices according to their respective element distances

        @param original_distance DistanceMatrix<float> containing the distances of the elements to be clustered, will be changed during clustering process, make sure to have a copy or be able to redo
        @param cluster_tree vector< BinaryTreeNode >, represents the clustering, each node contains the next merged clusters (not element indices) and their distance, strict order is kept: left_child < right_child
        @param threshold float value, the minimal distance from which on cluster merging is considered unrealistic. By default set to 1, i.e. complete clustering until only one cluster remains
        @throw ClusterFunctor::InsufficientInput thrown if input is <2
        The clustering method is average linkage, where the updated distances after merging two clusters are each the average distances between the elements of their clusters. After @p threshold is exceeded, @p cluster_tree is filled with dummy clusteringsteps (children: (0,1), distance: -1) to the root.
        @see ClusterFunctor , BinaryTreeNode
    */
    void operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold = 1) const override;

    /// creates a new instance of a AverageLinkage object
    static ClusterFunctor * create();

    /// get the identifier for this object
    static const String getProductName();

  };

}
#endif //OPENMS_COMPARISON_CLUSTERING_AVERAGELINKAGE_H
