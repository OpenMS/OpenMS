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
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERFUNCTOR_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERFUNCTOR_H

#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>

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

        @p original_distance is considered mirrored at the main diagonal, so only entrys up the main diagonal are used.
        The @p threshold can be taken from the maximal distance of two elements considered related and adapted in a way corresponding to the employed clustering method.
        The results are represented by @p cluster_tree, to get the actual clustering (with element indices) from a certain step of the clustering
        @see BinaryTreeNode , ClusterAnalyzer::cut
    */
    virtual void operator()(DistanceMatrix<float> & original_distance, std::vector<BinaryTreeNode> & cluster_tree, const float threshold = 1) const = 0;

    /// registers all derived products
    static void registerChildren();

  };

}
#endif // OPENMS_COMPARISON_CLUSTERFUNCTOR_H
