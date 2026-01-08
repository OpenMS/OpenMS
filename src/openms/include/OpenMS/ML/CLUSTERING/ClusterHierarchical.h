// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
//

//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/ML/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/ML/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/KERNEL/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>

namespace OpenMS
{

/**
    @brief Hierarchical clustering with generic clustering functions

    ClusterHierarchical clusters objects with corresponding distancemethod and clusteringmethod.
    @ingroup SpectraClustering
*/
class OPENMS_DLLAPI ClusterHierarchical
{
private:
  /// the threshold given to the ClusterFunctor
  double threshold_;

public:
  /// default constructor
  ClusterHierarchical(): threshold_(1.0)
  {
  }

  /// copy constructor
  ClusterHierarchical(const ClusterHierarchical& source): threshold_(source.threshold_)
  {
  }

  /// destructor
  virtual ~ClusterHierarchical()
  {
  }

  /**
      @brief Clustering function

      Cluster data using SimilarityComparator and ClusterFunctor.

      Creates a DistanceMatrix (if an empty matrix is passed) and the clustering is started.
      Clustering stops if the ClusterHierarchical::threshold_ is reached by the ClusterFunctor.

      First template parameter is the cluster object type,
      Second template parameter is the similarity functor applicable to the type.

      For example, @ref PeakSpectrum with a @ref PeakSpectrumCompareFunctor.

      The similarity functor must provide the similarity calculation with the ()-operator and
      yield normalized values in range of [0,1] for the type of < Data >.

      @param data Values to be clustered
      @param comparator Similarity functor which returns a similarity in [0, 1] for any pair of values in @p data
      @param clusterer A cluster method implementation, e.g. SingleLinkage or CompleteLinkage. See base class ClusterFunctor.
      @param cluster_tree The vector that will hold the BinaryTreeNodes representing the clustering (for further investigation with the
     ClusterAnalyzer methods)
      @param original_distance Precomputed DistanceMatrix holding the pairwise distances of the elements in @p data; if empty or wrong size, it will be re-created using @p comparator

      @see ClusterFunctor, BinaryTreeNode, ClusterAnalyzer
  */
  template<typename Data, typename SimilarityComparator>
  void cluster(std::vector<Data>& data,
               const SimilarityComparator& comparator,
               const ClusterFunctor& clusterer,
               std::vector<BinaryTreeNode>& cluster_tree,
               DistanceMatrix<float>& original_distance)
  {
    if (original_distance.dimensionsize() != data.size())
    {
      // create distance matrix for data using comparator
      original_distance.clear();
      original_distance.resize(data.size(), 1);
      for (Size i = 0; i < data.size(); i++)
      {
        for (Size j = 0; j < i; j++)
        {
          // distance value is 1-similarity value, since similarity is in range of [0,1]
          original_distance.setValueQuick(i, j, 1 - comparator(data[i], data[j]));
        }
      }
    }

    // create clustering with ClusterMethod, DistanceMatrix and Data
    clusterer(original_distance, cluster_tree, threshold_);
  }

  /**
    @brief clustering function for binned PeakSpectrum

    A version of the clustering function for PeakSpectra employing binned similarity methods. From the given PeakSpectrum BinnedSpectrum are
    generated, so the similarity functor @see BinnedSpectrumCompareFunctor can be applied.

    @param data vector of @ref PeakSpectrum s to be clustered
    @param comparator a BinnedSpectrumCompareFunctor
    @param sz Bin size for the @ref BinnedSpectrum s
    @param sp Bin spread for the @ref BinnedSpectrum s
    @param offset Bin offset for the @ref BinnedSpectrum s
    @param clusterer a clustermethod implementation, base class ClusterFunctor
    @param cluster_tree Vector of BinaryTreeNodes representing the clustering (for further investigation with the ClusterAnalyzer methods)
    @param original_distance The DistanceMatrix holding the pairwise distances of the elements in @p data, will be made newly if given size does not
                             fit to the number of elements given in @p data

    @see ClusterFunctor, BinaryTreeNode, ClusterAnalyzer, BinnedSpectrum, BinnedSpectrumCompareFunctor

    @ingroup SpectraClustering
  */
  void cluster(std::vector<PeakSpectrum>& data,
               const BinnedSpectrumCompareFunctor& comparator,
               double sz,
               UInt sp,
               float offset,
               const ClusterFunctor& clusterer,
               std::vector<BinaryTreeNode>& cluster_tree,
               DistanceMatrix<float>& original_distance) const
  {
    std::vector<BinnedSpectrum> binned_data;
    binned_data.reserve(data.size());

    // transform each PeakSpectrum to a corresponding BinnedSpectrum with given settings of size and spread
    for (Size i = 0; i < data.size(); i++)
    {
      // double sz(2), UInt sp(1);
      binned_data.emplace_back(data[i], sz, false, sp, offset);
    }

    // create distancematrix for data with comparator
    original_distance.clear();
    original_distance.resize(data.size(), 1);

    for (Size i = 0; i < binned_data.size(); i++)
    {
      for (Size j = 0; j < i; j++)
      {
        // distance value is 1-similarity value, since similarity is in range of [0,1]
        original_distance.setValue(i, j, 1 - comparator(binned_data[i], binned_data[j]));
      }
    }

    // create Clustering with ClusterMethod, DistanceMatrix and Data
    clusterer(original_distance, cluster_tree, threshold_);
  }

  /// get the threshold
  double getThreshold() const
  {
    return threshold_;
  }

  /// set the threshold (in terms of distance)
  /// The default is 1, i.e. only at similarity 0 the clustering stops.
  /// Warning: clustering is not supported by all methods yet (e.g. SingleLinkage does ignore it).
  void setThreshold(double x)
  {
    threshold_ = x;
  }
};

/** @brief Exception thrown if clustering is attempted without a normalized compare functor

        due to similarity - distance conversions that are mandatory in some context, compare functors
        must return values normalized in the range [0,1] to ensure a clean conversion
*/
class OPENMS_DLLAPI UnnormalizedComparator : public Exception::BaseException
{
public:
  UnnormalizedComparator(const char* file,
                         int line,
                         const char* function,
                         const char* message = "Clustering with unnormalized similarity measurement requested, normalized is mandatory") throw();
  ~UnnormalizedComparator() throw() override;
};

} // namespace OpenMS
