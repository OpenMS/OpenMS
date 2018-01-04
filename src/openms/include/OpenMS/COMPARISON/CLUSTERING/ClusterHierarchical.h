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
#ifndef OPENMS_COMPARISON_CLUSTERING_CLUSTERHIERARCHICAL_H
#define OPENMS_COMPARISON_CLUSTERING_CLUSTERHIERARCHICAL_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/CONCEPT/Exception.h>

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
    ClusterHierarchical() :
      threshold_(1.0)
    {
    }

    /// copy constructor
    ClusterHierarchical(const ClusterHierarchical & source) :
      threshold_(source.threshold_)
    {
    }

    /// destructor
    virtual ~ClusterHierarchical()
    {
    }

    /**
        @brief Clustering function

        Conducts the SimilarityComparator with a ClusterFunctor an produces a clustering.
        Will create a DistanceMatrix if not yet created and start the clustering up to the given ClusterHierarchical::threshold_ used for the ClusterFunctor.
        The type of the objects to be clustered has to be the first template argument, the
        similarity functor applicable to this type must be the second template argument, e.g.
        for @ref PeakSpectrum with a @ref PeakSpectrumCompareFunctor.
        The similarity functor must provide the similarity calculation with the ()-operator and
        yield normalized values in range of [0,1] for the type of < Data >.

        @param data vector of objects to be clustered
        @param comparator similarity functor fitting for types in data
        @param clusterer a clustermethod implementation, baseclass ClusterFunctor
        @param cluster_tree the vector that will hold the BinaryTreeNodes representing the clustering (for further investigation with the ClusterAnalyzer methods)
        @param original_distance the DistanceMatrix holding the pairwise distances of the elements in @p data, will be made newly if given size does not fit to the number of elements given in @ data
        @see ClusterFunctor, BinaryTreeNode, ClusterAnalyzer
    */
    template <typename Data, typename SimilarityComparator>
    void cluster(std::vector<Data> & data, const SimilarityComparator & comparator, const ClusterFunctor & clusterer, std::vector<BinaryTreeNode> & cluster_tree, DistanceMatrix<float> & original_distance)
    {
      if (original_distance.dimensionsize() != data.size())
      {
        //create distancematrix for data with comparator
        original_distance.clear();
        original_distance.resize(data.size(), 1);
        for (Size i = 0; i < data.size(); i++)
        {
          for (Size j = 0; j < i; j++)
          {
            //distance value is 1-similarity value, since similarity is in range of [0,1]
            original_distance.setValueQuick(i, j, 1 - comparator(data[i], data[j]));
          }
        }
      }

      //~ std::cout << "done" << std::endl; //maybe progress handler?
      // create clustering with ClusterMethod, DistanceMatrix and Data
      clusterer(original_distance, cluster_tree, threshold_);
    }

    /**
        @brief clustering function for binned PeakSpectrum

        A version of the clustering function for PeakSpectra employing binned similarity methods. From the given PeakSpectrum BinnedSpectrum are generated, so the similarity functor @see BinnedSpectrumCompareFunctor can be applied.

        @param data vector of @ref PeakSpectrum s to be clustered
        @param comparator a BinnedSpectrumCompareFunctor
        @param sz the desired binsize for the @ref BinnedSpectrum s
        @param sp the desired binspread for the @ref BinnedSpectrum s
        @param clusterer a clustermethod implementation, base class ClusterFunctor
        @param cluster_tree the vector that will hold the BinaryTreeNodes representing the clustering (for further investigation with the ClusterAnalyzer methods)
        @param original_distance the DistanceMatrix holding the pairwise distances of the elements in @p data, will be made newly if given size does not fit to the number of elements given in @p data
        @see ClusterFunctor, BinaryTreeNode, ClusterAnalyzer, BinnedSpectrum, BinnedSpectrumCompareFunctor

    @ingroup SpectraClustering
    */
    void cluster(std::vector<PeakSpectrum> & data, const BinnedSpectrumCompareFunctor & comparator, double sz, UInt sp, const ClusterFunctor & clusterer, std::vector<BinaryTreeNode> & cluster_tree, DistanceMatrix<float> & original_distance)
    {

      std::vector<BinnedSpectrum> binned_data;
      binned_data.reserve(data.size());

      //transform each PeakSpectrum to a corresponding BinnedSpectrum with given settings of size and spread
      for (Size i = 0; i < data.size(); i++)
      {
        //double sz(2), UInt sp(1);
        binned_data.push_back(BinnedSpectrum(sz, sp, data[i]));
      }

      //create distancematrix for data with comparator
      original_distance.clear();
      original_distance.resize(data.size(), 1);

      for (Size i = 0; i < binned_data.size(); i++)
      {
        for (Size j = 0; j < i; j++)
        {
          //distance value is 1-similarity value, since similarity is in range of [0,1]
          original_distance.setValue(i, j, 1 - comparator(binned_data[i], binned_data[j]));
        }
      }

      // create Clustering with ClusterMethod, DistanceMatrix and Data
      clusterer(original_distance, cluster_tree, threshold_);
    }

    /// get the threshold
    double getThreshold()
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
  class OPENMS_DLLAPI UnnormalizedComparator :
    public Exception::BaseException
  {
public:
    UnnormalizedComparator(const char * file, int line, const char * function, const char * message
                             = "Clustering with unnormalized similarity measurement requested, normalized is mandatory") throw();
    ~UnnormalizedComparator() throw() override;
  };

}
#endif //OPENMS_COMPARISON_CLUSTERING_CLUSTERHIERARCHICAL_H
