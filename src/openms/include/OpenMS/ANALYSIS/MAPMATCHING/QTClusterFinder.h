// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/DATASTRUCTURES/QTCluster.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>

#include <boost/unordered_map.hpp>
#include <boost/heap/fibonacci_heap.hpp>

#include <list>
#include <vector>
#include <unordered_set>
#include <utility> // for pair<>

namespace OpenMS
{

/**
   @brief A variant of QT clustering for the detection of feature groups.

   The algorithm accumulates all features from all input maps, then applies a
   variant of QT clustering to find groups of corresponding features. In more
   detail, every feature from every input map is considered as a potential
   cluster center. For every center, its nearest neighbors from the other input
   maps are detected and added to the potential cluster. Iteratively, the
   cluster with the highest quality is extracted and the clustering is updated.

   <b>Properties affecting the grouping</b>

   To be included in a particular cluster, a feature has to fulfill the following conditions:

   @li differences in RT and m/z from the cluster center must be below
       user-defined thresholds (@p distance_RT:max_difference and @p distance_MZ:max_difference),

   @li the charge state must match that of the cluster center (unless @p
       ignore_charge is set),

   @li if @p use_identifications is set and both the feature and the cluster
       center are annotated with peptide identifications, the identifications have
       to match.

   Every cluster contains at most one feature from each input map - namely the
   feature closest to the cluster center that meets the criteria and does not
   belong to a better cluster.

   The notion of "closeness" for features is defined by the distance function
   implemented in @ref FeatureDistance, the parameters of which can be set by
   the user.

   The quality of a cluster is computed from the number of elements in it and
   their distances to the cluster center. For more details see QTCluster.

   <b>Optimization</b>

   This algorithm includes a number of optimizations to reduce run-time:
   @li two-dimensional hashing of features,
   @li a look-up table for feature distances,
   @li a variant of QT clustering that requires only one round of clustering.

   @see FeatureGroupingAlgorithmQT

   @htmlinclude OpenMS_QTClusterFinder.parameters

   @ingroup FeatureGrouping
*/

  class OPENMS_DLLAPI QTClusterFinder :
    public BaseGroupFinder
  {
  public:

    /// Distances between pairs of grid features
    typedef OpenMSBoost::unordered_map< 
              std::pair<OpenMS::GridFeature*, OpenMS::GridFeature*>,
              double> PairDistances;

    /// Map to store which grid features are next to which clusters (saves the clusters ids)
    typedef OpenMSBoost::unordered_map<
              const OpenMS::GridFeature*, std::unordered_set<Size> > ElementMapping;

    /// Heap to efficiently find the best clusters
    typedef boost::heap::fibonacci_heap<QTCluster> Heap;

    typedef HashGrid<OpenMS::GridFeature*> Grid;

  private:
    /// Number of input maps
    Size num_maps_;

    /// Consider peptide identifications for grouping?
    bool use_IDs_;

    /// Maximum RT difference
    double max_diff_rt_;

    /// Maximum m/z difference
    double max_diff_mz_;

    /// Maximum m/z difference
    int nr_partitions_;

    /// Feature distance functor
    FeatureDistance feature_distance_;

    /// Set of features already used
    std::unordered_set<const OpenMS::GridFeature*> already_used_;

    /// set of used cluster ids during one run of run_internal_
    std::unordered_set<Size> deleted_ids_;

    /**
       @brief Calculates the distance between two grid features.
    */
    double getDistance_(const OpenMS::GridFeature* left, const
        OpenMS::GridFeature* right);

    /// Sets algorithm parameters
    void setParameters_(double max_intensity, double max_mz);

    /**
     * @brief Extract the best cluster from cluster_heads and turn it into a consensus feature
     * 
     * @param cluster_heads the heap where the clusters are stored, must not be empty
     * @param[out] feature The resulting consensus feature which is constructed here
     * @param element_mapping the element mapping is used to update clusters when features are removed
     * @param grid the grid is used to find new features for clusters that have to be updated
     * @param handles used to access clusters if we know their id from the element mapping
     * 
     * @return bool whether a consensus feature was made or not
     */
    bool makeConsensusFeature_(Heap& cluster_heads,
                               ConsensusFeature& feature,
                               ElementMapping& element_mapping,
                               const Grid& grid,
                               const std::vector<Heap::handle_type>& handles);

    /**
     * @brief Computes an initial QT clustering of the points in the hash grid
     * 
     * @param grid the grid is used to find new features for clusters that have to be updated
     * @param cluster_heads the heap where the QTClusters are inserted
     * @param cluster_data vector where the BulkData objects of the clusters are stored
     * @param handles vector where handles of the inserted clusters are stored
     * @param element_mapping the element mapping where all the clusters get registered for their features
     */
    void computeClustering_(const Grid& grid,
                            Heap& cluster_heads,
                            std::vector<QTCluster::BulkData>& cluster_data,
                            std::vector<Heap::handle_type>& handles,
                            ElementMapping& element_mapping);

    /** 
     * @brief Removes id of current top cluster in the heap from element mapping
     * 
     * @param cluster the current top cluster in the heap which id is removed
     * @param element_mapping the element mapping from which the ids are removed
     */
    void removeFromElementMapping_(const QTCluster& cluster,
                                   ElementMapping& element_mapping);
    
    /** 
     * @brief creates a consensus feature from the given elements
     * 
     * @param[out] feature The resulting consensus feature which is constructed here
     * @param quality the quality of the new consensus feature
     * @param elements the original features that from the new consensus feature
     * 
     * @note the features from elements are registered in the already_used_ member of this class
     */
    void createConsensusFeature_(ConsensusFeature& feature, const double quality, 
                                 const QTCluster::Elements& elements);

    /** 
     * @brief update the clustering:
     *
     * 1. remove current best cluster from the heap
     * 2. update all clusters accordingly by removing neighbors used by the current best
     * 3. invalidate clusters whose center has been used by the current best
     * 
     * @param element_mapping the element mapping is used to update clusters and updated itself
     * @param grid the grid is used to find new features for clusters that have to be updated
     * @param cluster_heads the heap is updated for changing clusters and popped in the end
     * @param elements the features that now have to be removed from other clusters than the current best
     * @param handles used to access clusters if we know their id from the element mapping
     * @param best_id id of the current best cluster, will be removed from the element mapping
     * 
     * @note The feature from elements are not deleted from the element mapping. 
     * After this function is called we dont't have any cluster with those features left and
     * therefore don't have to delete them.
     */
    void updateClustering_(ElementMapping& element_mapping,
                           const Grid& grid, 
                           const QTCluster::Elements& elements,
                           Heap& cluster_heads,
                           const std::vector<Heap::handle_type>& handles,
                           Size best_id);

    /// Runs the algorithm on feature maps or consensus maps
    template <typename MapType>
    void run_(const std::vector<MapType>& input_maps, ConsensusMap& result_map);

    /// Runs the algorithm on feature maps or consensus maps (internal)
    template <typename MapType>
    void run_internal_(const std::vector<MapType>& input_maps,
                       ConsensusMap& result_map, bool do_progress);

    /**
     * @brief Adds elements to the cluster based on the elements hashed in the grid
     * 
     * @param grid the grid is used to find neighboring features the cluster
     * @param cluster cluster to which the new elements are added
     */ 
    void addClusterElements_(const Grid& grid, QTCluster& cluster);

protected:

    enum
    {
      RT = Peak2D::RT,
      MZ = Peak2D::MZ
    };

public:

    /// Constructor
    QTClusterFinder();

    /// Destructor
    ~QTClusterFinder() override;

    /// Returns the name of the product
    static const String getProductName()
    {
      return "qt";
    }

    /**
       @brief Runs the algorithm on consensus maps

       @pre The data ranges of the input maps have to be up-to-date (use ConsensusMap::updateRanges).

       @exception Exception::IllegalArgument is thrown if the input data is not valid.
    */
    void run(const std::vector<ConsensusMap>& input_maps,
             ConsensusMap& result_map) override;

    /**
       @brief Runs the algorithm on feature maps

       @pre The data ranges of the input maps have to be up-to-date (use FeatureMap::updateRanges).

       @exception Exception::IllegalArgument is thrown if the input data is not valid.
    */
    void run(const std::vector<FeatureMap>& input_maps,
             ConsensusMap& result_map);

    /// Returns an instance of this class
    static BaseGroupFinder* create()
    {
      return new QTClusterFinder();
    }
  };
} // namespace OpenMS

