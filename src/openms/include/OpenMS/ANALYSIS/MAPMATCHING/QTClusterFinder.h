// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/ML/CLUSTERING/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/DATASTRUCTURES/QTCluster.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>

#include <boost/heap/fibonacci_heap.hpp>
#include <unordered_map>

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
    typedef std::unordered_map< 
              std::pair<OpenMS::GridFeature*, OpenMS::GridFeature*>,
              double> PairDistances;

    /// Map to store which grid features are next to which clusters (saves the clusters ids)
    typedef std::unordered_map<
              const OpenMS::GridFeature*, std::unordered_set<Size> > ElementMapping;

    /// Heap to efficiently find the best clusters
    typedef boost::heap::fibonacci_heap<QTCluster> Heap;

    typedef HashGrid<OpenMS::GridFeature*> Grid;

  private:
    /// Number of input maps
    Size num_maps_;

    /// Consider peptide identifications for grouping?
    bool use_IDs_;

    //TODO we could also bin by equal sized RT bins, or by a fixed RT size
    //TODO could be made dependent on nr. of maps (e.g. with 5 maps you get [around] 4 diffs per ID already)
    /// Min. nr. of differences from matched IDs requested to calculate a linking tolerance per RT bin
    Size min_nr_diffs_per_bin_;

    /// Min. score for an ID to be considered for tolerance estimation
    double min_score_;

    /// Distance penalty for unidentified features when finding best neighbor per map and for cluster quality calculation
    /// and therefore the order in which they are popped from the heap.
    /// Since distances are normalized, a penalty of 1.0 will always prefer IDed features.
    double noID_penalty_;

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

    /// Map of median RTs to allowed linking tolerances (on the same RT scale) for unIDed features.
    /// This should be interpreted as bins from the current median RT to the next.
    std::map<double, double> bin_tolerances_;

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
     * After this function is called we don't have any cluster with those features left and
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

    /**
     * @brief Looks up the matching bin for @p rt in bin_tolerances_ and checks if @p dist is in the allowed range.
     */
    bool distIsOutlier_(double dist, double rt);

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

  };
} // namespace OpenMS

