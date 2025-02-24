// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/config.h>

#include <unordered_map>

#include <map> // for multimap<>
#include <vector> // for vector<>
#include <set> // for set<>
#include <utility> // for pair<>

namespace OpenMS
{
  class GridFeature;

/**
     @brief A representation of a QT cluster used for feature grouping.

     Ultimately, a cluster represents a group of corresponding features (or
     consensus features) from different input maps (feature maps or consensus
     maps).

     Clusters are defined by their center points (one feature each). A cluster
     also stores a number of potential cluster elements (other features) from
     different input maps, together with their distances to the cluster center.
     Every feature that satisfies certain constraints with respect to the
     cluster center is a @e potential cluster element. However, since a feature
     group can only contain one feature from each input map, only the "best"
     (i.e. closest to the cluster center) such feature is considered a true
     cluster element. To save memory, only the "best" element for each map is
     stored inside a cluster.

     The QT clustering algorithm has the characteristic of initially producing
     all possible, overlapping clusters. Iteratively, the best cluster is then
     extracted and the clustering is recomputed for the remaining points.

     In our implementation, multiple rounds of clustering are not necessary.
     Instead, the clustering is updated in each iteration. This is the reason
     for temporarily storing all potential cluster elements: When a certain cluster is
     finalized, its elements have to be removed from the remaining clusters,
     and affected clusters change their composition. (Note that clusters can
     also be invalidated by this, if the cluster center is being removed.)

     The quality of a cluster is the normalized average distance to the cluster
     center for present and missing cluster elements. The distance value for
     missing elements (if the cluster contains no feature from a certain input
     map) is the user-defined threshold that marks the maximum allowed radius
     of a cluster.

     When adding elements to the cluster, the client needs to call
     initializeCluster first and the client needs to call finalizeCluster after
     adding the last element.  After finalizeCluster, the client may not add
     any more elements through the add function (the client must call
     initializeCluster again before adding new elements).

     If use_id_ is set, clusters are extended only with elements that have at least one
     matching ID. Quality is then computed as the best quality of all possible IDs and
     this ID is then used as the only (representative) ID of the cluster. The left-out
     alternative IDs might be added back later based on the original features though.

     @todo This implementation may benefit from two separate implementations (one considering IDs/annotations
        one without). The current implementation most likely hinders speed/memory of both by trying to do both in one.
        The ID-based implementation could additionally benefit from ID scores and make use of ConsensusID functions.

     @see QTClusterFinder

     @ingroup Datastructures
*/
  class OPENMS_DLLAPI QTCluster
  {
public:

    // need to store more than one
    typedef std::multimap<double, const GridFeature*> NeighborList;
    typedef std::unordered_map<Size, NeighborList> NeighborMapMulti;

    struct Neighbor
    {
      double distance;
      const GridFeature* feature;
    };

    typedef std::unordered_map<Size, Neighbor> NeighborMap;

    struct Element
    {
      Size map_index;
      const GridFeature* feature;
    };

    typedef std::vector<Element> Elements;

    /** 
     * @brief Class to store the bulk internal data (neighbors, annotations, etc.)
     * 
     * Has no functionality without a QTCluster pointing to it.
     * Create object of this class before calling constructor of QTCluster
     */
    class OPENMS_DLLAPI BulkData
    {
      friend class QTCluster;

      public:

        /**
         * @brief Detailed constructor of the cluster body
         * @param center_point Pointer to the center point
         * @param num_maps Number of input maps
         * @param max_distance Maximum allowed distance of two points
         * @param x_coord The x-coordinate in the grid cell
         * @param y_coord The y-coordinate in the grid cell
         * @param id Unique ID of this cluster
         */
        BulkData(const OpenMS::GridFeature* const center_point, 
                Size num_maps, double max_distance,
                Int x_coord, Int y_coord, Size id);
      
      private:

        /// Pointer to the cluster center
        const GridFeature* const center_point_;

        /// unique id of this cluster
        Size id_;

        /**
         * @brief Map that keeps track of the best current feature for each map
         *
         */
        NeighborMap neighbors_;

        /**
         * @brief Temporary map tracking *all* neighbors
         *
         * For each input run, a multimap which contains pointers to all
         * neighboring elements and the respective distance.
         *
         */
        NeighborMapMulti tmp_neighbors_;

        /// Maximum distance of a point that can still belong to the cluster
        double max_distance_;

        /// Number of input maps
        Size num_maps_;

        /// x coordinate in the grid cell
        Int x_coord_;

        /// y coordinate in the grid cell
        Int y_coord_;

        /**
         * @brief Set of annotations of the cluster
         *
         * The set of peptide sequences that is compatible to the cluster center
         * and results in the best cluster quality.
         */
        std::set<AASequence> annotations_;
    };

    /**
     * @brief Detailed constructor of the cluster head
     * @param data Pointer to internal data
     * @param use_IDs Use peptide annotations?
     */

    QTCluster(BulkData* const data, bool use_IDs);

    /** 
     * @brief Default constructor not accessible
     * Objects of this class should only exist with a valid BulkData* given.
     * Otherwise most of the member functions are undefined behavior or produce segfaults
     */
    QTCluster() = delete;

    /** 
     * @brief Cheap copy ctor because most of the data lies outside of this class (BulkData*)
     * Be very careful with this copy constructor. The copy will point to the same
     * BulkData object as the given QTCluster. The latter one shouldn't be used anymore.
     * This operation is only allowed because the boost::heap interface needs it.
     */
    QTCluster(const QTCluster& rhs) = default;

    /// Cheap copy assignment, see copy ctor for details
    QTCluster& operator=(const QTCluster& rhs) = default;

    /// cheap move ctor because most of the data lies outside of this class (BulkData*)
    QTCluster(QTCluster&& rhs) = default;

    /// cheap move assignment because most of the data lies outside of this class (BulkData*)
    QTCluster& operator=(QTCluster&& rhs) = default;

    ~QTCluster() = default;

    /// Returns the cluster center
    const GridFeature* getCenterPoint() const;

    /// returns the clusters id
    Size getId() const;

    /// Returns the RT value of the cluster
    double getCenterRT() const;

    /// Returns the m/z value of the cluster center
    double getCenterMZ() const;

    /// Returns the x coordinate in the grid
    Int getXCoord() const;

    /// Returns the y coordinate in the grid
    Int getYCoord() const;

    /// Returns the size of the cluster (number of elements, incl. center)
    Size size() const;

    /// Compare by quality
    bool operator<(const QTCluster& cluster) const;

    /**
     * @brief Adds a new element/neighbor to the cluster
     * @note There is no check whether the element/neighbor already exists in the cluster!
     * @param element The element to be added
     * @param distance Distance of the element to the center point
     */
    void add(const GridFeature* const element, double distance);

    /// Gets the clustered elements meaning neighbors + cluster center
    Elements getElements() const;

    /**
     * @brief Updates the cluster after the indicated data points are removed
     *
     * @param removed The datapoints to be removed from the cluster
     *
     * @return Whether the cluster composition has changed due to the update
     */
    bool update(const Elements& removed);

    /// Returns the cluster quality and recomputes if necessary
    double getQuality();

    /// Returns the cluster quality without recomputing
    double getCurrentQuality() const;

    /// Return the set of peptide sequences annotated to the cluster center
    const std::set<AASequence>& getAnnotations();

    /**
     * @brief Sets current cluster as invalid (also frees some memory)
     *
     * @note Do not attempt to use the cluster again once it is invalid, some
     * internal data structures have now been cleared
     *
     */
    void setInvalid();

    /// Whether current cluster is invalid
    inline bool isInvalid() const
    {
      return !valid_;
    }

    /// Has to be called before adding elements (calling QTCluster::add)
    void initializeCluster();

    /// Has to be called after adding elements (after calling QTCluster::add one or multiple times)
    void finalizeCluster();

    /// Get all current neighbors
    Elements getAllNeighbors() const;

    private:
      /// Computes the quality of the cluster
      void computeQuality_();

      /**
       * @brief Finds the optimal annotation (peptide sequences) for the cluster
       *
       * The optimal annotation is the one that results in the best quality. It
       * is stored in @p annotations_;
       *
       * This function is only needed when peptide ids are used and the current
       * center point does not have any peptide id associated with it. In this
       * case, it is not clear which peptide id the current cluster should use.
       * The function thus iterates through all possible peptide ids and selects
       * the one producing the best cluster.
       *
       * This function needs access to all possible neighbors for this cluster
       * and thus can only be run when tmp_neighbors_ is filled (which is during
       * the filling of a cluster). The function thus cannot be called after
       * finalizing the cluster.
       *
       * @returns The total distance between cluster elements and the center.
       */
      double optimizeAnnotations_();

      /// compute seq table, mapping: peptides -> best distance per input map
      void makeSeqTable_(std::map<AASequence, std::map<Size,double>>& seq_table) const;
      
      /// report elements that are compatible with the optimal annotation
      void recomputeNeighbors_();

      /// Quality of the cluster
      double quality_;

      /// Pointer to data members
      BulkData* data_;

      /// Whether current cluster is valid
      bool valid_;

      /// Has the cluster changed (if yes, quality needs to be recomputed)?
      bool changed_;

      /// Keep track of peptide IDs and use them for matching?
      bool use_IDs_;

      /** 
       * @brief Whether initial collection of all neighbors is needed
       *
       * This variable stores whether we need to collect all annotations first
       * before we can decide upon the best set of cluster points. This is
       * usually only necessary if the center point does not have an annotation
       * but we want to use ids.
       *
      */
      bool collect_annotations_;

      /// Whether current cluster is accepting new elements or not (if true, no more new elements allowed)
      bool finalized_;
  };

} // namespace OpenMS
