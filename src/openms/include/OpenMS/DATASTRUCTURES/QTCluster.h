// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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


#ifndef OPENMS_DATASTRUCTURES_QTCLUSTER_H
#define OPENMS_DATASTRUCTURES_QTCLUSTER_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>
#include <OpenMS/config.h>

#include <boost/unordered_map.hpp>

namespace OpenMS
{
  class GridFeature;
  class AASequence;
  
  // Boost switch since with 1.47 several classes got moved into a new
  // boost::unordered namespace (specifically unordered_map).
  namespace OpenMSBoost
  {
#if OPENMS_BOOST_VERSION_MINOR > 47
      using namespace boost::unordered;
#else
      using namespace boost;
#endif
  }


/**
     @brief A representation of a QT cluster used for feature grouping.

     Ultimately, a cluster represents a group of corresponding features (or consensus features) from different input maps (feature maps or consensus maps).

     Clusters are defined by their center points (one feature each). A cluster also stores a number of potential cluster elements (other features) from different input maps, together with their distances to the cluster center.
     Every feature that satisfies certain constraints with respect to the cluster center is a @e potential cluster element. However, since a feature group can only contain one feature from each input map, only the "best" (i.e. closest to the cluster center) such feature is considered a true cluster element.

     The QT clustering algorithm has the characteristic of initially producing all possible, overlapping clusters. Iteratively, the best cluster is then extracted and the clustering is recomputed for the remaining points.

     In our implementation, multiple rounds of clustering are not necessary. Instead, the clustering is updated in each iteration. This is the reason for storing all potential cluster elements: When a certain cluster is finalized, its elements have to be removed from the remaining clusters, and affected clusters change their composition. (Note that clusters can also be invalidated by this, if the cluster center is being removed.)

     The quality of a cluster is the normalized average distance to the cluster center for present and missing cluster elements. The distance value for missing elements (if the cluster contains no feature from a certain input map) is the user-defined threshold that marks the maximum allowed radius of a cluster.

     @see QTClusterFinder

     @ingroup Datastructures
*/

  class OPENMS_DLLAPI QTCluster
  {
private:

    /**
     * @brief Mapping: input map -> distance to center (ordered!) -> neighboring point
     * @note There should never be an empty sub-map! (When a sub-map becomes empty, it should be removed from the overall map.)
     */
    typedef OpenMSBoost::unordered_map<Size, std::multimap<double, GridFeature *> > NeighborMap;

    /// Pointer to the cluster center
    GridFeature * center_point_;

    /**
     * @brief Neighbors of the cluster center, sorted by distance, for different input maps.
     *
     * The first (best) point in each sub-map is considered a cluster element.
     */
    NeighborMap neighbors_;

    /// Maximum distance of a point that can still belong to the cluster
    double max_distance_;

    /// Number of input maps
    Size num_maps_;

    /// Quality of the cluster
    double quality_;

    /// Has the cluster changed (if yes, quality needs to be recomputed)?
    bool changed_;

    /// Keep track of peptide IDs and use them for matching?
    bool use_IDs_;

    /**
     * @brief Set of annotations of the cluster
     *
     * The set of peptide sequences that is compatible to the cluster center and results in the best cluster quality.
     */
    std::set<AASequence> annotations_;

    /// Base constructor (not accessible)
    QTCluster();

    /// Computes the quality of the cluster
    void computeQuality_();

    /**
     * @brief Finds the optimal annotation (peptide sequences) for the cluster
     *
     * The optimal annotation is the one that results in the best quality. It is stored in @p annotations_;
     *
     * @returns The total distance between cluster elements and the center.
     */
    double optimizeAnnotations_();

    bool valid_;
public:
    /**
     * @brief Detailed constructor
     * @param center_point Pointer to the center point
     * @param num_maps Number of input maps
     * @param max_distance Maximum allowed distance of two points
     * @param use_IDs Use peptide annotations?
     */
    QTCluster(GridFeature * center_point, Size num_maps,
              double max_distance, bool use_IDs);

    /// Destructor
    virtual ~QTCluster();

    /// Returns the RT value of the cluster
    double getCenterRT() const;

    /// Returns the m/z value of the cluster center
    double getCenterMZ() const;

    /// Returns the size of the cluster (number of elements, incl. center)
    Size size() const;

    /// Compare by quality
    bool operator<(QTCluster & cluster);

    /**
     * @brief Adds a new element/neighbor to the cluster
     * @note There is no check whether the element/neighbor already exists in the cluster!
     * @param element The element to be added
     * @param distance Distance of the element to the center point
     */
    void add(GridFeature * element, double distance);

    /// Gets the clustered elements
    void getElements(OpenMSBoost::unordered_map<Size, GridFeature *> & elements);

    /**
     * @brief Updates the cluster after data points were removed
     * @return Whether the cluster is still valid (it's not if the cluster center is among the removed points).
     */
    bool update(const OpenMSBoost::unordered_map<Size, GridFeature *> & removed);

    /// Returns the cluster quality
    double getQuality();

    /// Return the set of peptide sequences annotated to the cluster center
    const std::set<AASequence> & getAnnotations();

    inline void setInvalid() {valid_ = false;}

    inline bool isInvalid() {return !valid_;}

    NeighborMap getNeighbors() {return neighbors_;}

  };
}

#endif // OPENMS_DATASTRUCTURES_QTCLUSTER_H
