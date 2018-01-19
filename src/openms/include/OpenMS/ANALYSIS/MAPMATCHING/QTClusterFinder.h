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
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_QTCLUSTERFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_QTCLUSTERFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseGroupFinder.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/COMPARISON/CLUSTERING/HashGrid.h>
#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/DATASTRUCTURES/QTCluster.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureDistance.h>

#include <boost/unordered_map.hpp>

#include <list>
#include <vector>
#include <set>
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
private:

    /// Distances between pairs of grid features
    typedef OpenMSBoost::unordered_map< 
              std::pair<OpenMS::GridFeature*, OpenMS::GridFeature*>,
              double> PairDistances;

    /// Map to store which grid features are next to which clusters
    typedef OpenMSBoost::unordered_map<
              OpenMS::GridFeature*, std::vector<QTCluster*> > ElementMapping;

    typedef HashGrid<OpenMS::GridFeature*> Grid;

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
    std::set<OpenMS::GridFeature*> already_used_;

    /**
       @brief Calculates the distance between two grid features.
    */
    double getDistance_(const OpenMS::GridFeature* left, const
        OpenMS::GridFeature* right);

    /// Sets algorithm parameters
    void setParameters_(double max_intensity, double max_mz);

    /// Generates a consensus feature from the best cluster and updates the clustering
    void makeConsensusFeature_(std::list<QTCluster>& clustering,
                               ConsensusFeature& feature,
                               ElementMapping& element_mapping, Grid&);

    /// Computes an initial QT clustering of the points in the hash grid
    void computeClustering_(Grid& grid, std::list<QTCluster>& clustering);

    /// Runs the algorithm on feature maps or consensus maps
    template <typename MapType>
    void run_(const std::vector<MapType>& input_maps, ConsensusMap& result_map);

    /// Runs the algorithm on feature maps or consensus maps (internal)
    template <typename MapType>
    void run_internal_(const std::vector<MapType>& input_maps,
                       ConsensusMap& result_map, bool do_progress);

    /// Adds elements to the cluster based on the elements hashed in the grid
    void addClusterElements_(int x, int y, const Grid& grid, QTCluster& cluster,
      const OpenMS::GridFeature* center_feature);

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

#endif /* OPENMS_ANALYSIS_MAPMATCHING_QTCLUSTERFINDER_H */
