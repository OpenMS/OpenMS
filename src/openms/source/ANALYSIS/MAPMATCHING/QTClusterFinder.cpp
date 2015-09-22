// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace std;

namespace OpenMS
{

  QTClusterFinder::QTClusterFinder() :
    BaseGroupFinder(), feature_distance_(FeatureDistance())
  {
    setName(getProductName());

    defaults_.setValue("use_identifications", "false", "Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).");
    defaults_.setValidStrings("use_identifications", ListUtils::create<String>("true,false"));

    defaults_.insert("", feature_distance_.getDefaults());

    defaultsToParam_();
  }

  void QTClusterFinder::setParameters_(double max_intensity,
                                       double max_mz)
  {
    // don't check for low max. intensity, because intensities may be ignored:
    if ((max_mz < 1e-16) || (max_mz > 1e16) || (max_intensity > 1e16))
    {
      String msg = "Maximum m/z or intensity out of range (m/z: " + 
        String(max_mz) + ", intensity: " + String(max_intensity) + "). "
        "Has 'updateRanges' been called on the input maps?";
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       msg);
    }
    use_IDs_ = String(param_.getValue("use_identifications")) == "true";
    max_diff_rt_ = param_.getValue("distance_RT:max_difference");
    max_diff_mz_ = param_.getValue("distance_MZ:max_difference");
    // compute m/z tolerance in Da (if given in ppm; for the hash grid):
    if (param_.getValue("distance_MZ:unit") == "ppm")
    {
      max_diff_mz_ *= max_mz * 1e-6;
    }
    Param distance_params = param_.copy("");
    distance_params.remove("use_identifications");
    feature_distance_ = FeatureDistance(max_intensity, true);
    feature_distance_.setParameters(distance_params);
  }

  template <typename MapType>
  void QTClusterFinder::run_(const vector<MapType>& input_maps,
                             ConsensusMap& result_map)
  {
    num_maps_ = input_maps.size();
    if (num_maps_ < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "At least two input maps required");
    }

    // set up the distance functor (and set other parameters):
    double max_intensity = 0.0;
    double max_mz = 0.0;
    for (typename vector<MapType>::const_iterator map_it = input_maps.begin(); 
         map_it != input_maps.end(); ++map_it)
    {
      max_intensity = max(max_intensity, map_it->getMaxInt());
      max_mz = max(max_mz, map_it->getMax()[0]);
    }

    setParameters_(max_intensity, max_mz);

    // create the hash grid and fill it with features:
    //cout << "Hashing..." << endl;
    list<OpenMS::GridFeature> grid_features;
    Grid grid(Grid::ClusterCenter(max_diff_rt_, max_diff_mz_));
    for (Size map_index = 0; map_index < num_maps_; ++map_index)
    {
      for (Size feature_index = 0; feature_index < input_maps[map_index].size();
           ++feature_index)
      {
        grid_features.push_back(
          GridFeature(input_maps[map_index][feature_index], map_index, 
                      feature_index));
        GridFeature& gfeat = grid_features.back();
        // sort peptide hits once now, instead of multiple times later:
        BaseFeature& bfeat = const_cast<BaseFeature&>(gfeat.getFeature());
        for (vector<PeptideIdentification>::iterator pep_it =
               bfeat.getPeptideIdentifications().begin(); pep_it !=
               bfeat.getPeptideIdentifications().end(); ++pep_it)
        {
          pep_it->sort();
        }
        grid.insert(make_pair(Grid::ClusterCenter(gfeat.getRT(), gfeat.getMZ()),
                              &gfeat));
      }
    }

    // compute QT clustering:
    //cout << "Clustering..." << endl;
    list<QTCluster> clustering;
    computeClustering_(grid, clustering);
    // number of clusters == number of data points:
    Size size = clustering.size();

    // create a temp. map storing which grid features are next to which clusters
    ElementMapping element_mapping;
    for (list<QTCluster>::iterator it = clustering.begin();
         it != clustering.end(); ++it)
    {
      typedef OpenMSBoost::unordered_map<Size, std::vector<GridFeature*> > NeighborMap;
      OpenMSBoost::unordered_map<Size, std::vector<GridFeature*> > neigh = it->getAllNeighbors();
      for (NeighborMap::iterator n_it = neigh.begin(); n_it != neigh.end(); ++n_it)
      {
        for (std::vector<GridFeature*>::iterator i_it = n_it->second.begin(); i_it != n_it->second.end(); ++i_it)
        {
          // remember for each feature (gridfeature) all the cluster elements it belongs to
          element_mapping[*i_it].push_back(&(*it));
        }
      }
    }

    // ensure that all cluster centers are in the list
    for (list<QTCluster>::iterator it = clustering.begin();
         it != clustering.end(); ++it)
    {
      OpenMS::GridFeature* center_feature = it->getCenterPoint();
      element_mapping[center_feature].push_back(&(*it));
    }

    ProgressLogger logger;
    logger.setLogType(ProgressLogger::CMD);
    logger.startProgress(0, size, "linking features");
    Size progress = 0;
    result_map.clear(false);

    while (!clustering.empty())
    {
      // cout << "Clusters: " << clustering.size() << endl;
      ConsensusFeature consensus_feature;
      makeConsensusFeature_(clustering, consensus_feature, element_mapping, grid);
      if (!clustering.empty())
      {
        result_map.push_back(consensus_feature);
      }
      logger.setProgress(progress++);
    }

    logger.endProgress();
  }

  void QTClusterFinder::makeConsensusFeature_(list<QTCluster>& clustering,
                                              ConsensusFeature& feature,
                                              ElementMapping& element_mapping,
                                              Grid& grid)
  {
    // find the best cluster (a valid cluster with the highest score)
    list<QTCluster>::iterator best = clustering.begin();
    while (best != clustering.end() && best->isInvalid())
    {
      ++best;
    }
    for (list<QTCluster>::iterator it = best;
         it != clustering.end(); ++it)
    {
      if (!it->isInvalid())
      {
        if (it->getQuality() > best->getQuality())
        {
          best = it;
        }
      }
    }

    // no more clusters to process -> clear clustering and return
    if (best == clustering.end())
    {
      clustering.clear();
      return;
    }

    OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*> elements;
    best->getElements(elements);
    // cout << "Elements: " << elements.size() << " with best "
    //      << best->getQuality() << " invalid " << best->isInvalid() << endl;

    // create consensus feature from best cluster:
    feature.setQuality(best->getQuality());
    for (OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*>::const_iterator it = 
           elements.begin(); it != elements.end(); ++it)
    {
      feature.insert(it->first, it->second->getFeature());
    }
    feature.computeConsensus();

    // update the clustering:
    // 1. remove current "best" cluster
    // 2. update all clusters accordingly and invalidate elements whose central
    //    element is removed
    best->setInvalid();
    std::vector< QTCluster* > worklist; 

    // Store the id of already used features
    for (OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*>::const_iterator it = 
           elements.begin(); it != elements.end(); ++it)
    {
      already_used_.insert(it->second);
    }

    for (OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*>::const_iterator it = 
           elements.begin(); it != elements.end(); ++it)
    {
      // Identify all features that could potentially have been touched by this
      //  Get all clusters that may potentially need updating

      ElementMapping tmp_element_mapping; // modify copy, then update

      for (std::vector<QTCluster*>::iterator
           cluster  = element_mapping[&(*it->second)].begin();
           cluster != element_mapping[&(*it->second)].end(); ++cluster)
      {

        // we do not want to update invalid features (saves time and does not
        // recompute the quality)
        if (!(*cluster)->isInvalid())
        {

          // remove the elements of the new feature from the cluster
          if ((*cluster)->update(elements))
          {
            // If update returns true, it means that at least one element was
            // removed from the cluster and we need to update that cluster

            // Get the coordinates of the current cluster
            const Int x = (*cluster)->x_coord_; 
            const Int y = (*cluster)->y_coord_;

            (*cluster)->initializeCluster();

            ////////////////////////////////////////
            // Step 1: Iterate through all neighboring grid features and try to
            // add elements to the current cluster to replace the ones we just
            // removed
            OpenMS::GridFeature* center_feature = (*cluster)->getCenterPoint();
            // iterate over neighboring grid cells (1st dimension):
            for (int i = x - 1; i <= x + 1; ++i)
            {
              // iterate over neighboring grid cells (2nd dimension):
              for (int j = y - 1; j <= y + 1; ++j)
              {
                try
                {
                  const Grid::CellContent& act_pos = grid.grid_at(Grid::CellIndex(i, j));

                  for (Grid::const_cell_iterator it_cell = act_pos.begin();
                       it_cell != act_pos.end(); ++it_cell)
                  {

                    OpenMS::GridFeature* neighbor_feature = it_cell->second;

                    // Skip features that we have already used -> we cannot add them to be neighbors any more
                    // remember the ones we already used, we need to skip those 
                    if (already_used_.find(neighbor_feature) != already_used_.end() )
                    {
                      continue;
                    }

                    // consider only "real" neighbors, not the element itself:
                    if (center_feature != neighbor_feature)
                    {
                      // NOTE: this actually caches the distance -> memory problem
                      double dist = getDistance_(center_feature, neighbor_feature);

                      if (dist == FeatureDistance::infinity)
                      {
                        continue; // conditions not satisfied
                      }
                      // if neighbor point is a possible cluster point, add it:
                      if (!use_IDs_ || compatibleIDs_(*(*cluster), neighbor_feature))
                      {
                        (*cluster)->add(neighbor_feature, dist);
                      }
                    }
                  }
                }
                catch (std::out_of_range&)
                {
                }
              }
            }

            (*cluster)->finalizeCluster();

            ////////////////////////////////////////
            // Step 2: update element_mapping as the best feature for each
            // cluster may have changed
            typedef OpenMSBoost::unordered_map<Size, std::vector<GridFeature*> > NeighborMap;
            OpenMSBoost::unordered_map<Size, std::vector<GridFeature*> > neigh = (*cluster)->getAllNeighbors();
            for (NeighborMap::iterator n_it = neigh.begin(); n_it != neigh.end(); ++n_it)
            {
              for (std::vector<GridFeature*>::iterator i_it = n_it->second.begin(); i_it != n_it->second.end(); ++i_it)
              {
                // remember for each feature (gridfeature) all the cluster elements it belongs to
                tmp_element_mapping[*i_it].push_back(*cluster);
              }
            }

          }
        }
      }

      distances_.clear(); // clear memory 
      for (ElementMapping::iterator it = tmp_element_mapping.begin(); it != tmp_element_mapping.end(); it++ )
      {
        for(std::vector<QTCluster*>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++)
        {
          element_mapping[ it->first ].push_back( (*it2) );
        }
      }
    }

  }

  void QTClusterFinder::run(const vector<ConsensusMap>& input_maps,
                            ConsensusMap& result_map)
  {
    run_(input_maps, result_map);
  }

  void QTClusterFinder::run(const std::vector<FeatureMap>& input_maps,
                            ConsensusMap& result_map)
  {
    run_(input_maps, result_map);
  }

  void QTClusterFinder::computeClustering_(Grid& grid,
                                           list<QTCluster>& clustering)
  {
    clustering.clear();
    distances_.clear();
    // FeatureDistance produces normalized distances (between 0 and 1):
    const double max_distance = 1.0;

    // iterate over all grid cells:
    for (Grid::iterator it = grid.begin(); it != grid.end(); ++it)
    {
      const Grid::CellIndex& act_coords = it.index();
      const Int x = act_coords[0], y = act_coords[1];
      //cout << x << " " << y << endl;

      OpenMS::GridFeature* center_feature = it->second;
      QTCluster cluster(center_feature, num_maps_, max_distance, use_IDs_, x, y);

      cluster.initializeCluster();

      // iterate over neighboring grid cells (1st dimension):
      for (int i = x - 1; i <= x + 1; ++i)
      {
        // iterate over neighboring grid cells (2nd dimension):
        for (int j = y - 1; j <= y + 1; ++j)
        {
          try
          {
            const Grid::CellContent& act_pos = grid.grid_at(Grid::CellIndex(i, j));

            for (Grid::const_cell_iterator it_cell = act_pos.begin();
                 it_cell != act_pos.end(); ++it_cell)
            {
              OpenMS::GridFeature* neighbor_feature = it_cell->second;
              // consider only "real" neighbors, not the element itself:
              if (center_feature != neighbor_feature)
              {
                // NOTE: this actually caches the distance in memory (potential memory issue)
                double dist = getDistance_(center_feature, neighbor_feature);
                if (dist == FeatureDistance::infinity)
                {
                  continue; // conditions not satisfied
                }
                // if neighbor point is a possible cluster point, add it:
                if (!use_IDs_ || compatibleIDs_(cluster, neighbor_feature))
                {
                  cluster.add(neighbor_feature, dist);
                }
              }
            }
          }
          catch (std::out_of_range&)
          {
          }
        }
      }

      cluster.finalizeCluster();
      distances_.clear(); // reduces memory from 550 MB to ca 350 MB
      clustering.push_back(cluster);
    }

  }

  double QTClusterFinder::getDistance_(OpenMS::GridFeature* left,
                                       OpenMS::GridFeature* right)
  {
    // look-up in the distance map:
    const pair<OpenMS::GridFeature*, OpenMS::GridFeature*> key = make_pair(min(left, right),
                                                           max(left, right));
    PairDistances::const_iterator pos = distances_.find(key);
    if (pos != distances_.end()) // distance computed before
    {
      return pos->second;
    }
    else // compute distance now and store it for later
    {
      double dist = feature_distance_(left->getFeature(), 
                                      right->getFeature()).second;
      distances_[key] = dist;
      return dist;
    }
  }

  bool QTClusterFinder::compatibleIDs_(QTCluster& cluster,
                                       const OpenMS::GridFeature* neighbor)
  {
    if (cluster.getAnnotations().empty())
      return true;

    if (neighbor->getAnnotations().empty())
      return true;

    return cluster.getAnnotations() == neighbor->getAnnotations();
  }

  QTClusterFinder::~QTClusterFinder()
  {
  }

} // namespace OpenMS
