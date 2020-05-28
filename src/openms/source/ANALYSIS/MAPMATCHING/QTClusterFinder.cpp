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

#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/FeatureHandle.h>

// #define DEBUG_QTCLUSTERFINDER

using std::list;
using std::vector;
using std::max;
using std::make_pair;
using std::unordered_set;


namespace OpenMS
{
  QTClusterFinder::QTClusterFinder() :
    BaseGroupFinder(), feature_distance_(FeatureDistance())
  {
    setName(getProductName());

    defaults_.setValue("use_identifications", "false", "Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).");
    defaults_.setValidStrings("use_identifications", ListUtils::create<String>("true,false"));
    defaults_.setValue("nr_partitions", 100, "How many partitions in m/z space should be used for the algorithm (more partitions means faster runtime and more memory efficient execution )");
    defaults_.setMinInt("nr_partitions", 1);


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
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       msg);
    }
    use_IDs_ = String(param_.getValue("use_identifications")) == "true";
    nr_partitions_ = param_.getValue("nr_partitions");
    max_diff_rt_ = param_.getValue("distance_RT:max_difference");
    max_diff_mz_ = param_.getValue("distance_MZ:max_difference");
    // compute m/z tolerance in Da (if given in ppm; for the hash grid):
    if (param_.getValue("distance_MZ:unit") == "ppm")
    {
      max_diff_mz_ *= max_mz * 1e-6;
    }
    Param distance_params = param_.copy("");
    distance_params.remove("use_identifications");
    distance_params.remove("nr_partitions");
    feature_distance_ = FeatureDistance(max_intensity, true);
    feature_distance_.setParameters(distance_params);
  }

  template <typename MapType>
  void QTClusterFinder::run_(const vector<MapType>& input_maps,
                             ConsensusMap& result_map)
  {
    // update parameters (dummy)
    setParameters_(1, 1);

    result_map.clear(false);

    std::vector< double > massrange; 
    for (typename vector<MapType>::const_iterator map_it = input_maps.begin(); 
         map_it != input_maps.end(); ++map_it)
    {
      for (typename MapType::const_iterator feat_it = map_it->begin();
          feat_it != map_it->end(); ++feat_it)
      {
        massrange.push_back(feat_it->getMZ());
      }
    }
    std::sort(massrange.begin(), massrange.end());

    if (nr_partitions_ == 1)
    {
      // Only one partition 
      run_internal_(input_maps, result_map, true);
    }
    else
    {
      // partition at boundaries -> this should be safe because there cannot be
      // any cluster reaching across boundaries

      // minimal differences between two m/z values 
      double massrange_diff = max_diff_mz_;
      int pts_per_partition = massrange.size() / nr_partitions_;

      // if m/z tolerance is specified in ppm, we adapt massrange_diff
      // in each iteration below
      bool mz_ppm = param_.getValue("distance_MZ:unit") == "ppm";
      double mz_tol = param_.getValue("distance_MZ:max_difference");

      // compute partition boundaries
      std::vector< double > partition_boundaries; 
      partition_boundaries.push_back(massrange.front());
      for (size_t j = 0; j < massrange.size()-1; j++)
      {
        if (mz_ppm)
        {
          massrange_diff = mz_tol * 1e-6 * massrange[j+1];
        }
        if (fabs(massrange[j] - massrange[j+1]) > massrange_diff)
        {
          if (j >= (partition_boundaries.size() ) * pts_per_partition  )
          {
            partition_boundaries.push_back((massrange[j] + massrange[j+1])/2.0);
          }
        }
      }
      // add last partition (a bit more since we use "smaller than" below)
      partition_boundaries.push_back(massrange.back() + 1.0);

      ProgressLogger logger;
      Size progress = 0;
      logger.setLogType(ProgressLogger::CMD);
      logger.startProgress(0, partition_boundaries.size(), "Linking features");
      for (size_t j = 0; j < partition_boundaries.size()-1; j++)
      {
        double partition_start = partition_boundaries[j];
        double partition_end = partition_boundaries[j+1];

        std::vector<MapType> tmp_input_maps(input_maps.size());
        for (size_t k = 0; k < input_maps.size(); k++)
        {
          // iterate over all features in the current input map and append
          // matching features (within the current partition) to the temporary
          // map
          for (size_t m = 0; m < input_maps[k].size(); m++)
          {
            if (input_maps[k][m].getMZ() >= partition_start && 
                input_maps[k][m].getMZ() < partition_end)
            {
              tmp_input_maps[k].push_back(input_maps[k][m]);
            }
          }
          tmp_input_maps[k].updateRanges();
        }

        // run algo on current partition
        run_internal_(tmp_input_maps, result_map, false);
        logger.setProgress(progress++);
      }

      logger.endProgress();
    }
  }

  template <typename MapType>
  void QTClusterFinder::run_internal_(const vector<MapType>& input_maps,
                             ConsensusMap& result_map, bool do_progress)
  {
    // clear temporary data structures
    already_used_.clear();

    num_maps_ = input_maps.size();
    if (num_maps_ < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "At least two input maps required");
    }

    // set up the distance functor (and set other parameters)
    // for the current partition
    double max_intensity = 0.0;
    double max_mz = 0.0;
    for (typename vector<MapType>::const_iterator map_it = input_maps.begin(); 
         map_it != input_maps.end(); ++map_it)
    {
      max_intensity = max(max_intensity, map_it->getMaxInt());
      max_mz = max(max_mz, map_it->getMax().getY());
    }
    setParameters_(max_intensity, max_mz);

    // create the hash grid and fill it with features:
    // std::cout << "Hashing..." << std::endl;
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
    // std::cout << "Clustering..." << std::endl;

    // "hot" cluster heads, we can extract the best efficiently 
    Heap cluster_heads;

    // handles to cluster heads to reach them (index == cluster.id_) in cluster_heads for updating
    vector<Heap::handle_type> handles;

    // "cold" cluster bodies, where most of their data lies
    vector<QTCluster::BulkData> cluster_data;

    // map to get ids from clusters, who contain a certain grid feature
    ElementMapping element_mapping;

    computeClustering_(grid, cluster_heads, cluster_data, handles, element_mapping);

    // number of clusters == number of data points:
    Size size = cluster_heads.size();

    ProgressLogger logger;
    Size progress = 0;
    if (do_progress)
    {
      logger.setLogType(ProgressLogger::CMD);
      logger.startProgress(0, size, "Linking features");
    }

    while (!cluster_heads.empty())
    {
      // std::cout << "Clusters: " << clustering.size() << std::endl;

      ConsensusFeature consensus_feature;
      bool made_feature = makeConsensusFeature_(cluster_heads, consensus_feature, 
                                                element_mapping, grid, handles);

      if (made_feature)
      {
        result_map.push_back(consensus_feature);
      }
      if (do_progress) logger.setProgress(progress++);
    }

    if (do_progress) logger.endProgress();

    // clear the deleted ids for the next run_internal_ call
    deleted_ids_.clear();
  }

  bool QTClusterFinder::makeConsensusFeature_(Heap& cluster_heads,
                                              ConsensusFeature& feature,
                                              ElementMapping& element_mapping,
                                              const Grid& grid,
                                              const vector<Heap::handle_type>& handles)
  {
    // pop until the top is valid
    while (cluster_heads.top().isInvalid())
    {
      removeFromElementMapping_(cluster_heads.top(), element_mapping);

      // remember that we deleted this id
      deleted_ids_.insert(cluster_heads.top().getId());

      cluster_heads.pop();

      // if the last remaining cluster was invalid, no consensus feature is created
      if (cluster_heads.empty()) return false;
    }

    const QTCluster& best = cluster_heads.top();

    QTCluster::Elements const elements = best.getElements();

#ifdef DEBUG_QTCLUSTERFINDER
    std::cout << "Elements: " << elements.size() << " with best "
         << best->getQuality() << " invalid " << best->isInvalid() << std::endl;
#endif

    createConsensusFeature_(feature, best.getCurrentQuality(), elements);

#ifdef DEBUG_QTCLUSTERFINDER
    std::cout << " create new consensus feature " << feature.getRT() << " " << feature.getMZ() << " from " << best->getCenterPoint()->getFeature().getUniqueId() << std::endl;
    for (OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*>::const_iterator
         it = elements.begin(); it != elements.end(); ++it)
    {
      std::cout << "   = element id : " << it->second->getFeature().getUniqueId() << std::endl;
    }
#endif

    updateClustering_(element_mapping, grid, elements, cluster_heads, handles, best.getId());

    // made a consensus feature
    return true;
  }

  void QTClusterFinder::removeFromElementMapping_(const QTCluster& cluster,
                                                  ElementMapping& element_mapping)
  {
    /* We have to erase all references to this cluster from the element mapping
     * before it is popped from the heap and deleted.
     * This function is called in makeConsensusFeature() for clusters who were
     * invalidated because their center feature was used by a better cluster.
     * The neighbor features of this cluster have not necessarily been used
     * and might still be "active", i.e. their element mapping may still be accessed.
     * Therefore it should not contain references to a deleted cluster.
    */
    Size id = cluster.getId();
    for (const auto& element : cluster.getElements())
    {
      element_mapping[element.feature].erase(id);
    }
  }

void QTClusterFinder::createConsensusFeature_(ConsensusFeature& feature, 
                                              const double quality, 
                                              const QTCluster::Elements& elements)
  {
    feature.setQuality(quality);

    // the features of the current best cluster are inserted into the new consesus feature
    for (const auto& element : elements)
    {
      // Store the id of already used features (important: needs to be done
      // before updateClustering()) (not to be confused with the cluster id)
      already_used_.insert(element.feature);

      BaseFeature& elem_feat = const_cast<BaseFeature&>(element.feature->getFeature());
      feature.insert(element.map_index, elem_feat);
      if (elem_feat.metaValueExists("dc_charge_adducts"))
      {
        feature.setMetaValue(String(elem_feat.getUniqueId()), elem_feat.getMetaValue("dc_charge_adducts"));
      }
    }

    feature.computeConsensus();
  }

  void QTClusterFinder::updateClustering_(ElementMapping& element_mapping,
                                          const Grid& grid, 
                                          const QTCluster::Elements& elements,
                                          Heap& cluster_heads,
                                          const vector<Heap::handle_type>& handles,
                                          Size best_id)
  {
    for (const auto& element : elements)
    {
      const GridFeature* const curr_feature = element.feature;

      // ids of clusters the current feature belonged to
      unordered_set<Size>& cluster_ids = element_mapping[curr_feature];

      // delete the id of the current best cluster
      // we do not want to unnecessarily update it in the loop below
      cluster_ids.erase(best_id);

      // Identify all features that could potentially have been touched by this
      // Get all clusters that may potentially need updating

      ElementMapping tmp_element_mapping; // modify copy, then update

      for (const Size curr_id : cluster_ids)
      {
        // if we deleted the current id already, we shouldn't work with it
        // otherwise a segfault will happen
        // this is an ugly/quick fix and a better solution should be found in the future
        if (deleted_ids_.find(curr_id) != deleted_ids_.end()) continue;

        QTCluster& cluster = *handles[curr_id]; 

        // we do not want to update invalid features
        // (saves time and does not recompute the quality)
        if (!cluster.isInvalid())
        {
          // remove the elements of the new feature from the cluster

          if (cluster.update(elements))
          {
            // If update returns true, it means that at least one element was
            // removed from the cluster and we need to update that cluster

            ////////////////////////////////////////
            // Step 1: Iterate through all neighboring grid features and try to
            // add elements to the current cluster to replace the ones we just
            // removed
            addClusterElements_(grid, cluster);

            // update the heap, because the quality has changed
            cluster_heads.update_lazy(handles[curr_id]);

            ////////////////////////////////////////
            // Step 2: update element_mapping as the best feature for each
            // cluster may have changed

            for (const auto& neighbor : cluster.getAllNeighbors())
            {
              tmp_element_mapping[neighbor.feature].insert(curr_id);
            }
          }
        }
      }

      // we merge the tmp_element_mapping into the element_mapping after all clusters
      // that contained one feature of the current best cluster have been updated,
      // i.e. after every iteration of the outer loop
      for (const auto& cluster_ids: tmp_element_mapping)
      {
        for (const Size id : cluster_ids.second)
        {
          element_mapping[cluster_ids.first].insert(id);
        }
      }
    }

    // remove the current best from the heap and remember its id
    deleted_ids_.insert(best_id);
    cluster_heads.pop();
  }

  void QTClusterFinder::addClusterElements_(const Grid& grid, QTCluster& cluster)
  {
    cluster.initializeCluster();

#ifdef DEBUG_QTCLUSTERFINDER
    std::cout << " Compute Clustering: "<< x << " " << y << " with id " << center_feature->getFeature().getUniqueId() << std::endl;
    std::set<AASequence> a = cluster.getAnnotations();
    std::cout << " with annotations: ";
    for (std::set<AASequence>::iterator it = a.begin(); it != a.end(); ++it) std::cout << " " << *it;
    std::cout << std::endl;
#endif

    const int x = cluster.getXCoord(); 
    const int y = cluster.getYCoord(); 
    const GridFeature* center_feature = cluster.getCenterPoint();

    // iterate over neighboring grid cells (1st dimension):
    for (int i = x - 1; i <= x + 1; ++i)
    {
      // iterate over neighboring grid cells (2nd dimension):
      for (int j = y - 1; j <= y + 1; ++j)
      {
        auto act_pos = grid.grid_find(Grid::CellIndex(i, j));

        if (act_pos != grid.grid_end())
        {
          for (Grid::const_cell_iterator it_cell = act_pos->second.begin();
               it_cell != act_pos->second.end(); ++it_cell)
          {
            OpenMS::GridFeature* neighbor_feature = it_cell->second;

#ifdef DEBUG_QTCLUSTERFINDER
            std::cout << " considering to add feature " << neighbor_feature->getFeature().getUniqueId() << " to cluster " <<  center_feature->getFeature().getUniqueId()<< std::endl;
#endif

            // Skip features that we have already used -> we cannot add them to
            // be neighbors any more
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
              cluster.add(neighbor_feature, dist);
            }
          }
        }
      }
    }

    cluster.finalizeCluster();

#ifdef DEBUG_QTCLUSTERFINDER
    QTCluster::Elements elements = cluster.getElements();
    std::cout << " Done with cluster -> get quality " << cluster.getQuality() << " and nr elements " << elements.size() << std::endl;
    for (OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*>::const_iterator
         it = elements.begin(); it != elements.end(); ++it)
    {
      std::cout << "   = element id : " << it->second->getFeature().getUniqueId() << std::endl;
    }

    {
      std::set<AASequence> a = cluster.getAnnotations();
      std::cout << " FINAL with annotations: ";
      for (std::set<AASequence>::iterator it = a.begin(); it != a.end(); ++it) std::cout << " " << *it;
      std::cout << std::endl;
    }
#endif


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

  void QTClusterFinder::computeClustering_(const Grid& grid,
                                           Heap& cluster_heads,
                                           vector<QTCluster::BulkData>& cluster_data,
                                           vector<Heap::handle_type>& handles,
                                           ElementMapping& element_mapping)
  {
    cluster_heads.clear();
    already_used_.clear();
    cluster_data.clear();
    handles.clear();

    // do not remove this (will lead to segfault)
    // we need the pointers to cluster_data to stay valid,
    // therefore no reallocation is allowed to happen
    // we also reserve handles, because we don't know if we are allowed to move the handles
    // (the documentation of boost::heap does not tell us a lot about the handles)
    cluster_data.reserve(grid.size());
    handles.reserve(grid.size());

    Size id = 0;

    // FeatureDistance produces normalized distances (between 0 and 1):
    const double max_distance = 1.0;

    // iterate over all grid cells:
    for (Grid::const_iterator it = grid.begin(); it != grid.end(); ++it)
    {
      const Grid::CellIndex& act_coords = it.index();
      const Int x = act_coords[0], y = act_coords[1];

      const OpenMS::GridFeature* const center_feature = it->second;

      // construct empty data body for the new cluster and create the head afterwards
      cluster_data.emplace_back(center_feature, num_maps_, 
                                max_distance, x, y, id);
      
      QTCluster cluster(&cluster_data.back(), use_IDs_);

      addClusterElements_(grid, cluster);

      // push the cluster head of the new cluster into the heap
      // and the returned handle into our handle vector
      handles.push_back(cluster_heads.push(cluster));

      // register the new cluster for all its elements in the element mapping
      for (const auto& element : (*handles.back()).getElements())
      {
        element_mapping[element.feature].insert(id);
      }

      // next cluster gets the next id
      ++id;
    }
  }

  double QTClusterFinder::getDistance_(const OpenMS::GridFeature* left,
                                       const OpenMS::GridFeature* right)
  {
    return feature_distance_(left->getFeature(), right->getFeature()).second;
  }
  
  QTClusterFinder::~QTClusterFinder() = default;
} // namespace OpenMS
