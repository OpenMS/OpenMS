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
    vector<Handle> handles;

    // "cold" cluster bodies, where most of their data lies
    vector<QTCluster::Data_> cluster_data;

    // map to get ids from clusters, who contain a certain grid feature
    ElementMapping element_mapping;

    computeClustering_(grid, cluster_heads, cluster_data, handles);

    // number of clusters == number of data points:
    Size size = cluster_heads.size();

    // fill the element mapping
    for (Heap::iterator it = cluster_heads.begin();
         it != cluster_heads.end(); ++it)
    {
      NeighborMap const& neigh = it->getAllNeighborsDirect();
      for (NeighborMap::const_iterator n_it = neigh.begin(); n_it != neigh.end(); ++n_it)
      {
        GridFeature* gf_ptr = n_it->second.second;
        element_mapping[gf_ptr].insert(it->getId());
      }

      // ensure that the cluster center is in the list as well
      OpenMS::GridFeature* center_feature = it->getCenterPoint();
      element_mapping[center_feature].insert(it->getId());
    }

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

      if(made_feature)
      {
        result_map.push_back(consensus_feature);
      }
      if (do_progress) logger.setProgress(progress++);
    }

    if (do_progress) logger.endProgress();
  }

  bool QTClusterFinder::makeConsensusFeature_(Heap& cluster_heads,
                                              ConsensusFeature& feature,
                                              ElementMapping& element_mapping,
                                              Grid& grid,
                                              vector<Handle> const& handles)
  {

    // pop until the top is valid
    while (cluster_heads.top().isInvalid())
    {
      removeTopFromHeap_(cluster_heads, cluster_heads.top(), element_mapping);
      if (cluster_heads.empty()) return false;
    }

    QTCluster const& best = cluster_heads.top();
    NeighborMap const& elements = best.getElementsBeforeDestruction();

    #ifdef DEBUG_QTCLUSTERFINDER
    std::cout << "Elements: " << elements.size() << " with best "
         << best->getQuality() << " invalid " << best->isInvalid() << std::endl;
    #endif

    // create consensus feature from best cluster:
    createConsensusFeature_(feature, best.getCurrentQuality(), elements);

    #ifdef DEBUG_QTCLUSTERFINDER
    std::cout << " create new consensus feature " << feature.getRT() << " " << feature.getMZ() << " from " << best->getCenterPoint()->getFeature().getUniqueId() << std::endl;
    for (OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*>::const_iterator
         it = elements.begin(); it != elements.end(); ++it)
    {
      std::cout << "   = element id : " << it->second->getFeature().getUniqueId() << std::endl;
    }
    #endif

    // update the clustering:
    // 1. remove current "best" cluster from list
    // 2. update all clusters accordingly by removing already used elements
    // 3. Invalidate elements whose central has been used already

    updateClustering_(element_mapping, grid, elements, cluster_heads, handles, best.getId());

    return true;
  }

  void QTClusterFinder::removeTopFromHeap_(Heap& cluster_heads,
                                           QTCluster const& cluster,
                                           ElementMapping& element_mapping)
  {
    NeighborMap const& elements = cluster.getElementsBeforeDestruction();
    for(NeighborMap::const_iterator it = elements.begin(); it != elements.end(); ++it)
    {
      unordered_set<Size>& cluster_ids = element_mapping[it->second.second];
      cluster_ids.erase(cluster.getId());
    }

    cluster_heads.pop();
  }

  void QTClusterFinder::updateClustering_(ElementMapping &element_mapping,
                                          Grid const &grid, 
                                          NeighborMap const& elements,
                                          Heap& cluster_heads,
                                          vector<Handle> const& handles,
                                          Size best_id)
  {
    for (NeighborMap::const_iterator
        it = elements.begin(); it != elements.end(); ++it)
    {
      // delete id of the best cluster from element mapping
      unordered_set<Size>& cluster_ids = element_mapping[it->second.second];
      cluster_ids.erase(best_id);

      // Identify all features that could potentially have been touched by this
      // Get all clusters that may potentially need updating

      ElementMapping tmp_element_mapping; // modify copy, then update

      for (unordered_set<Size>::iterator
           id_it  = cluster_ids.begin();
           id_it != cluster_ids.end(); ++id_it)
      {
        QTCluster& cluster = *handles[*id_it]; 
        
        // we do not want to update invalid features (saves time and does not
        // recompute the quality)
        if (!cluster.isInvalid())
        {
          // remove the elements of the new feature from the cluster
          if (cluster.update(elements))
          {
            // If update returns true, it means that at least one element was
            // removed from the cluster and we need to update that cluster

            // Get the coordinates of the current cluster
            const Int x = cluster.getXCoord();
            const Int y = cluster.getYCoord();

            ////////////////////////////////////////
            // Step 1: Iterate through all neighboring grid features and try to
            // add elements to the current cluster to replace the ones we just
            // removed
            const GridFeature* center_feature = cluster.getCenterPoint();
            addClusterElements_(x, y, grid, cluster, center_feature);

            // update the heap, because the quality has changed
            cluster_heads.update_lazy(handles[*id_it]);

            ////////////////////////////////////////
            // Step 2: update element_mapping as the best feature for each
            // cluster may have changed

            NeighborMap neigh = cluster.getAllNeighborsDirect();
            for (NeighborMap::iterator n_it = neigh.begin(); n_it != neigh.end(); ++n_it)
            {
              GridFeature* gf_ptr = n_it->second.second;
              tmp_element_mapping[gf_ptr].insert(cluster.getId());
            }
          }
        }
      }

      for (ElementMapping::iterator it = tmp_element_mapping.begin();
          it != tmp_element_mapping.end(); ++it )
      {
        for (unordered_set<Size>::iterator it2 = it->second.begin();
            it2 != it->second.end(); ++it2)
        {
          element_mapping[ it->first ].insert(*it2);
        }
      }
    }

    // remove the current best from the heap
    cluster_heads.pop();
  }

  void QTClusterFinder::createConsensusFeature_(ConsensusFeature &feature, double quality, NeighborMap const& elements)
  {
    feature.setQuality(quality);
    for (NeighborMap::const_iterator
         it = elements.begin(); it != elements.end(); ++it)
    {
      // Store the id of already used features (important: needs to be done
      // before updateClustering)
      already_used_.insert(it->second.second);

      BaseFeature& elem_feat = const_cast<BaseFeature&>(it->second.second->getFeature());
      feature.insert(it->first, elem_feat);
      if (elem_feat.metaValueExists("dc_charge_adducts"))
      {
        feature.setMetaValue(String(elem_feat.getUniqueId()), elem_feat.getMetaValue("dc_charge_adducts"));
      }
    }

    feature.computeConsensus();
  }

  void QTClusterFinder::addClusterElements_(int x, int y, const Grid& grid, QTCluster& cluster,
    const OpenMS::GridFeature* center_feature)
  {
    cluster.initializeCluster();

    #ifdef DEBUG_QTCLUSTERFINDER
    std::cout << " Compute Clustering: "<< x << " " << y << " with id " << center_feature->getFeature().getUniqueId() << std::endl;
    std::set<AASequence> a = cluster.getAnnotations();
    std::cout << " with annotations: ";
    for (std::set<AASequence>::iterator it = a.begin(); it != a.end(); ++it) std::cout << " " << *it;
    std::cout << std::endl;
    #endif


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
    OpenMSBoost::unordered_map<Size, OpenMS::GridFeature*> elements;
    // this shouldn't be done anymore
    // cluster.getElements(elements);
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

  void QTClusterFinder::computeClustering_(Grid& grid,
                                           Heap& cluster_heads,
                                           vector<QTCluster::Data_>& cluster_data,
                                           vector<Handle>& handles)
  {
    cluster_heads.clear();
    already_used_.clear();
    cluster_data.clear();
    handles.clear();

    // do not remove this (will lead to segfault)
    cluster_data.reserve(grid.size());
    handles.reserve(grid.size());

    Size id = 0;

    // FeatureDistance produces normalized distances (between 0 and 1):
    const double max_distance = 1.0;

    // iterate over all grid cells:
    for (Grid::iterator it = grid.begin(); it != grid.end(); ++it)
    {
      const Grid::CellIndex& act_coords = it.index();
      const Int x = act_coords[0], y = act_coords[1];

      OpenMS::GridFeature* center_feature = it->second;

      // construct empty data body for the new cluster
      cluster_data.emplace_back();
      
      QTCluster cluster(&cluster_data.back(), center_feature, num_maps_, 
                                        max_distance, use_IDs_, x, y, id);

      addClusterElements_(x, y, grid, cluster, center_feature);

      // push cluster head into the heap
      handles[id] = cluster_heads.push(cluster);

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
