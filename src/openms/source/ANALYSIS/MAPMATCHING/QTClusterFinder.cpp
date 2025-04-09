// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Steffen Sass, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/QTClusterFinder.h>

#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/MATH/MathFunctions.h>

//#define DEBUG_QTCLUSTERFINDER_IDS

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
    setName("QTClusterFinder");

    defaults_.setValue("use_identifications", "false", "Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).");
    defaults_.setValidStrings("use_identifications", {"true","false"});
    defaults_.setValue("nr_partitions", 100, "How many partitions in m/z space should be used for the algorithm (more partitions means faster runtime and more memory efficient execution).");
    defaults_.setMinInt("nr_partitions", 1);
    defaults_.setValue("min_nr_diffs_per_bin", 50, "If IDs are used: How many differences from matching IDs should be used to calculate a linking tolerance for unIDed features in an RT region. RT regions will be extended until that number is reached.");
    defaults_.setMinInt("min_nr_diffs_per_bin", 5);
    defaults_.setValue("min_IDscore_forTolCalc", 1., "If IDs are used: What is the minimum score of an ID to assume a reliable match for tolerance calculation. Check your current score type!");
    defaults_.setValue("noID_penalty", 0.0, "If IDs are used: For the normalized distances, how high should the penalty for missing IDs be? 0 = no bias, 1 = IDs inside the max tolerances always preferred (even if much further away).");
    defaults_.setMinFloat("noID_penalty", 0.0);
    defaults_.setMaxFloat("noID_penalty", 1.0);

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
    use_IDs_ = param_.getValue("use_identifications").toBool();
    nr_partitions_ = param_.getValue("nr_partitions");
    min_nr_diffs_per_bin_ = param_.getValue("min_nr_diffs_per_bin");
    min_score_ = param_.getValue("min_IDscore_forTolCalc");
    noID_penalty_ = param_.getValue("noID_penalty");
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
    distance_params.remove("min_nr_diffs_per_bin");
    distance_params.remove("min_IDscore_forTolCalc");
    distance_params.remove("noID_penalty");
    feature_distance_ = FeatureDistance(max_intensity, true);
    feature_distance_.setParameters(distance_params);
  }

  template <typename MapType>
  void QTClusterFinder::run_(const vector<MapType>& input_maps,
                             ConsensusMap& result_map)
  {
    // update parameters (dummy)
    setParameters_(1, 1);

    if (use_IDs_)
    {
      // map string "modified sequence/charge" to all RTs the feature has been observed in the different maps
      std::unordered_map<String, std::vector<double>> ided_feat_rts;
      //std::unordered_map<String, std::vector<const typename MapType::FeatureType*>> ided_feats;
      double minRT = std::numeric_limits<double>::max();
      for (auto& map : input_maps)
      {
        for (auto feat : map) //OMS_CODING_TEST_EXCLUDE Note: needs copy to sort
        {
          if (feat.getRT() < minRT) minRT = feat.getRT();
          auto& pepIDs = feat.getPeptideIdentifications();
          if (!pepIDs.empty())
          {
            //TODO I think we sort in run_internal again. Could be avoided.
            feat.sortPeptideIdentifications();
            auto& hits = pepIDs[0].getHits();
            if (!hits.empty())
            {
              if ((hits[0].getScore() > min_score_ && pepIDs[0].isHigherScoreBetter()) ||
                  (hits[0].getScore() < min_score_ && !pepIDs[0].isHigherScoreBetter()))
              {
                //TODO we could loosen the score filtering by requiring only ONE IDed feature of a peptide to pass the threshold.
                // Would require a second pass though
                const String key = pepIDs[0].getHits()[0].getSequence().toString() + "/" + feat.getCharge();
                const auto [it, inserted] = ided_feat_rts.emplace(key, std::vector<double>{feat.getRT()});
                if (!inserted) // already present
                {
                  it->second.push_back(feat.getRT());
                }
                //TODO we could score the whole feature instead of just the RT to calculate tolerances based on
                // a combined score (RT/mz; using the scoring function of this class) instead of just RT
                /*const auto it_inserted_feat = ided_feats.emplace(key, std::vector<const typename MapType::FeatureType*>{&feat});
                if (!it_inserted_feat.second)
                {
                  it_inserted_feat.first->second.push_back(&feat);
                }*/
              }
            }
          }
        }
      }

      //Note: this does not differentiate between the variety of differences between distinct map pairs. E.g.
      // differences between map 1 and map 2 might be usually very small (e.g. they are replicates), while
      // differences between map 1 and map 3 are large, since they are different conditions. But we might lose
      // robust estimates and use more memory if we split them.
      std::vector<std::pair<double,std::vector<double>>> medians_diffs;
      medians_diffs.resize(ided_feat_rts.size());
      Size c = 0;
      // for every ID, calculate median RT and differences
      for (auto& id_rts : ided_feat_rts)
      {
        #ifdef DEBUG_QTCLUSTERFINDER_IDS
        std::cout << "Stats for " << id_rts.first << ": ";
        for (const auto& rt : id_rts.second)
        {
          std::cout << rt << ", ";
        }
        std::cout << std::endl;
        #endif

        auto& rts = id_rts.second;
        std::sort(rts.begin(),rts.end());
        medians_diffs[c].first = rts[rts.size()/2];
        medians_diffs[c].second.reserve(rts.size()-1);
        Size i = 0;
        for (const auto& rt : rts)
        {
          if (i++ == rts.size()/2) continue;
          //TODO would relative diffs solve the RT dependency issue sufficiently? Probably not with non-linear shifts
          //Note: One the one hand using abs. val. destroys the Gaussian distribution, on the other hand, if you only have
          // two RTs for an ID, you will always have a negative difference from the median (and therefore the distribution
          // would be biased anyway). We could use a stable median for evenly sized vectors but that would not reflect
          // real distances then. Or always add positive AND negative differences.
          medians_diffs[c].second.push_back(std::fabs(rt - medians_diffs[c].first));
        }
        c++;
      }

      //TODO check if we can assume sorted
      std::sort(medians_diffs.begin(), medians_diffs.end());

      Size cnt = 0;
      vector<double> tmp_diffs, last_tmp_diffs;
      // we calculate minRT (instead of starting first bin at 0) since RTs may start in the "negative" region after alignment.
      double start_rt = minRT;
      double min_tolerance = 20;
      double tol, q2, q3;
      OPENMS_LOG_INFO << "Calculating RT linking tolerance bins...\n";
      OPENMS_LOG_INFO << "RT_bin_start, Tolerance" << std::endl;

      // For every pair of median RT and differences, collect
      // differences until min_nr_diffs_per_bin_ is reached, then add the
      // start of the current bin and a tolerance based on Double MAD https://aakinshin.net/posts/harrell-davis-double-mad-outlier-detector/
      for (const auto& med_diffs : medians_diffs)
      {
        if (tmp_diffs.size() > min_nr_diffs_per_bin_ && !med_diffs.second.empty())
        {
          std::sort(tmp_diffs.begin(), tmp_diffs.end());
          // calculate allowed tolerance
          //q1 = quantile_(tmp_diffs, 0.25);
          q2 = Math::quantile(tmp_diffs, 0.5);
          q3 = Math::quantile(tmp_diffs, 0.75);
          //q95 = quantile_(tmp_diffs, 0.95);
          //iqr = q3 - q1;
          //tol = max(min_tolerance, max(fabs(q3 + iqr_mult * iqr), fabs(q1 - iqr_mult * iqr)));
          //tol = q95;
          //tol = 2. * 1.4826 * q2;
          tol = max(min_tolerance, q2 + 2. * 1.4826 * (q3-q2));
          bin_tolerances_.insert(make_pair(start_rt, tol));

          OPENMS_LOG_INFO << start_rt << ", " << tol << std::endl;
          #ifdef DEBUG_QTCLUSTERFINDER_IDS
          std::cout << "Differences used: ";
          for (const auto& diff : tmp_diffs)
          {
            std::cout << diff << ", ";
          }
          std::cout << std::endl;
          #endif
          std::swap(tmp_diffs, last_tmp_diffs);
          tmp_diffs.clear();
          if (cnt > 0) start_rt = (med_diffs.first + medians_diffs[cnt-1].first)/2;
        }
        else
        {
          tmp_diffs.insert(tmp_diffs.end(), med_diffs.second.begin(), med_diffs.second.end());
        }
        cnt++;
      }
      // calculate allowed tolerance
      std::sort(tmp_diffs.begin(), tmp_diffs.end());
      std::vector<double> last_and_before_diffs;
      last_and_before_diffs.reserve(tmp_diffs.size() + last_tmp_diffs.size());
      std::merge(tmp_diffs.begin(), tmp_diffs.end(), last_tmp_diffs.begin(), last_tmp_diffs.end(), std::back_inserter(last_and_before_diffs));
      if (!last_and_before_diffs.empty())
      {
        q2 = Math::quantile(last_and_before_diffs, 0.5);
        q3 = Math::quantile(last_and_before_diffs, 0.75);
        tol = max(min_tolerance, q2 + 2. * 1.4826 * (q3-q2));
        bin_tolerances_.insert(make_pair(start_rt, tol));

        OPENMS_LOG_INFO << start_rt << ", " << tol << std::endl;
        #ifdef DEBUG_QTCLUSTERFINDER_IDS
        std::cout << "Differences used: ";
        for (const auto& diff : last_and_before_diffs)
        {
          std::cout << diff << ", ";
        }
        std::cout << std::endl;
        #endif

        #ifdef DEBUG_QTCLUSTERFINDER_IDS
        std::cout << "size of last bin: " << last_and_before_diffs.size() << std::endl;
        #endif
      }
      last_and_before_diffs.clear();
      last_tmp_diffs.clear();
      tmp_diffs.clear();
    }

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
      int pts_per_partition = int(massrange.size()) / nr_partitions_;

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
      max_intensity = max(max_intensity, map_it->getMaxIntensity());
      max_mz = max(max_mz, map_it->getMaxMZ());
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
        auto& bfeat = const_cast<BaseFeature&>(gfeat.getFeature());
        for (auto& pep : bfeat.getPeptideIdentifications())
        {
          pep.sort();
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
      // pops heap until a valid best cluster or empty, makes a consensusFeature and updates
      // other clusters affected by the inclusion of this cluster
      bool made_feature = makeConsensusFeature_(cluster_heads, consensus_feature, 
                                                element_mapping, grid, handles);

      if (made_feature)
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
                                              const Grid& grid,
                                              const vector<Heap::handle_type>& handles)
  {
    // pop until the top is valid
    while (cluster_heads.top().isInvalid())
    {
      removeFromElementMapping_(cluster_heads.top(), element_mapping);
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

    Adduct adduct;
    // determine best quality feature for adduct ion annotation (Constanst::UserParam::IIMN_BEST_ION)
    float best_quality = 0;
    size_t best_quality_index = 0;
    // collect the "Group" MetaValues of Features in a ConsensusFeature MetaValue (Constanst::UserParam::IIMN_LINKED_GROUPS)
    vector<String> linked_groups;
    // the features of the current best cluster are inserted into the new consensus feature
    for (const auto& element : elements)
    {
      // Store the id of already used features (important: needs to be done
      // before updateClustering()) (not to be confused with the cluster id)
      already_used_.insert(element.feature);

      BaseFeature& elem_feat = const_cast<BaseFeature&>(element.feature->getFeature());
      feature.insert(element.map_index, elem_feat);
      if (elem_feat.metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
      {
        feature.setMetaValue(String(elem_feat.getUniqueId()), elem_feat.getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS));
      }
      if (elem_feat.metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS) && (elem_feat.getQuality() > best_quality))
      {
        feature.setMetaValue(Constants::UserParam::IIMN_BEST_ION, elem_feat.getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS));
        best_quality = elem_feat.getQuality();
      }
      if (elem_feat.metaValueExists(Constants::UserParam::ADDUCT_GROUP))
      {
        linked_groups.emplace_back(elem_feat.getMetaValue(Constants::UserParam::ADDUCT_GROUP));
      }
    }
    if (elements[best_quality_index].feature->getFeature().metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
    {
      feature.setMetaValue(Constants::UserParam::IIMN_BEST_ION, 
                      adduct.toAdductString(elements[best_quality_index].feature->getFeature().getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS),
                                            elements[best_quality_index].feature->getFeature().getCharge()));
    }
    if (!linked_groups.empty())
    {
      feature.setMetaValue(Constants::UserParam::IIMN_LINKED_GROUPS, linked_groups);
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
    // remove the current best from the heap and consolidate the heap from previous lazy updates
    // we cannot pop at the end since update_lazy may theoretically change top_element immediately.
    cluster_heads.pop();

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

            /*
            ////////////////////////////////////////
            Step 1: Iterate through all neighboring grid features and try to
            add elements to the current cluster to replace the ones we just
            removed

            Before that we must delete this clusters id from the element mapping. (important!)
            It is possible that addClusterElements_() removes features from the cluster 
            we are updating. (Through finalizeCluster_ -> computeQuality_ -> optimizeAnnotations).
            These are not to be confused with the features we removed
            because they are part of the current best cluster. Those are removed in 
            QTCluster::update (above).

            If this happens, the element mapping for the additionally removed features 
            (which are valid and unused!) still contains the id of the cluster which 
            we are currently updating. But the cluster does not contain the feature anymore. 
            When the cluster is deleted, the element mapping for the removed feature doesn't 
            get updated. The element mapping for the feature then contains an id of a 
            deleted cluster, which will surely lead to a segfault when the feature is actually 
            used in another cluster later.

            TODO Check guarantee that addClusterElements does not add a feature that was removed
             earlier in the loop. Should not happen because they are in the already_used set by now.
            */
            removeFromElementMapping_(cluster, element_mapping);

            // re-add closest cluster elements that were not used yet.
            addClusterElements_(grid, cluster);

            // update the heap, because the quality has changed
            // compares with top_element to see if a different node needs to be popped now.
            // for comparison getQuality() is called for the clusters here
            // TODO check if we can guarantee cluster_heads.increase/decrease since they may have
            //  better theoretical runtimes although a lazy update until the next pop is probably not bad
            cluster_heads.update_lazy(handles[curr_id]);

            ////////////////////////////////////////
            // Step 2: reinsert the updated cluster's features into a temporary element mapping.
            // This can be merged later since the methods called in the loop here seem not to access the mapping.
            for (const auto& neighbor : cluster.getElements())
            {
              tmp_element_mapping[neighbor.feature].insert(curr_id);
            }
          }
        }
      }

      // we merge the tmp_element_mapping into the element_mapping after all clusters
      // that contained one feature of the current best cluster have been updated,
      // i.e. after every iteration of the outer loop
      for (const auto& feat_clusterids : tmp_element_mapping)
      {
        for (const Size id : feat_clusterids.second)
        {
          element_mapping[feat_clusterids.first].insert(id);
        }
      }
    }
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
              // if IDs are used during linking, check if unidentified features are too far off from "usual" RT shifts
              // in that region
              if (use_IDs_ && neighbor_feature->getAnnotations().empty())
              {
                double rt_dist = std::fabs(neighbor_feature->getRT() - center_feature->getRT());
                if (distIsOutlier_(rt_dist, center_feature->getRT())) continue;
                dist += noID_penalty_;
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
      std::set<AASequence> ax = cluster.getAnnotations();
      std::cout << " FINAL with annotations: ";
      for (std::set<AASequence>::iterator it = ax.begin(); it != ax.end(); ++it) std::cout << " " << *it;
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

    // FeatureDistance produces normalized distances (between 0 and 1 plus a possible noID penalty):
    const double max_distance = 1.0 + noID_penalty_;

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

  bool QTClusterFinder::distIsOutlier_(double dist, double rt)
  {
    if (bin_tolerances_.empty()) return false;
    auto it = bin_tolerances_.upper_bound(rt);
    if (it == bin_tolerances_.begin()) return dist >= it->second;
    return dist >= (--it)->second;
  }
  
  QTClusterFinder::~QTClusterFinder() = default;
} // namespace OpenMS
