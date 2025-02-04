// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>
#include <OpenMS/ANALYSIS/ID/IonIdentityMolecularNetworking.h>
#include <OpenMS/DATASTRUCTURES/Adduct.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace std;

namespace OpenMS
{

  FeatureGroupingAlgorithmKD::FeatureGroupingAlgorithmKD() :
    ProgressLogger(),
    feature_distance_(FeatureDistance())
  {
    setName("FeatureGroupingAlgorithmKD");

    defaults_.setValue("warp:enabled", "true", "Whether or not to internally warp feature RTs using LOWESS transformation before linking (reported RTs in results will always be the original RTs)");
    defaults_.setValidStrings("warp:enabled", {"true","false"});
    defaults_.setValue("warp:rt_tol", 100.0, "Width of RT tolerance window (sec)");
    defaults_.setMinFloat("warp:rt_tol", 0.0);
    defaults_.setValue("warp:mz_tol", 5.0, "m/z tolerance (in ppm or Da)");
    defaults_.setMinFloat("warp:mz_tol", 0.0);
    defaults_.setValue("warp:max_pairwise_log_fc", 0.5, "Maximum absolute log10 fold change between two compatible signals during compatibility graph construction. Two signals from different maps will not be connected by an edge in the compatibility graph if absolute log fold change exceeds this limit (they might still end up in the same connected component, however). Note: this does not limit fold changes in the linking stage, only during RT alignment, where we try to find high-quality alignment anchor points. Setting this to a value < 0 disables the FC check.", {"advanced"});
    defaults_.setValue("warp:min_rel_cc_size", 0.5, "Only connected components containing compatible features from at least max(2, (warp_min_occur * number_of_input_maps)) input maps are considered for computing the warping function", {"advanced"});
    defaults_.setMinFloat("warp:min_rel_cc_size", 0.0);
    defaults_.setMaxFloat("warp:min_rel_cc_size", 1.0);
    defaults_.setValue("warp:max_nr_conflicts", 0, "Allow up to this many conflicts (features from the same map) per connected component to be used for alignment (-1 means allow any number of conflicts)", {"advanced"});
    defaults_.setMinInt("warp:max_nr_conflicts", -1);

    defaults_.setValue("link:rt_tol", 30.0, "Width of RT tolerance window (sec)");
    defaults_.setMinFloat("link:rt_tol", 0.0);
    defaults_.setValue("link:mz_tol", 10.0, "m/z tolerance (in ppm or Da)");
    defaults_.setMinFloat("link:mz_tol", 0.0);
    defaults_.setValue("link:charge_merging","With_charge_zero","whether to disallow charge mismatches (Identical), allow to link charge zero (i.e., unknown charge state) with every charge state, or disregard charges (Any).");
    defaults_.setValidStrings("link:charge_merging", {"Identical", "With_charge_zero", "Any"});
    defaults_.setValue("link:adduct_merging","Any","whether to only allow the same adduct for linking (Identical), also allow linking features with adduct-free ones, or disregard adducts (Any).");
    defaults_.setValidStrings("link:adduct_merging", {"Identical", "With_unknown_adducts", "Any"});

    defaults_.setValue("mz_unit", "ppm", "Unit of m/z tolerance");
    defaults_.setValidStrings("mz_unit", {"ppm","Da"});
    defaults_.setValue("nr_partitions", 100, "Number of partitions in m/z space");
    defaults_.setMinInt("nr_partitions", 1);

    // FeatureDistance defaults
    defaults_.insert("", feature_distance_.getDefaults());

    // override some of them
    defaults_.setValue("distance_intensity:weight", 1.0);
    defaults_.setValue("distance_intensity:log_transform", "enabled");
    defaults_.addTag("distance_intensity:weight", "advanced");
    defaults_.addTag("distance_intensity:log_transform", "advanced");
    defaults_.remove("distance_RT:max_difference");
    defaults_.remove("distance_MZ:max_difference");
    defaults_.remove("distance_MZ:unit");
    defaults_.remove("ignore_charge");
    defaults_.remove("ignore_adduct");      

    // LOWESS defaults
    Param lowess_defaults;
    TransformationModelLowess::getDefaultParameters(lowess_defaults);
    for (Param::ParamIterator it = lowess_defaults.begin(); it != lowess_defaults.end(); ++it)
    {
      const_cast<Param::ParamEntry&>(*it).tags.insert("advanced");
    }
    defaults_.insert("LOWESS:", lowess_defaults);
    defaults_.setSectionDescription("LOWESS", "LOWESS parameters for internal RT transformations (only relevant if 'warp:enabled' is set to 'true')");

    defaultsToParam_();
    setLogType(CMD);
  }

  FeatureGroupingAlgorithmKD::~FeatureGroupingAlgorithmKD() = default;

  template <typename MapType>
  void FeatureGroupingAlgorithmKD::group_(const vector<MapType>& input_maps,
                                          ConsensusMap& out)
  {
    // set parameters
    String mz_unit(param_.getValue("mz_unit").toString());
    mz_ppm_ = mz_unit == "ppm";
    mz_tol_ = (double)(param_.getValue("link:mz_tol"));
    rt_tol_secs_ = (double)(param_.getValue("link:rt_tol"));

    // check that the number of maps is ok:
    if (input_maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "At least two maps must be given!");
    }

    out.clear(false);

    // collect all m/z values for partitioning, find intensity maximum
    vector<double> massrange;
    double max_intensity(0.0);
    for (typename vector<MapType>::const_iterator map_it = input_maps.begin();
         map_it != input_maps.end(); ++map_it)
    {
      for (typename MapType::const_iterator feat_it = map_it->begin();
          feat_it != map_it->end(); feat_it++)
      {
        massrange.push_back(feat_it->getMZ());
        double inty = feat_it->getIntensity();
        if (inty > max_intensity)
        {
          max_intensity = inty;
        }
      }
    }

    // set up distance functor
    Param distance_params;
    distance_params.insert("", param_.copy("distance_RT:"));
    distance_params.insert("", param_.copy("distance_MZ:"));
    distance_params.insert("", param_.copy("distance_intensity:"));
    distance_params.setValue("distance_RT:max_difference", rt_tol_secs_);
    distance_params.setValue("distance_MZ:max_difference", mz_tol_);
    distance_params.setValue("distance_MZ:unit", (mz_ppm_ ? "ppm" : "Da"));
    feature_distance_ = FeatureDistance(max_intensity, false);
    feature_distance_.setParameters(distance_params);

    // partition at boundaries -> this should be safe because there cannot be
    // any cluster reaching across boundaries

    sort(massrange.begin(), massrange.end());
    int pts_per_partition = massrange.size() / (int)(param_.getValue("nr_partitions"));

    double warp_mz_tol = (double)(param_.getValue("warp:mz_tol"));
    double max_mz_tol = max(mz_tol_, warp_mz_tol);

    // compute partition boundaries
    vector<double> partition_boundaries;
    partition_boundaries.push_back(massrange.front());
    for (size_t j = 0; j < massrange.size()-1; j++)
    {
      // minimal differences between two m/z values
      double massrange_diff = mz_ppm_ ? max_mz_tol * 1e-6 * massrange[j+1] : max_mz_tol;

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

    // ------------ compute RT transformation models ------------

    MapAlignmentAlgorithmKD aligner(input_maps.size(), param_);
    bool align = param_.getValue("warp:enabled").toString() == "true";
    if (align)
    {
      Size progress = 0;
      startProgress(0, partition_boundaries.size(), "computing RT transformations");
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

        // set up kd-tree
        KDTreeFeatureMaps kd_data(tmp_input_maps, param_);
        aligner.addRTFitData(kd_data);
        setProgress(progress++);
      }

      // fit LOWESS on RT fit data collected across all partitions
      try
      {
        aligner.fitLOWESS();
      }
      catch (Exception::BaseException& e)
      {
        OPENMS_LOG_ERROR << "Error: " << e.what() << endl;
        return;
      }

      endProgress();
    }

    // ------------ run alignment + feature linking on individual partitions ------------
    Size progress = 0;
    startProgress(0, partition_boundaries.size(), "linking features");
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

      // set up kd-tree
      KDTreeFeatureMaps kd_data(tmp_input_maps, param_);

      // alignment
      if (align)
      {
        aligner.transform(kd_data);
      }

      // link features
      runClustering_(kd_data, out);
      setProgress(progress++);
    }
    endProgress();
    
    postprocess_(input_maps, out);
  }

  void FeatureGroupingAlgorithmKD::group(const std::vector<FeatureMap>& maps,
                                         ConsensusMap& out)
  {
    group_(maps, out);
  }

  void FeatureGroupingAlgorithmKD::group(const std::vector<ConsensusMap>& maps,
                                         ConsensusMap& out)
  {
    group_(maps, out);
  }

  void FeatureGroupingAlgorithmKD::runClustering_(const KDTreeFeatureMaps& kd_data, ConsensusMap& out)
  {
    Size n = kd_data.size();

    // pass 1: initialize best potential clusters for all possible cluster centers
    set<Size> update_these;
    for (Size i = 0; i < kd_data.size(); ++i)
    {
      update_these.insert(i);
    }
    set<ClusterProxyKD> potential_clusters;
    vector<ClusterProxyKD> cluster_for_idx(n);
    vector<Int> assigned(n, false);
    updateClusterProxies_(potential_clusters, cluster_for_idx, update_these, assigned, kd_data);

    // pass 2: construct consensus features until all points assigned.
    while (!potential_clusters.empty())
    {
      // get index of current best cluster center (as defined by ClusterProxyKD::operator<)
      Size i = potential_clusters.begin()->getCenterIndex();

      // compile the actual list of sub feature indices for cluster with center i
      vector<Size> cf_indices;
      computeBestClusterForCenter_(i, cf_indices, assigned, kd_data);

      // add consensus feature
      addConsensusFeature_(cf_indices, kd_data, out);

      // mark selected sub features assigned and delete them from potential_clusters
      for (vector<Size>::const_iterator f_it = cf_indices.begin(); f_it != cf_indices.end(); ++f_it)
      {
        assigned[*f_it] = true;
        potential_clusters.erase(cluster_for_idx[*f_it]);
      }

      // compile set of all points whose neighborhoods will need updating
      update_these = set<Size>();
      for (vector<Size>::const_iterator f_it = cf_indices.begin(); f_it != cf_indices.end(); ++f_it)
      {
        vector<Size> f_neighbors;
        kd_data.getNeighborhood(*f_it, f_neighbors, rt_tol_secs_, mz_tol_, mz_ppm_, true);
        for (vector<Size>::const_iterator it = f_neighbors.begin(); it != f_neighbors.end(); ++it)
        {
          if (!assigned[*it])
          {
            update_these.insert(*it);
          }
        }
      }

      // now that the points are marked assigned, update the neighborhoods of their neighbors
      updateClusterProxies_(potential_clusters, cluster_for_idx, update_these, assigned, kd_data);
    }
  }



  void FeatureGroupingAlgorithmKD::updateClusterProxies_(set<ClusterProxyKD>& potential_clusters,
                                                         vector<ClusterProxyKD>& cluster_for_idx,
                                                         const set<Size>& update_these,
                                                         const vector<Int>& assigned,
                                                         const KDTreeFeatureMaps& kd_data)
  {
    for (set<Size>::const_iterator it = update_these.begin(); it != update_these.end(); ++it)
    {
      Size i = *it;
      const ClusterProxyKD& old_proxy = cluster_for_idx[i];
      vector<Size> unused;
      ClusterProxyKD new_proxy = computeBestClusterForCenter_(i, unused, assigned, kd_data);

      // only need to update if size and/or average distance have changed
      if (new_proxy != old_proxy)
      {
        potential_clusters.erase(old_proxy);
        cluster_for_idx[i] = new_proxy;
        potential_clusters.insert(new_proxy);
      }
    }
  }

  ClusterProxyKD FeatureGroupingAlgorithmKD::computeBestClusterForCenter_(Size i, vector<Size>& cf_indices, const vector<Int>& assigned, const KDTreeFeatureMaps& kd_data) const
  {
    //Parameters how to use charge/adduct information
    String merge_charge(param_.getValue("link:charge_merging").toString());
    String merge_adduct(param_.getValue("link:adduct_merging").toString());

    // compute i's neighborhood, together with a look-up table
    // map index -> corresponding points
    map<Size, vector<Size> > points_for_map_index;
    vector<Size> neighbors;
    kd_data.getNeighborhood(i, neighbors, rt_tol_secs_, mz_tol_, mz_ppm_, true);
    Int charge_i = kd_data.charge(i);
    const BaseFeature* f_i = kd_data.feature(i);
    for (vector<Size>::const_iterator it = neighbors.begin(); it != neighbors.end(); ++it)
    {
      // If the feature was already assigned, don't consider it at all!
      if (assigned[*it])
      {
        continue;
      }

      if (merge_charge == "Identical")
      {
        if (kd_data.charge(*it) != charge_i)
        {
          continue;
        }
      }
      // what to consider for linking with existing features _that have charge_. This ensures that we won't collect different non-zero charges.
      else if (merge_charge == "With_charge_zero")
      {
        if ((kd_data.charge(*it) != charge_i) && (kd_data.charge(*it) != 0))
        {
          continue;
        }
      }
      // else if (merge_charge == "Any")
      //{
      //  //we allow to merge all
      //}

      // analogous adduct block
      if (merge_adduct == "Identical")
      {
        // subcase 1: one has adduct, other not
        if (kd_data.feature(*it)->metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS) != f_i->metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
        {
          continue;
        }
        // subcase 2: both have adduct, but is it the same?
        if (kd_data.feature(*it)->metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
        {
          if (EmpiricalFormula(kd_data.feature(*it)->getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)) != EmpiricalFormula(f_i->getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)))
          {
            continue;
          }  
        }
      }
      // what to consider for linking with existing features _that have adduct_. If one has no adduct, it's fine
      // anyway. If one has an adduct we have to compare.
      else if (merge_adduct == "With_unknown_adducts")
      {
        // subcase1: *it has adduct, but i not. don't want to collect potentially different adducts to previous without adduct 
        if ((kd_data.feature(*it)->metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS)) && (!f_i->metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS)))
        {
          continue;
        }
        // subcase2: both have adduct
        if ((kd_data.feature(*it)->metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS)) && (f_i->metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS)))
        {
          // cheaper string check first, only check EF extensively if strings differ (might be just different element orders)
          if ((kd_data.feature(*it)->getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS) != f_i->getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)) &&
              (EmpiricalFormula(kd_data.feature(*it)->getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS)) != EmpiricalFormula(f_i->getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS))))
          {
            continue;
          }
        }
      }
      // else if (merge_adduct == "Any")
      //{
      //  //we allow to merge all
      //}

      // if everything is OK, add feature
      points_for_map_index[kd_data.mapIndex(*it)].push_back(*it);
    }
    // center i is always part of CF, no other points from i's map can be contained
    points_for_map_index[kd_data.mapIndex(i)] = vector<Size>(1, i);

    // compile list of sub feature indices and corresponding average distance
    double avg_distance = 0.0;
    for (map<Size, vector<Size> >::const_iterator it = points_for_map_index.begin();
         it != points_for_map_index.end();
         ++it)
    {
      const vector<Size>& candidates = it->second;

      // choose a point j with minimal distance to center i
      double min_dist = numeric_limits<double>::max();
      Size best_index = numeric_limits<Size>::max();
      for (vector<Size>::const_iterator c_it = candidates.begin(); c_it != candidates.end(); ++c_it)
      {
        double dist = const_cast<FeatureDistance&>(feature_distance_)(*(kd_data.feature(*c_it)),
                                                                      *(kd_data.feature(i))).second;

        if (dist < min_dist)
        {
          min_dist = dist;
          best_index = *c_it;
        }
      }
      cf_indices.push_back(best_index);
      avg_distance += min_dist;
    }
    avg_distance /= cf_indices.size();

    return ClusterProxyKD(cf_indices.size(), avg_distance, i);
  }

  void FeatureGroupingAlgorithmKD::addConsensusFeature_(const vector<Size>& indices, const KDTreeFeatureMaps& kd_data, ConsensusMap& out) const
  {
    ConsensusFeature cf;
    Adduct adduct;
    float avg_quality = 0;
    // determine best quality feature for adduct ion annotation (Constanst::UserParam::IIMN_BEST_ION)
    float best_quality = 0;
    size_t best_quality_index = 0;
    // collect the "Group" MetaValues of Features in a ConsensusFeature MetaValue (Constant::UserParam::IIMN_LINKED_GROUPS)
    vector<String> linked_groups;
    for (vector<Size>::const_iterator it = indices.begin(); it != indices.end(); ++it)
    {
      Size i = *it;
      cf.insert(kd_data.mapIndex(i), *(kd_data.feature(i)));
      avg_quality += kd_data.feature(i)->getQuality();
      if (kd_data.feature(i)->metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS) &&
         (kd_data.feature(i)->getQuality() > best_quality) &&
         (kd_data.feature(i)->getCharge()))
      {
       best_quality = kd_data.feature(i)->getQuality();
       best_quality_index = i;
      }
      if (kd_data.feature(i)->metaValueExists(Constants::UserParam::ADDUCT_GROUP))
      {
        linked_groups.emplace_back(kd_data.feature(i)->getMetaValue(Constants::UserParam::ADDUCT_GROUP));
      }
    }
    if (kd_data.feature(best_quality_index)->metaValueExists(Constants::UserParam::DC_CHARGE_ADDUCTS))
    {
      cf.setMetaValue(Constants::UserParam::IIMN_BEST_ION, 
                      adduct.toAdductString(kd_data.feature(best_quality_index)->getMetaValue(Constants::UserParam::DC_CHARGE_ADDUCTS),
                                            kd_data.feature(best_quality_index)->getCharge()));
    }
    if (!linked_groups.empty())
    {
      cf.setMetaValue(Constants::UserParam::IIMN_LINKED_GROUPS, linked_groups);
    }
    avg_quality /= indices.size();
    cf.setQuality(avg_quality);
    cf.computeConsensus();
    out.push_back(cf);
  }

} // namespace OpenMS
