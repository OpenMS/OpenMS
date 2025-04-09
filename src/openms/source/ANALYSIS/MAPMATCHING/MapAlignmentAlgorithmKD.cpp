// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <queue>

using namespace std;

namespace OpenMS
{

MapAlignmentAlgorithmKD::MapAlignmentAlgorithmKD(Size num_maps, const Param& param) :
  fit_data_(num_maps),
  transformations_(num_maps),
  param_(param),
  max_pairwise_log_fc_(-1)
{
  updateMembers_();
}

MapAlignmentAlgorithmKD::~MapAlignmentAlgorithmKD()
{
  for (vector<TransformationModelLowess*>::iterator it = transformations_.begin();
       it != transformations_.end(); ++it)
  {
    delete *it;
  }
}

void MapAlignmentAlgorithmKD::addRTFitData(const KDTreeFeatureMaps& kd_data)
{
  // compute connected components
  map<Size, vector<Size> > ccs;
  getCCs_(kd_data, ccs);

  // keep only conflict-free CCs of sufficient size
  map<Size, vector<Size> > filtered_ccs;
  filterCCs_(kd_data, ccs, filtered_ccs);
  // save some memory
  ccs.clear();

  // compute average RTs for all CCs
  map<Size, double> avg_rts;
  for (map<Size, vector<Size> >::const_iterator it = filtered_ccs.begin(); it != filtered_ccs.end(); ++it)
  {
    double avg_rt = 0;
    Size cc_index = it->first;
    const vector<Size>& cc = it->second;
    for (vector<Size>::const_iterator cc_it = cc.begin(); cc_it != cc.end(); ++cc_it)
    {
      Size i = *cc_it;
      avg_rt += kd_data.rt(i);
    }
    avg_rt /= cc.size();
    avg_rts[cc_index] = avg_rt;
  }

  // generate fit data for each map, add to fit_data_
  for (map<Size, vector<Size> >::const_iterator it = filtered_ccs.begin(); it != filtered_ccs.end(); ++it)
  {
    Size cc_index = it->first;
    const vector<Size>& cc = it->second;
    for (vector<Size>::const_iterator cc_it = cc.begin(); cc_it != cc.end(); ++cc_it)
    {
      Size i = *cc_it;
      double rt = kd_data.rt(i);
      double avg_rt = avg_rts[cc_index];
      fit_data_[kd_data.mapIndex(i)].push_back(make_pair(rt, avg_rt));
    }
  }
}

void MapAlignmentAlgorithmKD::fitLOWESS()
{
  Size num_maps = fit_data_.size();
  for (Size i = 0; i < num_maps; ++i)
  {
    Size n = fit_data_[i].size();
    const Param& lowess_param = param_.copy("LOWESS:", true);
    if (n < 50)
    {
      OPENMS_LOG_WARN << "Warning: Only " << n << " data points for LOWESS fit of map " << i << ". Consider adjusting RT or m/z tolerance or max_pairwise_log_fc, decreasing min_rel_cc_size, or increasing max_nr_conflicts." << endl;
      TransformationModel::DataPoints identity = {{0,0}, {1,1}, {1e6,1e6}};
      transformations_[i] = new TransformationModelLowess(identity, lowess_param);
    }
    else
    {
      transformations_[i] = new TransformationModelLowess(fit_data_[i], lowess_param);
    }
  }
}

void MapAlignmentAlgorithmKD::transform(KDTreeFeatureMaps& kd_data) const
{
  // apply transformations to kd_data
  kd_data.applyTransformations(transformations_);

  // re-optimize kd-tree
  kd_data.optimizeTree();
}

Size MapAlignmentAlgorithmKD::computeCCs_(const KDTreeFeatureMaps& kd_data, vector<Size>& result) const
{
  //compute CCs by means of repeated BFS (without actually storing the graph (edges) in memory)

  Size num_nodes = kd_data.size();

  //clear CC indices
  result.clear();
  result.resize(num_nodes, numeric_limits<Size>::max());

  //set up data structures
  queue<Size> bfs_queue;
  vector<Int> bfs_visited(num_nodes, false);
  Size search_pos = 0;
  Size cc_index = 0;

  //BFS until every node has been visited
  while (true)
  {
    bool finished = true;
    for (Size i = search_pos; i < num_nodes; ++i)
    {
      if (!bfs_visited[i])
      {
        bfs_queue.push(i);
        bfs_visited[i] = true;
        finished = false;
        search_pos = i + 1;
        break;
      }
    }
    if (finished) break;

    while (!bfs_queue.empty())
    {
      Size i = bfs_queue.front();
      bfs_queue.pop();
      result[i] = cc_index;

      vector<Size> compatible_features;
      kd_data.getNeighborhood(i, compatible_features, rt_tol_secs_, mz_tol_, mz_ppm_, false, max_pairwise_log_fc_);
      for (vector<Size>::const_iterator it = compatible_features.begin();
           it != compatible_features.end();
           ++it)
      {
        Size j = *it;
        if (!bfs_visited[j])
        {
          bfs_queue.push(j);
          bfs_visited[j] = true;
        }
      }
    }
    ++cc_index;
  }

  return cc_index;
}

void MapAlignmentAlgorithmKD::getCCs_(const KDTreeFeatureMaps& kd_data, map<Size, vector<Size> >& result) const
{
  vector<Size> cc_index;
  computeCCs_(kd_data, cc_index);

  result.clear();
  for (Size i = 0; i < kd_data.size(); ++i)
  {
    result[cc_index[i]].push_back(i);
  }
}

void MapAlignmentAlgorithmKD::filterCCs_(const KDTreeFeatureMaps& kd_data, const map<Size, vector<Size> >& ccs, map<Size, vector<Size> >& filtered_ccs) const
{
  Size num_maps = fit_data_.size();
  Size min_size = max(2.0, (double)(param_.getValue("warp:min_rel_cc_size")) * (double)num_maps);
  int max_nr_conflicts = (int)param_.getValue("warp:max_nr_conflicts");
  filtered_ccs.clear();

  for (map<Size, vector<Size> >::const_iterator it = ccs.begin(); it != ccs.end(); ++it)
  {
    const vector<Size>& cc = it->second;

    // size OK?
    if (cc.size() < min_size)
    {
      // nope
      continue;
    }

    // charges compatible?
    set<int> charges;
    for (vector<Size>::const_iterator idx_it = cc.begin(); idx_it != cc.end(); ++idx_it)
    {
      int z = kd_data.charge(*idx_it);
      if (z != 0)
      {
        charges.insert(z);
      }
      if (charges.size() > 1)
      {
        // nope
        continue;
      }
    }

    // check for conflicts
    bool passes = true;
    if (max_nr_conflicts != -1)
    {
      set<Size> map_indices;
      int nr_conflicts = 0;
      for (vector<Size>::const_iterator idx_it = cc.begin(); idx_it != cc.end(); ++idx_it)
      {
        // filter out if too many features from same map
        Size map_idx = kd_data.mapIndex(*idx_it);
        if (map_indices.find(map_idx) != map_indices.end())
        {
          if (++nr_conflicts > max_nr_conflicts)
          {
            passes = false;
            break;
          }
        }
        else
        {
          map_indices.insert(map_idx);
        }
      }
    }

    if (passes)
    {
      filtered_ccs[it->first] = cc;
    }
  }
}

void MapAlignmentAlgorithmKD::updateMembers_()
{
  if (param_.empty()) return;

  rt_tol_secs_ = (double)(param_.getValue("warp:rt_tol"));
  mz_tol_ = (double)(param_.getValue("warp:mz_tol"));
  mz_ppm_ = (param_.getValue("mz_unit").toString() == "ppm");
  max_pairwise_log_fc_ = param_.getValue("warp:max_pairwise_log_fc");
}

}
