// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Michal Startek $
// $Authors: Michal Startek $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmWNet.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <wnetalign/aligner.hpp>
#include <wnetalign/spectrum.hpp>

#include <numeric>
#include <set>

using namespace std;

namespace OpenMS
{

  FeatureGroupingAlgorithmWNet::FeatureGroupingAlgorithmWNet() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmWNet");

    defaults_.setValue("distance_metric", "LINF", "Distance metric for comparing feature positions (after m/z scaling)");
    defaults_.setValidStrings("distance_metric", {"L1", "L2", "LINF"});

    defaults_.setValue("max_distance", 100.0,
      "Maximum distance between features to consider a match. "
      "After auto m/z scaling, both dimensions are in RT-like units (seconds), "
      "so this threshold is effectively in seconds. "
      "Comparable to distance_RT:max_difference in other linkers.");
    defaults_.setMinFloat("max_distance", 0.0);

    defaults_.setValue("trash_cost", 100.0,
      "Cost of leaving a feature unmatched. When equal to max_distance, "
      "a match at the maximum allowed distance is as expensive as no match. "
      "Increase to prefer matching over leaving features unmatched.");
    defaults_.setMinFloat("trash_cost", 0.0);

    defaults_.setValue("mz_scale_factor", 0.0,
      "Scale factor applied to m/z values before distance computation. "
      "If 0 (default), auto-computed as max_RT_range / max_MZ_range to bring "
      "both dimensions to comparable numeric scale.");
    defaults_.setMinFloat("mz_scale_factor", 0.0);

    defaults_.setValue("normalize_intensities", "true",
      "Normalize feature intensities per map before alignment so that "
      "total intensity is equal across maps.");
    defaults_.setValidStrings("normalize_intensities", {"true", "false"});

    defaultsToParam_();
  }

  FeatureGroupingAlgorithmWNet::~FeatureGroupingAlgorithmWNet() = default;

  namespace
  {
    /// Convert a map (FeatureMap or ConsensusMap) to a Spectrum<2> with (scaled_mz, RT) positions.
    template <typename MapType>
    Spectrum<2> mapToSpectrum(const MapType& map, double mz_scale)
    {
      vector<array<double, 2>> positions;
      vector<double> intensities;
      positions.reserve(map.size());
      intensities.reserve(map.size());
      for (const auto& f : map)
      {
        positions.push_back({f.getMZ() * mz_scale, f.getRT()});
        intensities.push_back(f.getIntensity());
      }
      return Spectrum<2>(positions, intensities);
    }

    /// Compute mz_scale as max_RT_range / max_MZ_range across all maps.
    template <typename MapType>
    double computeMzScale(const vector<MapType>& maps)
    {
      double mz_min = numeric_limits<double>::max();
      double mz_max = numeric_limits<double>::lowest();
      double rt_min = numeric_limits<double>::max();
      double rt_max = numeric_limits<double>::lowest();
      for (const auto& map : maps)
      {
        for (const auto& f : map)
        {
          mz_min = min(mz_min, f.getMZ());
          mz_max = max(mz_max, f.getMZ());
          rt_min = min(rt_min, f.getRT());
          rt_max = max(rt_max, f.getRT());
        }
      }
      double mz_range = mz_max - mz_min;
      double rt_range = rt_max - rt_min;
      if (mz_range <= 0)
      {
        return 1.0; // all features at same m/z
      }
      if (rt_range <= 0)
      {
        return 1.0; // all features at same RT
      }
      return rt_range / mz_range;
    }

    DistanceMetric parseDistanceMetric(const string& s)
    {
      if (s == "L1") return DistanceMetric::L1;
      if (s == "L2") return DistanceMetric::L2;
      return DistanceMetric::LINF;
    }

  } // anonymous namespace

  template <typename MapType>
  void FeatureGroupingAlgorithmWNet::group_(const vector<MapType>& maps,
                                            ConsensusMap& out)
  {
    if (maps.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "At least two maps must be given!");
    }

    // read parameters
    DistanceMetric metric = parseDistanceMetric(param_.getValue("distance_metric").toString());
    double max_distance = param_.getValue("max_distance");
    double trash_cost = param_.getValue("trash_cost");
    double mz_scale = param_.getValue("mz_scale_factor");
    bool normalize = param_.getValue("normalize_intensities").toString() == "true";

    if (mz_scale <= 0)
    {
      mz_scale = computeMzScale(maps);
    }
    OPENMS_LOG_INFO << "FeatureGroupingAlgorithmWNet: using mz_scale_factor = " << mz_scale << endl;

    // convert all maps to Spectrum<2>
    Size n = maps.size();
    vector<Spectrum<2>> spectra;
    spectra.reserve(n);
    for (Size i = 0; i < n; ++i)
    {
      auto spec = mapToSpectrum(maps[i], mz_scale);
      if (normalize && spec.size() > 0)
      {
        spec = spec.normalized();
      }
      spectra.push_back(std::move(spec));
    }

    // Union-Find for grouping features across maps.
    // Each feature is identified as (map_index, feature_index).
    // We flatten to a single index: feature_offsets[map] + feature_idx.
    vector<Size> feature_offsets(n);
    Size total_features = 0;
    for (Size i = 0; i < n; ++i)
    {
      feature_offsets[i] = total_features;
      total_features += maps[i].size();
    }

    // map_index for each flat feature index
    vector<Size> feat_map(total_features);
    for (Size i = 0; i < n; ++i)
    {
      for (Size fi = 0; fi < maps[i].size(); ++fi)
      {
        feat_map[feature_offsets[i] + fi] = i;
      }
    }

    // Union-Find data structure
    vector<Size> parent(total_features);
    iota(parent.begin(), parent.end(), Size(0));
    vector<Size> uf_rank(total_features, 0);
    // track which map indices are present in each group (keyed by root)
    vector<set<Size>> group_maps(total_features);
    for (Size idx = 0; idx < total_features; ++idx)
    {
      group_maps[idx].insert(feat_map[idx]);
    }

    // Path-halving find (avoids std::function overhead)
    auto find = [&parent](Size x) -> Size
    {
      while (parent[x] != x)
      {
        parent[x] = parent[parent[x]]; // path halving
        x = parent[x];
      }
      return x;
    };

    // Unite two features, but only if it won't create a group
    // with multiple features from the same map.
    auto unite = [&](Size a, Size b) -> bool
    {
      a = find(a);
      b = find(b);
      if (a == b) return true; // already in same group

      // check for map conflicts
      for (Size m : group_maps[b])
      {
        if (group_maps[a].count(m))
        {
          return false; // would create duplicate map entry
        }
      }

      // merge smaller into larger (by rank)
      if (uf_rank[a] < uf_rank[b]) swap(a, b);
      parent[b] = a;
      if (uf_rank[a] == uf_rank[b]) ++uf_rank[a];
      group_maps[a].insert(group_maps[b].begin(), group_maps[b].end());
      group_maps[b].clear();
      return true;
    };

    // Pairwise alignment
    for (Size i = 0; i < n; ++i)
    {
      for (Size j = i + 1; j < n; ++j)
      {
        if (spectra[i].size() == 0 || spectra[j].size() == 0)
        {
          continue;
        }

        vector<Spectrum<2>*> theoretical = {&spectra[j]};
        WNetAligner<2> aligner(spectra[i], theoretical, metric, max_distance, trash_cost);
        aligner.set_point({1.0});

        auto [emp_ids, theo_ids] = aligner.consensus_for_target(0);

        Size accepted = 0;
        for (size_t k = 0; k < emp_ids.size(); ++k)
        {
          Size flat_i = feature_offsets[i] + static_cast<Size>(emp_ids[k]);
          Size flat_j = feature_offsets[j] + static_cast<Size>(theo_ids[k]);
          if (unite(flat_i, flat_j))
          {
            ++accepted;
          }
        }

        OPENMS_LOG_INFO << "FeatureGroupingAlgorithmWNet: maps " << i << " vs " << j
                        << ": " << emp_ids.size() << " pairs matched, "
                        << accepted << " accepted (cost = "
                        << aligner.total_cost() << ")" << endl;
      }
    }

    // Build ConsensusMap from union-find groups
    map<Size, vector<pair<Size, Size>>> groups; // root -> [(map_idx, feature_idx), ...]
    for (Size i = 0; i < n; ++i)
    {
      for (Size fi = 0; fi < maps[i].size(); ++fi)
      {
        Size root = find(feature_offsets[i] + fi);
        groups[root].emplace_back(i, fi);
      }
    }

    for (const auto& [root, members] : groups)
    {
      // Only create consensus features for groups with features from >1 map
      // (singleton = unmatched feature)
      if (members.size() < 2)
      {
        continue;
      }

      ConsensusFeature cf;
      for (const auto& [map_idx, feat_idx] : members)
      {
        cf.insert(map_idx, maps[map_idx][feat_idx]);
      }
      cf.computeConsensus();
      out.push_back(std::move(cf));
    }

    postprocess_(maps, out);
  }

  void FeatureGroupingAlgorithmWNet::group(const vector<FeatureMap>& maps,
                                           ConsensusMap& out)
  {
    group_(maps, out);
  }

  void FeatureGroupingAlgorithmWNet::group(const vector<ConsensusMap>& maps,
                                           ConsensusMap& out)
  {
    group_(maps, out);
  }

} // namespace OpenMS
