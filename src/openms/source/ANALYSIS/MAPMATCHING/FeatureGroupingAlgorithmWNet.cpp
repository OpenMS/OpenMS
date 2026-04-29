// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Michal Startek $
// $Authors: Michal Startek $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmWNet.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <wnetalign/aligner.hpp>
#include <wnetalign/spectrum.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>

using namespace std;

namespace OpenMS
{

  FeatureGroupingAlgorithmWNet::FeatureGroupingAlgorithmWNet() :
    FeatureGroupingAlgorithm()
  {
    setName("FeatureGroupingAlgorithmWNet");

    defaults_.setValue("distance_metric", "LINF", "Distance metric for comparing feature positions");
    defaults_.setValidStrings("distance_metric", {"L1", "L2", "LINF"});
    defaults_.addTag("distance_metric", "advanced");

    defaults_.setValue("max_rt_shift", 100.0,
      "Maximum allowed RT difference (in seconds) for matching features.");
    defaults_.setMinFloat("max_rt_shift", 0.0);

    defaults_.setValue("mz_unit", "ppm",
      "Unit for the m/z tolerance. 'Da' uses max_mz_shift_da; "
      "'ppm' uses max_mz_shift_ppm and log-transforms the m/z axis.");
    defaults_.setValidStrings("mz_unit", {"Da", "ppm"});

    defaults_.setValue("max_mz_shift_da", 0.3,
      "Maximum allowed m/z difference in Daltons. Used when mz_unit is 'Da'.");
    defaults_.setMinFloat("max_mz_shift_da", 0.0);

    defaults_.setValue("max_mz_shift_ppm", 10.0,
      "Maximum allowed m/z difference in ppm. Used when mz_unit is 'ppm'.");
    defaults_.setMinFloat("max_mz_shift_ppm", 0.0);

    defaults_.setValue("trash_cost", 0.0,
      "Cost of leaving a feature unmatched (in seconds). Set to 0 to "
      "auto-derive from max_rt_shift (a match at the boundary costs the same "
      "as no match). Increase to prefer matching over leaving features unmatched.");
    defaults_.setMinFloat("trash_cost", 0.0);
    defaults_.addTag("trash_cost", "advanced");

    defaults_.setValue("normalize_intensities", "true",
      "Normalize feature intensities per map before alignment so that "
      "total intensity is equal across maps.");
    defaults_.setValidStrings("normalize_intensities", {"true", "false"});

    defaults_.setValue("decharge_mz", "false",
      "Convert observed m/z values to their singly-charged (z=1) equivalent before "
      "alignment. Features with charge z>1 are projected to [M+H]+ m/z, making "
      "features of the same compound but different charge states directly comparable. "
      "If a feature carries a 'dc_charge_adduct_mass' MetaValue (set by "
      "MetaboliteFeatureDeconvolution), the stored adduct mass is used; otherwise "
      "pure H+ adducts are assumed. Features with unknown charge (z=0) are left "
      "unchanged.");
    defaults_.setValidStrings("decharge_mz", {"true", "false"});

    defaultsToParam_();
  }

  FeatureGroupingAlgorithmWNet::~FeatureGroupingAlgorithmWNet() = default;

  namespace
  {
    /// Transform m/z value: identity for Da, log for ppm.
    /// In ppm mode, a multiplicative shift becomes additive in log-space:
    /// log(mz * factor) = log(mz) + log(factor), so the distance is
    /// independent of the absolute m/z value.
    inline double transformMz(double mz, bool use_ppm)
    {
      return use_ppm ? log(mz) : mz;
    }

    /// Convert observed m/z (at charge z) to the equivalent singly-charged [M+H]+ m/z.
    /// If the feature carries a dc_charge_adduct_mass MetaValue (from
    /// MetaboliteFeatureDeconvolution), that adduct mass is used to compute the neutral
    /// mass; otherwise pure H+ adducts are assumed.
    /// Features with charge 0 (unknown) or charge 1 are returned unchanged.
    template <typename FeatureType>
    double dechargeMz(const FeatureType& f)
    {
      const int z = f.getCharge();
      if (z <= 1) return f.getMZ();
      const double mz = f.getMZ();
      double neutral_mass;
      if (f.metaValueExists("dc_charge_adduct_mass"))
      {
        neutral_mass = mz * z - (double)f.getMetaValue("dc_charge_adduct_mass");
      }
      else
      {
        neutral_mass = mz * z - z * Constants::PROTON_MASS_U;
      }
      return neutral_mass + Constants::PROTON_MASS_U;
    }

    /// Convert a map (FeatureMap or ConsensusMap) to a Spectrum<2> with (scaled_mz, RT) positions.
    template <typename MapType>
    Spectrum<2> mapToSpectrum(const MapType& map, double mz_scale, bool use_ppm, bool decharge)
    {
      vector<array<double, 2>> positions;
      vector<double> intensities;
      positions.reserve(map.size());
      intensities.reserve(map.size());
      for (const auto& f : map)
      {
        double mz = decharge ? dechargeMz(f) : f.getMZ();
        positions.push_back({transformMz(mz, use_ppm) * mz_scale, f.getRT()});
        intensities.push_back(f.getIntensity());
      }
      return Spectrum<2>(positions, intensities);
    }

    DistanceMetric parseDistanceMetric(const string& s)
    {
      if (s == "L1") return DistanceMetric::L1;
      if (s == "L2") return DistanceMetric::L2;
      return DistanceMetric::LINF;
    }

    struct AlignParams
    {
      DistanceMetric metric;
      double max_rt_shift;
      double trash_cost;
      double mz_scale;
      double max_mz_shift; // for logging
      double max_distance;
      bool use_ppm;
      bool normalize;
      bool decharge;
    };

    AlignParams buildAlignParams(const Param& param)
    {
      AlignParams p;
      p.metric    = parseDistanceMetric(param.getValue("distance_metric").toString());
      p.max_rt_shift = param.getValue("max_rt_shift");
      p.trash_cost   = param.getValue("trash_cost");
      p.normalize    = param.getValue("normalize_intensities").toString() == "true";
      p.use_ppm      = param.getValue("mz_unit").toString() == "ppm";
      p.decharge     = param.getValue("decharge_mz").toString() == "true";

      if (p.trash_cost <= 0) p.trash_cost = p.max_rt_shift;

      // Compute mz_scale so that the m/z tolerance maps to max_rt_shift in the
      // scaled space. After scaling, both dimensions are in seconds and
      // max_distance = max_rt_shift.
      double mz_shift_native;
      if (p.use_ppm)
      {
        p.max_mz_shift = param.getValue("max_mz_shift_ppm");
        mz_shift_native = std::log1p(p.max_mz_shift * 1e-6);
      }
      else
      {
        p.max_mz_shift = param.getValue("max_mz_shift_da");
        mz_shift_native = p.max_mz_shift;
      }
      // When max_mz_shift == 0 the m/z axis must be a strict constraint (no
      // m/z difference is acceptable). Setting mz_scale = ∞ maps any non-zero
      // m/z difference to an infinite scaled distance, rejecting such pairs.
      p.mz_scale = (mz_shift_native > 0) ? p.max_rt_shift / mz_shift_native
                                          : std::numeric_limits<double>::infinity();

      // After mz_scale, both axes are in "RT-equivalent seconds".
      // max_distance must reflect the chosen norm so that a feature exactly at
      // the per-axis threshold [max_rt_shift, max_rt_shift] lands on the boundary:
      //   LINF: max(|Δmz|, |Δrt|) ≤ max_rt_shift  →  max_distance = max_rt_shift
      //   L1:   |Δmz| + |Δrt| ≤ 2·max_rt_shift    →  max_distance = 2·max_rt_shift
      //   L2:   √(Δmz²+Δrt²) ≤ √2·max_rt_shift    →  max_distance = √2·max_rt_shift
      switch (p.metric)
      {
        case DistanceMetric::L1: p.max_distance = 2.0 * p.max_rt_shift; break;
        case DistanceMetric::L2: p.max_distance = std::sqrt(2.0) * p.max_rt_shift; break;
        default:                 p.max_distance = p.max_rt_shift; break; // LINF
      }
      return p;
    }

    template <typename MapType>
    vector<Spectrum<2>> convertAndNormalizeSpectra(const vector<MapType>& maps,
                                                    double mz_scale, bool use_ppm,
                                                    bool decharge, bool normalize)
    {
      vector<Spectrum<2>> spectra;
      spectra.reserve(maps.size());
      for (const auto& m : maps)
      {
        auto spec = mapToSpectrum(m, mz_scale, use_ppm, decharge);
        if (normalize && spec.size() > 0) spec = spec.normalized();
        spectra.push_back(std::move(spec));
      }
      return spectra;
    }

    struct Candidate { Size flat_i, flat_j; double cost; };

    // spectra is non-const because WNetAligner requires non-const Spectrum<2>*
    vector<Candidate> collectCandidates(vector<Spectrum<2>>& spectra,
                                         const vector<Size>& feature_offsets,
                                         DistanceMetric metric,
                                         double max_distance,
                                         double trash_cost)
    {
      vector<Candidate> candidates;
      const Size n = spectra.size();
      for (Size i = 0; i < n; ++i)
      {
        for (Size j = i + 1; j < n; ++j)
        {
          if (spectra[i].size() == 0 || spectra[j].size() == 0) continue;

          vector<Spectrum<2>*> theoretical = {&spectra[j]};
          WNetAligner<2> aligner(spectra[i], theoretical, metric, max_distance, trash_cost);
          aligner.set_point({1.0}); // weight=1: solve a single unit-weight consensus

          auto [emp_ids, theo_ids] = aligner.consensus_for_target(0);
          // Average per-pair cost: globally cheaper matches are merged first
          // regardless of map order, so this serves as a confidence proxy.
          double pair_cost = emp_ids.empty() ? 0.0
                           : aligner.total_cost() / static_cast<double>(emp_ids.size());

          OPENMS_LOG_INFO << "FeatureGroupingAlgorithmWNet: maps " << i << " vs " << j
                          << ": " << emp_ids.size() << " candidate pairs (total cost = "
                          << aligner.total_cost() << ")" << '\n';

          for (Size k = 0; k < emp_ids.size(); ++k)
          {
            candidates.push_back({
              feature_offsets[i] + static_cast<Size>(emp_ids[k]),
              feature_offsets[j] + static_cast<Size>(theo_ids[k]),
              pair_cost
            });
          }
        }
      }
      return candidates;
    }

    /// Union-Find with map-conflict prevention.
    /// Ensures no consensus group contains two features from the same map.
    class UnionFind
    {
    public:
      UnionFind(Size total, const vector<Size>& feat_map)
        : parent_(total), rank_(total, 0), group_maps_(total)
      {
        iota(parent_.begin(), parent_.end(), Size(0));
        for (Size idx = 0; idx < total; ++idx)
          group_maps_[idx].insert(feat_map[idx]);
      }

      Size find(Size x)
      {
        while (parent_[x] != x)
        {
          parent_[x] = parent_[parent_[x]]; // path halving
          x = parent_[x];
        }
        return x;
      }

      /// Returns true if merged (or already in same group), false on map conflict.
      bool unite(Size a, Size b)
      {
        a = find(a);
        b = find(b);
        if (a == b) return true;

        for (Size m : group_maps_[b])
        {
          if (group_maps_[a].count(m)) return false; // would create duplicate map entry
        }

        if (rank_[a] < rank_[b]) swap(a, b);
        parent_[b] = a;
        if (rank_[a] == rank_[b]) ++rank_[a];
        group_maps_[a].insert(group_maps_[b].begin(), group_maps_[b].end());
        group_maps_[b].clear();
        return true;
      }

    private:
      vector<Size> parent_;
      vector<Size> rank_;
      vector<set<Size>> group_maps_;
    };

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

    AlignParams p = buildAlignParams(param_);

    OPENMS_LOG_INFO << "FeatureGroupingAlgorithmWNet: max_mz_shift = " << p.max_mz_shift
                    << (p.use_ppm ? " ppm" : " Da")
                    << ", max_rt_shift = " << p.max_rt_shift << " s"
                    << ", mz_scale = " << p.mz_scale << '\n';

    Size n = maps.size();
    auto spectra = convertAndNormalizeSpectra(maps, p.mz_scale, p.use_ppm, p.decharge, p.normalize);

    // Flat feature indexing: feature_offsets[map_i] + feature_idx_in_map
    vector<Size> feature_offsets(n);
    Size total_features = 0;
    for (Size i = 0; i < n; ++i)
    {
      feature_offsets[i] = total_features;
      total_features += maps[i].size();
    }

    vector<Size> feat_map(total_features);
    for (Size i = 0; i < n; ++i)
    {
      for (Size fi = 0; fi < maps[i].size(); ++fi)
      {
        feat_map[feature_offsets[i] + fi] = i;
      }
    }

    // Collect all pairwise candidates, sort best-first, merge greedily.
    auto candidates = collectCandidates(spectra, feature_offsets,
                                        p.metric, p.max_distance, p.trash_cost);

    sort(candidates.begin(), candidates.end(),
         [](const Candidate& a, const Candidate& b)
         {
           // Primary: ascending cost. Tie-break by ordered endpoint pair so
           // equal-cost candidates are always processed in the same order.
           auto key = [](const Candidate& c) {
             return std::make_tuple(c.cost,
                                    std::min(c.flat_i, c.flat_j),
                                    std::max(c.flat_i, c.flat_j));
           };
           return key(a) < key(b);
         });

    UnionFind uf(total_features, feat_map);
    Size accepted = 0;
    for (const auto& c : candidates)
    {
      if (uf.unite(c.flat_i, c.flat_j)) ++accepted;
    }
    OPENMS_LOG_INFO << "FeatureGroupingAlgorithmWNet: " << candidates.size()
                    << " total candidates, " << accepted
                    << " accepted after globally sorted merge\n";

    // Build ConsensusMap from union-find groups
    map<Size, vector<pair<Size, Size>>> groups; // root -> [(map_idx, feature_idx), ...]
    for (Size i = 0; i < n; ++i)
    {
      for (Size fi = 0; fi < maps[i].size(); ++fi)
      {
        Size root = uf.find(feature_offsets[i] + fi);
        groups[root].emplace_back(i, fi);
      }
    }

    for (const auto& [root, members] : groups)
    {
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
