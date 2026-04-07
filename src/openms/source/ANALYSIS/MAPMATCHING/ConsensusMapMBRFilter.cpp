// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapMBRFilter.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace OpenMS
{
  namespace
  {
    /// Identity key derived from a peptide hit: modified sequence + charge.
    using IdKey = std::pair<String, Int>;

    /// Returns the highest-scoring PeptideHit from a PeptideIdentification (score-direction aware).
    /// Caller must ensure pid.getHits() is not empty.
    const PeptideHit& bestHit_(const PeptideIdentification& pid)
    {
      const auto& hits = pid.getHits();
      if (pid.isHigherScoreBetter())
      {
        return *std::max_element(hits.begin(), hits.end(), PeptideHit::ScoreLess());
      }
      else
      {
        return *std::min_element(hits.begin(), hits.end(), PeptideHit::ScoreLess());
      }
    }

    /// True if @p a is "better" than @p b under the PID's score direction.
    bool scoreBetter_(double a, double b, bool higher_is_better)
    {
      return higher_is_better ? (a > b) : (a < b);
    }

    /// Build map<map_index, (modified_seq, charge)> for a ConsensusFeature.
    /// For map indices that appear in multiple PIDs, retains the best-scoring hit.
    /// Skips PIDs with no map_index meta value or no hits.
    std::map<UInt64, IdKey>
    buildMapIndexToIdKey_(const ConsensusFeature& cf)
    {
      std::map<UInt64, IdKey> result;
      // remember the score of the currently-stored best hit per map_index
      std::map<UInt64, std::pair<double, bool>> best_score;

      for (const PeptideIdentification& pid : cf.getPeptideIdentifications())
      {
        if (pid.getHits().empty()) { continue; }
        if (!pid.metaValueExists("map_index")) { continue; }
        const UInt64 map_index = static_cast<UInt64>(pid.getMetaValue("map_index"));

        const PeptideHit& hit = bestHit_(pid);
        const double score = hit.getScore();
        const bool higher = pid.isHigherScoreBetter();

        auto it = best_score.find(map_index);
        if (it == best_score.end() || scoreBetter_(score, it->second.first, higher))
        {
          result[map_index] = std::make_pair(hit.getSequence().toString(), hit.getCharge());
          best_score[map_index] = std::make_pair(score, higher);
        }
      }
      return result;
    }
  } // unnamed namespace

  ConsensusMapMBRFilter::ConsensusMapMBRFilter() :
    DefaultParamHandler("ConsensusMapMBRFilter")
  {
    defaults_.setValue("k_sigma", 3.0,
                       "Multiplier on the per-pair robust sigma (1.4826*MAD) for the rejection threshold.");
    defaults_.setMinFloat("k_sigma", 0.0);

    defaults_.setValue("min_anchors_per_pair", 5,
                       "Minimum number of anchor pairs (same modified sequence + charge) required to use a per-pair distribution. "
                       "Pairs below this threshold use the pooled fallback if enabled.");
    defaults_.setMinInt("min_anchors_per_pair", 2);

    defaults_.setValue("reject_fraction", 0.5,
                       "Fraction of ID-bearing handles in the same ConsensusFeature that must reject the unidentified handle "
                       "for it to be removed.");
    defaults_.setMinFloat("reject_fraction", 0.0);
    defaults_.setMaxFloat("reject_fraction", 1.0);

    defaults_.setValue("min_pair_sigma", 1.0,
                       "Floor (in seconds) for the per-pair sigma derived from MAD. Avoids zero-MAD pathologies.");
    defaults_.setMinFloat("min_pair_sigma", 0.0);

    defaults_.setValue("use_pooled_fallback", "true",
                       "If 'true', map pairs with fewer than 'min_anchors_per_pair' anchors fall back to the pooled distribution.");
    defaults_.setValidStrings("use_pooled_fallback", {"true", "false"});

    defaultsToParam_();
  }

  ConsensusMapMBRFilter::Report ConsensusMapMBRFilter::filter(ConsensusMap& cmap)
  {
    Report report;

    const double k_sigma              = (double)param_.getValue("k_sigma");
    const Size   min_anchors_per_pair = (Size)(int)param_.getValue("min_anchors_per_pair");
    const double reject_fraction      = (double)param_.getValue("reject_fraction");
    const double min_pair_sigma       = (double)param_.getValue("min_pair_sigma");
    const bool   use_pooled_fallback  = param_.getValue("use_pooled_fallback").toBool();

    constexpr double MAD_TO_SIGMA = 1.4826;

    // Cache map_index → (mod_seq, charge) per ConsensusFeature so both passes share the same lookup.
    std::vector<std::map<UInt64, IdKey>> cf_id_maps;
    cf_id_maps.reserve(cmap.size());
    for (const ConsensusFeature& cf : cmap)
    {
      cf_id_maps.push_back(buildMapIndexToIdKey_(cf));
    }

    // Computes median + robust sigma (1.4826 * MAD, floored to min_pair_sigma) on @p vec (modified in-place).
    auto computeStats = [&](std::vector<double>& vec, PairStats& s)
    {
      s.n        = vec.size();
      s.median   = Math::median(vec.begin(), vec.end(), false);
      const double mad = Math::MAD(vec.begin(), vec.end(), s.median);
      s.mad_sigma = std::max(min_pair_sigma, MAD_TO_SIGMA * mad);
    };

    // ----- Pass 1: collect anchor residuals per (smaller, larger) map pair -----
    std::map<PairKey, std::vector<double>> raw_residuals;

    for (Size cf_idx = 0; cf_idx < cmap.size(); ++cf_idx)
    {
      const ConsensusFeature& cf = cmap[cf_idx];
      if (cf.size() < 2) { continue; }
      const auto& map_to_idkey = cf_id_maps[cf_idx];
      if (map_to_idkey.size() < 2) { continue; }

      // The handle set is ordered by IndexLess (map_index ascending), so the
      // outer/inner double-loop with j>i naturally gives map_i <= map_j.
      const auto& handles = cf.getFeatures();
      for (auto it_i = handles.begin(); it_i != handles.end(); ++it_i)
      {
        const UInt64 mi = it_i->getMapIndex();
        auto k_i = map_to_idkey.find(mi);
        if (k_i == map_to_idkey.end()) { continue; }

        auto it_j = it_i;
        ++it_j;
        for (; it_j != handles.end(); ++it_j)
        {
          const UInt64 mj = it_j->getMapIndex();
          if (mj <= mi) { continue; } // belt-and-braces: equal map_index would be a duplicate handle
          auto k_j = map_to_idkey.find(mj);
          if (k_j == map_to_idkey.end()) { continue; }
          if (k_i->second != k_j->second) { continue; }

          // signed: rt(map_i) - rt(map_j), with map_i < map_j
          raw_residuals[{mi, mj}].push_back(it_i->getRT() - it_j->getRT());
        }
      }
    }

    // Compute per-pair stats
    for (auto& kv : raw_residuals)
    {
      PairStats s;
      computeStats(kv.second, s);
      report.pair_stats[kv.first] = s;
    }

    // Compute pooled stats from pairs that already have enough anchors. Centering each
    // residual by its own pair median measures spread independent of per-pair bias.
    {
      std::vector<double> pooled_centered;
      for (const auto& kv : raw_residuals)
      {
        if (kv.second.size() < min_anchors_per_pair) { continue; }
        const double pm = report.pair_stats[kv.first].median;
        for (double r : kv.second)
        {
          pooled_centered.push_back(r - pm);
        }
      }
      if (!pooled_centered.empty())
      {
        computeStats(pooled_centered, report.pooled);
      }
    }

    // ----- Pass 2: filter unidentified handles -----
    for (Size cf_idx = 0; cf_idx < cmap.size(); ++cf_idx)
    {
      ConsensusFeature& cf = cmap[cf_idx];
      if (cf.size() < 2) { continue; }
      const auto& map_to_idkey = cf_id_maps[cf_idx];
      if (map_to_idkey.empty()) { continue; }

      auto isAnchor = [&](UInt64 mi) { return map_to_idkey.find(mi) != map_to_idkey.end(); };

      Size id_count = 0;
      for (const auto& fh : cf.getFeatures())
      {
        if (isAnchor(fh.getMapIndex())) { ++id_count; }
      }
      if (id_count == 0) { continue; }
      if (id_count == cf.size()) { continue; }

      std::set<UInt64> remove_map_indices;
      std::map<PairKey, Size> per_pair_remove;

      for (const auto& h_t : cf.getFeatures())
      {
        const UInt64 mt = h_t.getMapIndex();
        if (isAnchor(mt)) { continue; }

        Size rejections = 0;
        Size comparisons = 0;
        PairKey worst_pair{0, 0};
        double  worst_z = -std::numeric_limits<double>::infinity();

        for (const auto& h_a : cf.getFeatures())
        {
          const UInt64 ma = h_a.getMapIndex();
          if (!isAnchor(ma)) { continue; }

          PairKey key = (mt < ma) ? PairKey{mt, ma} : PairKey{ma, mt};
          // signed_dt aligned to (smaller, larger): rt(smaller) - rt(larger)
          double signed_dt = (mt < ma) ? (h_t.getRT() - h_a.getRT())
                                       : (h_a.getRT() - h_t.getRT());

          double median = 0.0;
          double sigma  = 0.0;
          auto it = report.pair_stats.find(key);
          if (it != report.pair_stats.end() && it->second.n >= min_anchors_per_pair)
          {
            median = it->second.median;
            sigma  = it->second.mad_sigma;
          }
          else if (use_pooled_fallback && report.pooled.n >= min_anchors_per_pair)
          {
            median = 0.0;
            sigma  = report.pooled.mad_sigma;
          }
          else
          {
            continue;
          }

          ++comparisons;
          const double z = std::abs(signed_dt - median) / sigma;
          if (z > k_sigma) { ++rejections; }
          if (z > worst_z) { worst_z = z; worst_pair = key; }
        }

        if (comparisons == 0) { continue; }
        const double frac = static_cast<double>(rejections) / static_cast<double>(comparisons);
        if (frac >= reject_fraction && rejections > 0)
        {
          remove_map_indices.insert(mt);
          per_pair_remove[worst_pair] += 1;
        }
      }

      if (remove_map_indices.empty()) { continue; }

      ConsensusFeature::HandleSetType kept;
      for (const auto& fh : cf.getFeatures())
      {
        if (remove_map_indices.find(fh.getMapIndex()) == remove_map_indices.end())
        {
          kept.insert(fh);
        }
      }
      // Should never become empty (id_count >= 1 always survives)
      if (kept.empty()) { continue; }

      cf.setFeatures(std::move(kept));
      cf.computeConsensus();

      report.handles_removed += remove_map_indices.size();
      for (const auto& kv : per_pair_remove)
      {
        report.pair_stats[kv.first].removals += kv.second;
      }
    }

    return report;
  }

} // namespace OpenMS
