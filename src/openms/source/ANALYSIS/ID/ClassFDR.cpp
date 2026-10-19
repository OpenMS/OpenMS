// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ClassFDR.h>

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <numeric>
#include <tuple>

namespace OpenMS
{
  const std::string ClassFDR::CLASS_Q_VALUE = "class_q-value";

  namespace
  {
    constexpr Size MIN_GUARD_DECOYS = 20;

    struct Ref
    {
      Size id;
      Size hit;
      double rank;
    };

    /// q = (decoys + 1) / targets within a class, made monotone; at equal rank a decoy counts first
    std::vector<double> classTDC(const std::vector<Ref>& targets, const std::vector<double>& decoy_ranks)
    {
      // rows: (rank, is_target, index into targets)
      std::vector<std::tuple<double, bool, Size>> rows;
      rows.reserve(targets.size() + decoy_ranks.size());
      for (Size i = 0; i < targets.size(); ++i) rows.emplace_back(targets[i].rank, true, i);
      for (double d : decoy_ranks) rows.emplace_back(d, false, 0);
      // descending rank; at equal rank decoys first (false < true)
      std::sort(rows.begin(), rows.end(), [](const auto& a, const auto& b)
      {
        if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) > std::get<0>(b);
        return std::get<1>(a) < std::get<1>(b);
      });
      std::vector<double> fdr(rows.size());
      Size t = 0, d = 0;
      for (Size r = 0; r < rows.size(); ++r)
      {
        if (std::get<1>(rows[r])) ++t; else ++d;
        fdr[r] = std::min(1.0, static_cast<double>(d + 1) / static_cast<double>(std::max<Size>(t, 1)));
      }
      std::vector<double> q(targets.size(), 1.0);
      double running = 1.0;
      for (Size r = rows.size(); r-- > 0;)
      {
        running = std::min(running, fdr[r]);
        if (std::get<1>(rows[r])) q[std::get<2>(rows[r])] = running;
      }
      return q;
    }
  }

  std::vector<double> ClassFDR::storeyQValues(const std::vector<double>& p, double& pi0)
  {
    const Size m = p.size();
    pi0 = 1.0;
    if (m == 0) return {};
    const Size above = static_cast<Size>(std::count_if(p.begin(), p.end(), [](double x) { return x > 0.5; }));
    pi0 = std::min(1.0, static_cast<double>(above) / (0.5 * static_cast<double>(m)));
    std::vector<Size> order(m);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&p](Size a, Size b) { return p[a] < p[b]; });
    std::vector<double> q(m, 1.0);
    double running = 1.0;
    for (Size r = m; r-- > 0;)
    {
      const Size i = order[r];
      running = std::min(running, pi0 * p[i] * static_cast<double>(m) / static_cast<double>(r + 1));
      q[i] = std::min(running, 1.0);
    }
    return q;
  }

  std::vector<ClassFDR::ClassSummary> ClassFDR::annotateClassQValues(PeptideIdentificationList& ids,
                                                                     const std::function<std::string(const PeptideHit&)>& class_of,
                                                                     bool use_all_hits)
  {
    std::map<std::string, std::pair<std::vector<Ref>, std::vector<double>>> by_class;
    std::vector<double> pooled;
    for (Size i = 0; i < ids.size(); ++i)
    {
      const double sign = ids[i].isHigherScoreBetter() ? 1.0 : -1.0;
      const auto& hits = ids[i].getHits();
      if (hits.empty()) continue;
      Size best = 0;
      for (Size j = 1; j < hits.size(); ++j)
      {
        if (sign * hits[j].getScore() > sign * hits[best].getScore()) best = j;
      }
      for (Size j = 0; j < hits.size(); ++j)
      {
        if (!use_all_hits && j != best) continue;
        const double rank = sign * hits[j].getScore();
        auto& entry = by_class[class_of(hits[j])];
        if (hits[j].isDecoy())
        {
          entry.second.push_back(rank);
          pooled.push_back(rank);
        }
        else
        {
          entry.first.push_back({i, j, rank});
        }
      }
    }
    std::sort(pooled.begin(), pooled.end());
    const double n = static_cast<double>(pooled.size());
    // share of pooled decoys ranking at least s, with +1
    auto p_of = [&pooled, n](double s)
    {
      const double ge = static_cast<double>(pooled.end() - std::lower_bound(pooled.begin(), pooled.end(), s));
      return (ge + 1.0) / (n + 1.0);
    };

    const bool single = by_class.size() < 2;
    std::vector<ClassSummary> summaries;
    for (auto& [name, rows] : by_class)
    {
      auto& [targets, decoys] = rows;
      ClassSummary summary;
      summary.name = name;
      summary.targets = targets.size();
      summary.decoys = decoys.size();
      bool heavier = false;
      if (decoys.size() >= MIN_GUARD_DECOYS)
      {
        const double m = static_cast<double>(decoys.size());
        const double k = static_cast<double>(std::count_if(decoys.begin(), decoys.end(), [&](double s) { return p_of(s) <= 0.05; }));
        heavier = (k - 0.05 * m) / std::sqrt(m * 0.05 * 0.95) > 2.33;
      }
      summary.method = (decoys.size() < MIN_GUARD_DECOYS || heavier) ? Method::CLASS_TDC : Method::POOLED_NULL;
      if (!single)
      {
        std::vector<double> q;
        if (summary.method == Method::POOLED_NULL)
        {
          std::vector<double> p;
          p.reserve(targets.size());
          for (const Ref& t : targets) p.push_back(p_of(t.rank));
          q = storeyQValues(p, summary.pi0);
        }
        else
        {
          q = classTDC(targets, decoys);
        }
        for (Size k = 0; k < targets.size(); ++k)
        {
          ids[targets[k].id].getHits()[targets[k].hit].setMetaValue(CLASS_Q_VALUE, q[k]);
        }
      }
      summaries.push_back(std::move(summary));
    }
    std::sort(summaries.begin(), summaries.end(), [](const ClassSummary& a, const ClassSummary& b)
    {
      if (a.targets != b.targets) return a.targets > b.targets;
      return a.name < b.name;
    });
    return summaries;
  }

  Size ClassFDR::raiseToClassQValues(PeptideIdentificationList& ids)
  {
    Size raised = 0;
    for (auto& id : ids)
    {
      for (auto& hit : id.getHits())
      {
        if (!hit.metaValueExists(CLASS_Q_VALUE)) continue;
        const double class_q = hit.getMetaValue(CLASS_Q_VALUE);
        if (class_q > hit.getScore())
        {
          hit.setScore(class_q);
          ++raised;
        }
      }
    }
    return raised;
  }

  std::string ClassFDR::classFromMetaValue(const PeptideHit& hit, const std::string& name)
  {
    return hit.metaValueExists(name) ? hit.getMetaValue(name).toString() : std::string("NA");
  }

  std::string ClassFDR::classFromModifications(const PeptideHit& hit, const std::vector<std::string>& ignore)
  {
    std::vector<std::string> mods;
    const AASequence& seq = hit.getSequence();
    auto add = [&](const ResidueModification* mod)
    {
      if (mod == nullptr) return;
      if (std::find(ignore.begin(), ignore.end(), mod->getFullId()) != ignore.end()) return;
      mods.push_back(mod->getId());
    };
    add(seq.getNTerminalModification());
    add(seq.getCTerminalModification());
    for (Size i = 0; i < seq.size(); ++i) add(seq[i].getModification());
    if (mods.empty()) return "unmodified";
    std::sort(mods.begin(), mods.end());
    mods.erase(std::unique(mods.begin(), mods.end()), mods.end());
    std::string cls;
    for (const auto& m : mods) cls += (cls.empty() ? "" : "+") + m;
    return cls;
  }
} // namespace OpenMS
