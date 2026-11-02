// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "NuXLRescoring.h"

#include <OpenMS/METADATA/PeptideIdentification.h>

namespace OpenMS
{
  std::array<std::vector<Size>, 2> NuXLRescoring::confidentHitsPerHalf(const PeptideIdentificationList& peptide_ids, double fdr,
                                                                       std::array<std::array<Size, 2>, 2>& counts)
  {
    struct Row { double score; bool decoy; Size id; };
    std::array<std::array<std::vector<Row>, 2>, 2> rows; // [half][class]
    for (Size i = 0; i < peptide_ids.size(); ++i)
    {
      const PeptideIdentification& pi = peptide_ids[i];
      if (pi.getHits().empty() || !pi.metaValueExists("scan_index")) continue;
      const PeptideHit& top = pi.getHits()[0];
      const Size half = static_cast<Size>(static_cast<int>(pi.getMetaValue("scan_index"))) % 2;
      const Size cls = static_cast<int>(top.getMetaValue("NuXL:isXL", 0)) != 0 ? 1 : 0;
      rows[half][cls].push_back({top.getScore(), top.isDecoy(), i});
    }
    std::array<std::vector<Size>, 2> training;
    counts = {};
    for (Size half = 0; half < 2; ++half)
    {
      for (Size cls = 0; cls < 2; ++cls)
      {
        std::vector<Row>& r = rows[half][cls];
        std::sort(r.begin(), r.end(), [](const Row& a, const Row& b)
        {
          if (a.score != b.score) return a.score > b.score;
          if (a.decoy != b.decoy) return a.decoy;
          return a.id < b.id;
        });
        std::vector<double> q(r.size(), 1.0);
        Size targets = 0, decoys = 0;
        for (Size k = 0; k < r.size(); ++k)
        {
          if (r[k].decoy) { ++decoys; } else { ++targets; }
          q[k] = static_cast<double>(decoys + 1) / static_cast<double>(std::max<Size>(targets, 1));
        }
        for (Size k = r.size(); k-- > 1;) { q[k - 1] = std::min(q[k - 1], q[k]); }
        for (Size k = 0; k < r.size(); ++k)
        {
          if (r[k].decoy || q[k] > fdr) continue;
          training[half].push_back(r[k].id);
          ++counts[half][cls];
        }
      }
    }
    return training;
  }

  void NuXLRescoring::rescoreHits(PeptideIdentificationList& peptide_ids)
  {
    for (auto& pi : peptide_ids)
    {
      for (auto& ph : pi.getHits())
      {
        double score = ph.getMetaValue("NuXL:score", ph.getScore());
        const bool is_xl = static_cast<int>(ph.getMetaValue("NuXL:isXL", 0)) != 0;
        if (is_xl) { score += DELTA_WEIGHT * static_cast<double>(ph.getMetaValue("NuXL:delta_score", 0.0)); }
        score += ION_LLR_WEIGHT * static_cast<double>(ph.getMetaValue("NuXL:ion_llr", 0.0));
        ph.setScore(score);
      }
      pi.sort();
      std::vector<PeptideHit>& phs = pi.getHits();
      for (Size r = 0; r != phs.size(); ++r) { phs[r].setRank(static_cast<int>(r)); }
    }
  }
} // namespace OpenMS
