// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIACNNRTWindow.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

namespace OpenMS::OpenDIA
{
  std::vector<bool> crossfitRTWindowMask(const std::vector<CNNQValueInput>& inputs, const std::vector<double>& apex_rt,
                                         const std::vector<double>& library_rt, std::vector<RTWindowFold>& report,
                                         std::size_t anchors_per_fold)
  {
    if (apex_rt.size() != inputs.size() || library_rt.size() != inputs.size())
    {
      throw std::invalid_argument("RT window inputs must have one apex and library RT per window");
    }
    // Best window per precursor (ties: smaller sample ID), as in the q-value winners.
    std::map<std::string, std::size_t> best;
    for (std::size_t i = 0; i < inputs.size(); ++i)
    {
      if (!std::isfinite(apex_rt[i]) || !std::isfinite(library_rt[i]) || !std::isfinite(inputs[i].score))
      {
        throw std::invalid_argument("RT window inputs must be finite");
      }
      const auto [it, inserted] = best.emplace(inputs[i].precursor_id, i);
      const auto& previous = inputs[it->second];
      if (!inserted && (inputs[i].score > previous.score ||
                        (inputs[i].score == previous.score && inputs[i].sample_id < previous.sample_id)))
      {
        it->second = i;
      }
    }
    // Top target precursors per fold by (score, precursor ID).
    std::map<int, std::vector<std::size_t>> targets;
    std::map<int, std::size_t> precursors;
    for (const auto& entry : best)
    {
      const auto& row = inputs[entry.second];
      ++precursors[row.fold];
      if (!row.is_decoy)
      {
        targets[row.fold].push_back(entry.second);
      }
    }
    if (precursors.size() < 2)
    {
      throw std::invalid_argument("Out-of-fold RT calibration needs at least two folds");
    }
    for (auto& entry : targets)
    {
      auto& indices = entry.second;
      std::sort(indices.begin(), indices.end(), [&](std::size_t a, std::size_t b)
      {
        if (inputs[a].score != inputs[b].score)
        {
          return inputs[a].score > inputs[b].score;
        }
        return inputs[a].precursor_id > inputs[b].precursor_id;
      });
      indices.resize(std::min(indices.size(), anchors_per_fold));
    }

    report.clear();
    std::map<int, std::size_t> fold_index;
    for (const auto& entry : precursors)
    {
      RTWindowFold summary;
      summary.fold = entry.first;
      summary.precursors = entry.second;
      std::vector<std::pair<double, double>> anchors;
      for (const auto& other : targets)
      {
        if (other.first == entry.first)
        {
          continue;
        }
        for (const std::size_t i : other.second)
        {
          anchors.emplace_back(library_rt[i], apex_rt[i]);
        }
      }
      summary.map = fitRTMap(std::move(anchors));
      summary.half_window = std::max(4.0 * summary.map.sigma, 30.0);
      fold_index[entry.first] = report.size();
      report.push_back(std::move(summary));
    }

    std::vector<bool> mask(inputs.size(), false);
    std::map<std::string, bool> kept_precursor;
    for (std::size_t i = 0; i < inputs.size(); ++i)
    {
      const auto& fold = report[fold_index.at(inputs[i].fold)];
      mask[i] = std::fabs(apex_rt[i] - fold.map(library_rt[i])) <= fold.half_window;
      if (mask[i] && kept_precursor.emplace(inputs[i].precursor_id, true).second)
      {
        ++report[fold_index.at(inputs[i].fold)].kept_precursors;
      }
    }
    return mask;
  }
}
