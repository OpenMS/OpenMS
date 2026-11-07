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
  namespace
  {
    /// Median of a nonempty range; even sizes average the two middle values.
    double median_(std::vector<double> values)
    {
      const std::size_t mid = values.size() / 2;
      std::nth_element(values.begin(), values.begin() + mid, values.end());
      const double upper = values[mid];
      if (values.size() % 2 == 1)
      {
        return upper;
      }
      return (*std::max_element(values.begin(), values.begin() + mid) + upper) / 2.0;
    }

    /// Knots from 20 equal-count bins of anchors sorted by (library RT, run RT).
    RTMap binnedKnots_(const std::vector<std::pair<double, double>>& points)
    {
      RTMap map;
      const std::size_t n = points.size();
      for (std::size_t bin = 0; bin < 20; ++bin)
      {
        const std::size_t begin = bin * n / 20;
        const std::size_t end = (bin + 1) * n / 20;
        if (begin == end)
        {
          continue;
        }
        std::vector<double> x, y;
        for (std::size_t i = begin; i < end; ++i)
        {
          x.push_back(points[i].first);
          y.push_back(points[i].second);
        }
        map.library_rt.push_back(median_(std::move(x)));
        map.run_rt.push_back(median_(std::move(y)));
      }
      return map;
    }
  }

  double RTMap::operator()(double x) const
  {
    if (library_rt.empty())
    {
      throw std::logic_error("Empty RT map");
    }
    if (x <= library_rt.front())
    {
      return run_rt.front();
    }
    if (x >= library_rt.back())
    {
      return run_rt.back();
    }
    // First knot strictly above x; the previous knot is at or below x.
    const std::size_t i = std::upper_bound(library_rt.begin(), library_rt.end(), x) - library_rt.begin();
    const double t = (x - library_rt[i - 1]) / (library_rt[i] - library_rt[i - 1]);
    return run_rt[i - 1] + t * (run_rt[i] - run_rt[i - 1]);
  }

  RTMap fitRTMap(std::vector<std::pair<double, double>> anchors)
  {
    if (anchors.size() < 20)
    {
      throw std::invalid_argument("RT calibration needs at least 20 anchors");
    }
    for (const auto& anchor : anchors)
    {
      if (!std::isfinite(anchor.first) || !std::isfinite(anchor.second))
      {
        throw std::invalid_argument("RT anchors must be finite");
      }
    }
    std::sort(anchors.begin(), anchors.end());
    RTMap map;
    // One fit plus up to three rejection refits.
    for (int round = 0; round < 4; ++round)
    {
      map = binnedKnots_(anchors);
      std::vector<double> residuals;
      residuals.reserve(anchors.size());
      for (const auto& anchor : anchors)
      {
        residuals.push_back(std::fabs(anchor.second - map(anchor.first)));
      }
      map.sigma = 1.4826 * median_(residuals);
      map.anchors = anchors.size();
      const double cut = std::max(3.0 * map.sigma, 3.0);
      std::vector<std::pair<double, double>> kept;
      for (std::size_t i = 0; i < anchors.size(); ++i)
      {
        if (residuals[i] <= cut)
        {
          kept.push_back(anchors[i]);
        }
      }
      if (kept.size() == anchors.size())
      {
        break;
      }
      if (kept.size() < 20)
      {
        throw std::invalid_argument("RT calibration rejected too many anchors");
      }
      anchors = std::move(kept);
    }
    return map;
  }

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
