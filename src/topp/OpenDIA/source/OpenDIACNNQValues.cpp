// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenDIACNNQValues.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>

namespace OpenMS::OpenDIA
{
  namespace
  {
    /// Assign tie-grouped (D+1)*ratio/T q-values to competing indices, then take reverse cumulative minima.
    void assignFoldQValues_(const std::vector<CNNQValueInput>& inputs, std::vector<CNNQValueResult>& results,
                            std::vector<std::size_t>& indices, double target_decoy_ratio)
    {
      std::sort(indices.begin(), indices.end(), [&](std::size_t a, std::size_t b)
      {
        if (inputs[a].score != inputs[b].score)
        {
          return inputs[a].score > inputs[b].score;
        }
        return inputs[a].sample_id < inputs[b].sample_id;
      });

      std::size_t targets = 0;
      std::size_t decoys = 0;
      for (std::size_t begin = 0; begin < indices.size();)
      {
        std::size_t end = begin;
        // Include the whole exact-score tie before estimating its FDR.
        while (end < indices.size() && inputs[indices[end]].score == inputs[indices[begin]].score)
        {
          if (inputs[indices[end]].is_decoy)
          {
            ++decoys;
          }
          else
          {
            ++targets;
          }
          ++end;
        }
        // Multiply after dividing to avoid needless overflow for large ratios.
        const double fdr = targets == 0 ? 1.0 : std::min(1.0,
          ((static_cast<double>(decoys) + 1.0) / static_cast<double>(targets)) * target_decoy_ratio);
        for (std::size_t j = begin; j < end; ++j)
        {
          results[indices[j]].q_value = fdr;
        }
        begin = end;
      }
      double minimum = 1.0;
      for (auto it = indices.rbegin(); it != indices.rend(); ++it)
      {
        minimum = std::min(minimum, results[*it].q_value);
        results[*it].q_value = minimum;
      }
    }

    /// Validate inputs and select one winner per precursor; marks winners in results.
    std::map<std::string, std::size_t> selectWinners_(const std::vector<CNNQValueInput>& inputs,
                                                      std::vector<CNNQValueResult>& results)
    {
      results.reserve(inputs.size());
      std::set<std::string> sample_ids;
      std::map<std::string, std::size_t> winners;
      std::map<std::string, int> family_folds;
      for (std::size_t i = 0; i < inputs.size(); ++i)
      {
        const auto& row = inputs[i];
        for (const auto* id : {&row.sample_id, &row.precursor_id, &row.group_id})
        {
          if (id->find_first_not_of(" \t\r\n\v\f") == std::string::npos ||
              id->find_first_of("\t\r\n") != std::string::npos || id->find('\0') != std::string::npos)
          {
            throw std::invalid_argument("IDs must be nonblank and contain no tabs, newlines or NUL bytes");
          }
        }
        if (!sample_ids.insert(row.sample_id).second)
        {
          throw std::invalid_argument("Sample IDs must be unique");
        }
        if (!std::isfinite(row.score) || row.fold < 0)
        {
          throw std::invalid_argument("Scores must be finite and folds nonnegative");
        }
        const auto [family, new_family] = family_folds.emplace(row.group_id, row.fold);
        if (!new_family && family->second != row.fold)
        {
          throw std::invalid_argument("A peptide family must belong to one fold");
        }
        const auto [winner, new_precursor] = winners.emplace(row.precursor_id, i);
        if (!new_precursor)
        {
          const auto& previous = inputs[winner->second];
          if (row.is_decoy != previous.is_decoy || row.group_id != previous.group_id || row.fold != previous.fold ||
              row.pair_key != previous.pair_key)
          {
            throw std::invalid_argument("Precursor label, group, fold and pair key must be consistent");
          }
          if (row.score > previous.score || (row.score == previous.score && row.sample_id < previous.sample_id))
          {
            winner->second = i;
          }
        }
        results.push_back({row, false, 1.0});
      }
      for (const auto& entry : winners)
      {
        results[entry.second].precursor_best = true;
      }
      return winners;
    }
  }

  std::vector<CNNQValueResult> aggregateCNNQValues(
    const std::vector<CNNQValueInput>& inputs, double target_decoy_ratio)
  {
    if (!std::isfinite(target_decoy_ratio) || target_decoy_ratio <= 0.0)
    {
      throw std::invalid_argument("Target-to-decoy ratio must be positive and finite");
    }

    std::vector<CNNQValueResult> results;
    const auto winners = selectWinners_(inputs, results);

    // Only precursor winners contribute to fold-local assay counts.
    std::map<int, std::vector<std::size_t>> folds;
    for (const auto& entry : winners)
    {
      folds[inputs[entry.second].fold].push_back(entry.second);
    }
    for (auto& entry : folds)
    {
      auto& indices = entry.second;
      const auto decoy_count = std::count_if(indices.begin(), indices.end(), [&](std::size_t i)
      {
        return inputs[i].is_decoy;
      });
      if (decoy_count == 0 || static_cast<std::size_t>(decoy_count) == indices.size())
      {
        throw std::invalid_argument("Each represented fold needs target and decoy precursor winners");
      }
      assignFoldQValues_(inputs, results, indices, target_decoy_ratio);
    }
    return results;
  }

  std::vector<CNNQValueResult> aggregatePairedCNNQValues(const std::vector<CNNQValueInput>& inputs)
  {
    std::vector<CNNQValueResult> results;
    const auto winners = selectWinners_(inputs, results);

    // Group winners by pair key; each pair holds at most one target and one decoy in one fold.
    struct Pair
    {
      std::size_t target = SIZE_MAX;
      std::size_t decoy = SIZE_MAX;
    };
    std::map<std::string, Pair> pairs;
    std::map<int, std::vector<std::size_t>> folds;
    for (const auto& entry : winners)
    {
      const std::size_t index = entry.second;
      const auto& row = inputs[index];
      if (row.pair_key.empty())
      {
        folds[row.fold].push_back(index);
        continue;
      }
      auto& pair = pairs[row.pair_key];
      auto& slot = row.is_decoy ? pair.decoy : pair.target;
      if (slot != SIZE_MAX)
      {
        throw std::invalid_argument("A pair key may join at most one target and one decoy precursor");
      }
      slot = index;
    }
    for (const auto& entry : pairs)
    {
      const auto& pair = entry.second;
      if (pair.target == SIZE_MAX || pair.decoy == SIZE_MAX)
      {
        const std::size_t single = pair.target == SIZE_MAX ? pair.decoy : pair.target;
        folds[inputs[single].fold].push_back(single);
        continue;
      }
      if (inputs[pair.target].fold != inputs[pair.decoy].fold)
      {
        throw std::invalid_argument("Paired target and decoy precursors must share a fold");
      }
      // Exact ties go to the decoy, the conservative choice.
      const std::size_t competitor = inputs[pair.target].score > inputs[pair.decoy].score ? pair.target : pair.decoy;
      folds[inputs[competitor].fold].push_back(competitor);
    }
    for (auto& entry : folds)
    {
      assignFoldQValues_(inputs, results, entry.second, 1.0);
    }
    return results;
  }
}
