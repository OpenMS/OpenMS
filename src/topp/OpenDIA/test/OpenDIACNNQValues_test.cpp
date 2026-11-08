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
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>

using namespace OpenMS::OpenDIA;

namespace
{
  void check(bool condition)
  {
    if (!condition)
    {
      throw std::runtime_error("CNN q-value test failed");
    }
  }

  void near(double actual, double expected)
  {
    check(std::abs(actual - expected) < 1e-12);
  }

  void invalid(const std::vector<CNNQValueInput>& rows, double ratio = 1.0)
  {
    try
    {
      aggregateCNNQValues(rows, ratio);
    }
    catch (const std::invalid_argument&)
    {
      return;
    }
    throw std::runtime_error("Expected invalid_argument");
  }

  void paired_invalid_pooled_(const std::vector<CNNQValueInput>& rows)
  {
    try
    {
      aggregatePairedCNNQValues(rows, true);
    }
    catch (const std::invalid_argument&)
    {
      return;
    }
    throw std::runtime_error("Expected pooled paired invalid_argument");
  }

  std::map<std::string, std::pair<bool, double>> keyed(const std::vector<CNNQValueResult>& rows)
  {
    std::map<std::string, std::pair<bool, double>> result;
    for (const auto& row : rows)
    {
      result.emplace(row.input.sample_id, std::make_pair(row.precursor_best, row.q_value));
    }
    return result;
  }
}

int main()
{
  try
  {
    // The two target winners at score 10 must both get 1/2, not 1.
    std::vector<CNNQValueInput> rows{
      {"z", "p1", "g1", false, 0, 10},
      {"a", "p1", "g1", false, 0, 10},
      {"low", "p1", "g1", false, 0, 1},
      {"t2", "p2", "g1", false, 0, 10},
      {"d1", "d1", "gd1", true, 0, 9},
      {"t3", "p3", "g3", false, 0, 8},
      {"d2", "d2", "gd2", true, 0, 8},
      {"other_t", "other_p", "other_g", false, 7, -100},
      {"other_d", "other_d", "other_gd", true, 7, -200}
    };
    const auto result = aggregateCNNQValues(rows);
    check(result.size() == rows.size());
    for (std::size_t i = 0; i < rows.size(); ++i)
    {
      check(result[i].input.sample_id == rows[i].sample_id);
      check(result[i].input.precursor_id == rows[i].precursor_id);
      check(result[i].input.group_id == rows[i].group_id);
      check(result[i].input.is_decoy == rows[i].is_decoy);
      check(result[i].input.fold == rows[i].fold);
      check(result[i].input.score == rows[i].score);
    }
    check(!result[0].precursor_best && result[1].precursor_best && !result[2].precursor_best);
    near(result[0].q_value, 1.0);
    near(result[2].q_value, 1.0);
    near(result[1].q_value, 0.5);
    near(result[3].q_value, 0.5);
    near(result[4].q_value, 1.0);
    near(result[5].q_value, 1.0);
    near(result[6].q_value, 1.0);
    near(result[7].q_value, 1.0);
    near(aggregateCNNQValues(rows, 0.5)[1].q_value, 0.25);

    // Independent positive affine model scales cannot change fold-local output.
    auto transformed = rows;
    for (auto& row : transformed)
    {
      row.score = row.fold == 0 ? row.score * 3 + 500 : row.score * 2 + 10000;
    }
    check(keyed(aggregateCNNQValues(transformed)) == keyed(result));
    std::reverse(rows.begin(), rows.end());
    check(keyed(aggregateCNNQValues(rows)) == keyed(result));
    std::reverse(rows.begin(), rows.end());
    for (const auto& a : result)
    {
      for (const auto& b : result)
      {
        if (a.precursor_best && b.precursor_best && a.input.fold == b.input.fold && a.input.score >= b.input.score)
        {
          check(a.q_value <= b.q_value);
        }
      }
    }

    // Reverse minima can lower even a leading decoy's initial T=0 FDR of 1.
    const std::vector<CNNQValueInput> leading_decoy{
      {"d", "d", "d", true, 0, 10},
      {"t1", "t1", "t1", false, 0, 9},
      {"t2", "t2", "t2", false, 0, 8},
      {"t3", "t3", "t3", false, 0, 7}
    };
    const auto reverse_min = aggregateCNNQValues(leading_decoy);
    for (const auto& row : reverse_min)
    {
      near(row.q_value, 2.0 / 3.0);
    }
    auto all_tied = leading_decoy;
    for (auto& row : all_tied)
    {
      row.score = 0;
    }
    for (const auto& row : aggregateCNNQValues(all_tied))
    {
      near(row.q_value, 2.0 / 3.0);
    }

    check(aggregateCNNQValues({}).empty());
    for (double ratio : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
      invalid(rows, ratio);
      invalid({}, ratio);
    }
    auto bad = rows;
    bad[0].sample_id = bad[1].sample_id;
    invalid(bad);
    for (double score : {std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()})
    {
      bad = rows;
      bad[0].score = score;
      invalid(bad);
    }
    for (const char separator : {'\t', '\n', '\r'})
    {
      for (int field = 0; field < 3; ++field)
      {
        bad = rows;
        auto& row = bad[0];
        (field == 0 ? row.sample_id : field == 1 ? row.precursor_id : row.group_id) += separator;
        invalid(bad);
      }
    }
    for (const auto& invalid_id : {std::string{}, std::string("   "), std::string("\v\f"), std::string("id\0suffix", 9)})
    {
      for (int field = 0; field < 3; ++field)
      {
        bad = rows;
        auto& row = bad[0];
        (field == 0 ? row.sample_id : field == 1 ? row.precursor_id : row.group_id) = invalid_id;
        invalid(bad);
      }
    }
    bad = rows;
    bad[0].fold = -1;
    invalid(bad);
    bad = rows;
    bad[0].is_decoy = true;
    invalid(bad);
    bad = rows;
    bad[0].group_id = "inconsistent";
    invalid(bad);
    bad = rows;
    bad[0].fold = 1;
    invalid(bad);
    bad = rows;
    bad[3].fold = 1; // Different precursors in the same family must share a fold.
    invalid(bad);
    invalid({rows[0]});
    invalid({rows[4]});
    bad = rows;
    bad.pop_back(); // Fold 7 lacks decoys, even though fold 0 has both classes.
    invalid(bad);
    bad = rows;
    bad.erase(bad.begin() + 7); // Fold 7 lacks targets.
    invalid(bad);

    // Paired competition: pair "a" target wins, pair "b" decoy wins (tie goes to decoy),
    // pair "c" target wins, unpaired target "u" competes alone.
    // Two target windows for precursor ta check that only the best window competes.
    const std::vector<CNNQValueInput> paired_rows = {
      {"ta1", "ta", "fa", false, 0, 9.0, "a"}, {"ta2", "ta", "fa", false, 0, 1.0, "a"},
      {"da", "dA", "fa", true, 0, 8.0, "a"},
      {"tb", "tb", "fb", false, 0, 5.0, "b"}, {"db", "dB", "fb", true, 0, 5.0, "b"},
      {"tc", "tc", "fc", false, 0, 4.0, "c"}, {"dc", "dC", "fc", true, 0, 0.5, "c"},
      {"tu", "tu", "fu", false, 0, 3.0, ""}};
    const auto paired = keyed(aggregatePairedCNNQValues(paired_rows));
    // Competitors by score: ta 9 (T), db 5 (D), tc 4 (T), tu 3 (T).
    // FDR: ta (0+1)/1 = 1, db (1+1)/1 -> 1, tc (1+1)/2 = 1, tu (1+1)/3 = 2/3.
    near(paired.at("ta1").second, 2.0 / 3.0);
    check(paired.at("ta1").first && !paired.at("ta2").first);
    near(paired.at("ta2").second, 1.0);
    near(paired.at("tb").second, 1.0); // Lost its tied pair.
    near(paired.at("da").second, 1.0); // Lost to its target.
    near(paired.at("tc").second, 2.0 / 3.0);
    near(paired.at("tu").second, 2.0 / 3.0);
    // Input order does not change paired results.
    auto reversed_rows = paired_rows;
    std::reverse(reversed_rows.begin(), reversed_rows.end());
    check(keyed(aggregatePairedCNNQValues(reversed_rows)) == paired);

    // Pooled paired competition: fold 1 logits are fold 0 logits x 10 + 100, so after per-fold decoy
    // standardisation both folds coincide and one pooled list doubles every count.
    {
      std::vector<CNNQValueInput> pooled_rows;
      const std::vector<double> target_scores = {9.0, 7.0, 6.0, 1.0}, decoy_scores = {8.0, 2.0, 0.0, -1.0};
      for (int fold = 0; fold < 2; ++fold)
      {
        for (std::size_t k = 0; k < target_scores.size(); ++k)
        {
          const std::string key = std::to_string(fold) + "p" + std::to_string(k);
          const auto shift = [&](double s) { return fold == 0 ? s : 10.0 * s + 100.0; };
          pooled_rows.push_back({"t" + key, "t" + key, "f" + key, false, fold, shift(target_scores[k]), key});
          pooled_rows.push_back({"d" + key, "d" + key, "f" + key, true, fold, shift(decoy_scores[k]), key});
        }
      }
      const auto per_fold = keyed(aggregatePairedCNNQValues(pooled_rows));
      const auto pooled = keyed(aggregatePairedCNNQValues(pooled_rows, true));
      // Per fold: competitors 9 T, 7 T, 6 T, 1 T (all targets win) -> q = 1/4 at the bottom, (0+1)/1 at the top.
      near(per_fold.at("t0p0").second, 0.25);
      near(per_fold.at("t1p0").second, 0.25);
      // Pooled: 8 target competitors, no decoy winners -> (0+1)/8.
      near(pooled.at("t0p0").second, 0.125);
      near(pooled.at("t1p3").second, 0.125);
      // Winner selection and reported scores are unchanged.
      for (const auto& [id, value] : per_fold)
      {
        check(pooled.at(id).first == value.first);
      }
      const auto raw = aggregatePairedCNNQValues(pooled_rows, true);
      for (std::size_t i = 0; i < raw.size(); ++i)
      {
        near(raw[i].input.score, pooled_rows[i].score);
      }
      // A fold without decoy winners cannot be standardised.
      paired_invalid_pooled_({{"t1", "t1", "f1", false, 0, 1.0, "k1"}, {"d1", "d1", "f1", true, 0, 0.0, "k1"},
                              {"t2", "t2", "f2", false, 1, 1.0, ""}});
    }
    // Pairing rejects duplicate classes in one key and split folds.
    const auto paired_invalid = [](const std::vector<CNNQValueInput>& rows)
    {
      try
      {
        aggregatePairedCNNQValues(rows);
      }
      catch (const std::invalid_argument&)
      {
        return;
      }
      throw std::runtime_error("Expected paired invalid_argument");
    };
    paired_invalid({{"t1", "t1", "f1", false, 0, 1.0, "k"}, {"t2", "t2", "f2", false, 0, 2.0, "k"}});
    paired_invalid({{"t1", "t1", "f1", false, 0, 1.0, "k"}, {"d1", "d1", "f2", true, 1, 2.0, "k"}});
    paired_invalid({{"t1", "t1", "f1", false, 0, 1.0, "k"}, {"t1b", "t1", "f1", false, 0, 2.0, "other"}});
    // Without pair keys, every winner competes alone like the unpaired estimator at ratio 1.
    auto unkeyed = paired_rows;
    for (auto& row : unkeyed)
    {
      row.pair_key.clear();
    }
    check(keyed(aggregatePairedCNNQValues(unkeyed)) == keyed(aggregateCNNQValues(unkeyed, 1.0)));
    std::cout << "OpenDIACNNQValues tests passed\n";
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
