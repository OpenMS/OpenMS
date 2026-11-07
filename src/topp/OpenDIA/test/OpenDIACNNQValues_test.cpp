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
    std::cout << "OpenDIACNNQValues tests passed\n";
    return 0;
  }
  catch (const std::exception& error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
