// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/ID/PercolatorTypes.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

START_TEST(Percolator, "$Id$")

START_SECTION((RescoreOutput rescore(const RescoreInput& input)))
{
  // Trivially-separable two-class problem: target rows have feature[0] = +1,
  // decoy rows have feature[0] = -1.
  RescoreInput input;
  const size_t n_per_class = 200;
  input.feature_names = StringList{"separator", "noise"};
  input.features.reserve(2 * n_per_class);
  input.is_decoy.reserve(2 * n_per_class);

  for (size_t i = 0; i < n_per_class; ++i)
  {
    input.features.push_back({+1.0 + 0.01 * static_cast<double>(i),
                              static_cast<double>(i % 7)});
    input.is_decoy.push_back(false);
  }
  for (size_t i = 0; i < n_per_class; ++i)
  {
    input.features.push_back({-1.0 + 0.01 * static_cast<double>(i),
                              static_cast<double>(i % 11)});
    input.is_decoy.push_back(true);
  }

  Percolator p;
  Param param = p.getDefaults();
  param.setValue("seed", 42);
  p.setParameters(param);

  RescoreOutput out = p.rescore(input);

  TEST_EQUAL(out.scores.size(), input.features.size())
  TEST_EQUAL(out.q_values.size(), input.features.size())
  TEST_EQUAL(out.peps.size(), input.features.size())

  // Median target score should exceed median decoy score.
  vector<double> target_scores, decoy_scores;
  for (size_t i = 0; i < out.scores.size(); ++i)
  {
    (input.is_decoy[i] ? decoy_scores : target_scores).push_back(out.scores[i]);
  }
  std::sort(target_scores.begin(), target_scores.end());
  std::sort(decoy_scores.begin(), decoy_scores.end());
  const double target_median = target_scores[target_scores.size() / 2];
  const double decoy_median  = decoy_scores[decoy_scores.size() / 2];
  TEST_EQUAL(target_median > decoy_median, true)
}
END_SECTION

END_TEST
