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
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cstdlib>
#include <numeric>

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

START_SECTION([EXTRA] reproducibility with fixed seed)
{
  // Same input + same seed on two instances → identical scores.
  RescoreInput input;
  input.feature_names = StringList{"f"};
  const size_t n = 400;
  input.features.reserve(n);
  input.is_decoy.reserve(n);
  for (size_t i = 0; i < n; ++i)
  {
    const bool is_decoy = (i % 2 == 1);
    input.features.push_back({(is_decoy ? -1.0 : +1.0) + 0.01 * static_cast<double>(i)});
    input.is_decoy.push_back(is_decoy);
  }

  Percolator p1;
  Percolator p2;
  for (auto* p : {&p1, &p2})
  {
    Param par = p->getDefaults();
    par.setValue("seed", 17);
    p->setParameters(par);
  }

  RescoreOutput out1 = p1.rescore(input);
  RescoreOutput out2 = p2.rescore(input);

  TEST_EQUAL(out1.scores.size(), out2.scores.size())
  bool all_equal = out1.scores.size() == out2.scores.size();
  for (size_t i = 0; all_equal && i < out1.scores.size(); ++i)
  {
    if (std::abs(out1.scores[i] - out2.scores[i]) > 1e-9) all_equal = false;
  }
  TEST_EQUAL(all_equal, true)
}
END_SECTION

START_SECTION([EXTRA] q-values in expected range)
{
  // Moderately separable data. Verify q-values are bounded [0, 1] and that
  // targets scoring high have low q-values on average.
  //
  // Deliberately NOT checking strict monotonicity across score-sorted rows:
  // Percolator computes q-values per CV fold and merges by score — the
  // merged sequence can have small out-of-order runs where per-fold q-value
  // estimates differ near tied scores. That's an upstream merge quirk,
  // not a property the wrapper is obliged to enforce.
  RescoreInput input;
  input.feature_names = StringList{"f", "noise"};
  std::srand(42);
  auto rand01 = []() { return static_cast<double>(std::rand()) / RAND_MAX; };
  const size_t n_per_class = 500;
  for (size_t i = 0; i < n_per_class; ++i)
  {
    const double t = 1.0 + 0.6 * (rand01() - 0.5) * 2.0;
    const double d = 0.0 + 0.6 * (rand01() - 0.5) * 2.0;
    input.features.push_back({t, rand01()});
    input.is_decoy.push_back(false);
    input.features.push_back({d, rand01()});
    input.is_decoy.push_back(true);
  }

  Percolator p;
  RescoreOutput out = p.rescore(input);

  // Every q-value must be in [0, 1].
  for (double q : out.q_values) TEST_TRUE(q >= 0.0 && q <= 1.0 + 1e-9)
  // Every PEP must be in [0, 1].
  for (double pep : out.peps)   TEST_TRUE(pep >= 0.0 && pep <= 1.0 + 1e-9)

  // Top-scoring targets should have low mean q-value; bottom-scoring decoys
  // should have high mean q-value.
  std::vector<size_t> idx(out.scores.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
    [&](size_t a, size_t b){ return out.scores[a] > out.scores[b]; });

  double mean_q_top = 0, mean_q_bot = 0;
  const size_t window = 50;
  for (size_t k = 0; k < window; ++k) mean_q_top += out.q_values[idx[k]];
  for (size_t k = out.q_values.size() - window; k < out.q_values.size(); ++k)
    mean_q_bot += out.q_values[idx[k]];
  mean_q_top /= window;
  mean_q_bot /= window;
  TEST_EQUAL(mean_q_top < mean_q_bot, true)
}
END_SECTION

START_SECTION((void rescore(std::vector<PeptideIdentification>& peptide_ids, const StringList& feature_names)))
{
  // Build a set of PeptideIdentifications with target/decoy meta values and
  // two numeric features; run rescore; verify percolator_score / _q_value /
  // _pep meta values landed on each hit and the target mean score > decoy mean.
  std::vector<PeptideIdentification> peps;
  std::srand(7);
  auto rand01 = []() { return static_cast<double>(std::rand()) / RAND_MAX; };

  const size_t n = 800;
  for (size_t i = 0; i < n; ++i)
  {
    const bool is_decoy = (i % 2 == 1);
    PeptideIdentification pid;
    pid.setRT(static_cast<double>(i) * 0.1);
    pid.setIdentifier("run1");
    PeptideHit hit;
    hit.setScore(0.0);
    hit.setMetaValue("target_decoy", is_decoy ? "decoy" : "target");
    const double f = (is_decoy ? 0.0 : 1.0) + 0.6 * (rand01() - 0.5) * 2.0;
    hit.setMetaValue("feat_sep", f);
    hit.setMetaValue("feat_noise", rand01());
    pid.insertHit(hit);
    peps.push_back(pid);
  }

  Percolator p;
  p.rescore(peps, StringList{"feat_sep", "feat_noise"});

  // Every hit has the three new meta values.
  double target_sum = 0, decoy_sum = 0;
  size_t n_targets = 0, n_decoys = 0;
  for (const auto& pid : peps)
  {
    TEST_EQUAL(pid.getHits().size(), 1)
    const auto& hit = pid.getHits()[0];
    TEST_TRUE(hit.metaValueExists("percolator_score"))
    TEST_TRUE(hit.metaValueExists("percolator_q_value"))
    TEST_TRUE(hit.metaValueExists("percolator_pep"))
    const double s = hit.getMetaValue("percolator_score");
    if (hit.getMetaValue("target_decoy").toString() == "decoy")
    { decoy_sum += s; ++n_decoys; }
    else
    { target_sum += s; ++n_targets; }
  }
  TEST_EQUAL(n_targets > 0, true)
  TEST_EQUAL(n_decoys > 0, true)
  TEST_EQUAL(target_sum / n_targets > decoy_sum / n_decoys, true)
}
END_SECTION

START_SECTION([EXTRA] missing target_decoy meta throws)
{
  std::vector<PeptideIdentification> peps;
  for (size_t i = 0; i < 20; ++i)
  {
    PeptideIdentification pid;
    PeptideHit h;
    h.setMetaValue("feat_x", 1.0);
    pid.insertHit(h);
    peps.push_back(pid);
  }
  Percolator perc;
  TEST_EXCEPTION(Exception::InvalidValue, perc.rescore(peps, StringList{"feat_x"}))
}
END_SECTION

START_SECTION([EXTRA] malformed input throws InvalidValue)
{
  // Empty input
  {
    Percolator p;
    RescoreInput input;
    TEST_EXCEPTION(Exception::InvalidValue, p.rescore(input))
  }
  // Length mismatch between features and is_decoy
  {
    Percolator p;
    RescoreInput input;
    input.features.push_back({1.0});
    input.features.push_back({1.0});
    input.is_decoy.push_back(false);
    TEST_EXCEPTION(Exception::InvalidValue, p.rescore(input))
  }
  // Uneven feature row lengths
  {
    Percolator p;
    RescoreInput input;
    input.features.push_back({1.0, 2.0});
    input.features.push_back({1.0});
    input.is_decoy.push_back(false);
    input.is_decoy.push_back(true);
    TEST_EXCEPTION(Exception::InvalidValue, p.rescore(input))
  }
}
END_SECTION

END_TEST
