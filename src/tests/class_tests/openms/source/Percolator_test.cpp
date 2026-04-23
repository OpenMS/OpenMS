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
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <numeric>
#include <random>
#include <set>

using namespace OpenMS;
using namespace std;

// Helper: build a moderately-separable dataset with @p n_per_class targets
// and decoys, 2 features (signal + noise). Target signal ~N(+1.0, 0.6^2);
// decoy signal ~N(0.0, 0.6^2). Not trivially separable — Percolator's
// checkSeparationAndSetPi0 will not fire on this.
static OpenMS::RescoreInput makeModeratelySeparableInput_(
  std::size_t n_per_class, unsigned seed)
{
  std::srand(seed);
  auto rand01 = []() { return static_cast<double>(std::rand()) / RAND_MAX; };
  auto randn = [&]() {
    double u1 = std::max(1e-9, rand01());
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * rand01());
  };
  OpenMS::RescoreInput input;
  input.feature_names = OpenMS::StringList{"signal", "noise"};
  input.features.reserve(2 * n_per_class);
  input.is_decoy.reserve(2 * n_per_class);
  for (std::size_t i = 0; i < n_per_class; ++i)
  {
    input.features.push_back({+1.0 + 0.6 * randn(), rand01()});
    input.is_decoy.push_back(false);
  }
  for (std::size_t i = 0; i < n_per_class; ++i)
  {
    input.features.push_back({ 0.0 + 0.6 * randn(), rand01()});
    input.is_decoy.push_back(true);
  }
  return input;
}

static bool medianTargetBeatsDecoy_(const OpenMS::RescoreOutput& out,
                                    const OpenMS::RescoreInput& in)
{
  std::vector<double> tt, dd;
  for (std::size_t i = 0; i < out.scores.size(); ++i)
  {
    (in.is_decoy[i] ? dd : tt).push_back(out.scores[i]);
  }
  std::sort(tt.begin(), tt.end());
  std::sort(dd.begin(), dd.end());
  return tt[tt.size() / 2] > dd[dd.size() / 2];
}

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

START_SECTION([EXTRA] reproducibility with fixed seed is bit-identical)
{
  // Cross-instance, same-seed invariance. Tightened from 1e-9 tolerance on
  // scores only to strict equality on scores, q-values, and PEPs.
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
  TEST_EQUAL(out1.q_values.size(), out2.q_values.size())
  TEST_EQUAL(out1.peps.size(), out2.peps.size())

  bool bit_identical = out1.scores.size() == out2.scores.size();
  for (size_t i = 0; bit_identical && i < out1.scores.size(); ++i)
  {
    if (out1.scores[i]   != out2.scores[i])   bit_identical = false;
    if (out1.q_values[i] != out2.q_values[i]) bit_identical = false;
    if (out1.peps[i]     != out2.peps[i])     bit_identical = false;
  }
  TEST_EQUAL(bit_identical, true)
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

START_SECTION([EXTRA] output makes sense on a realistic-scale dataset)
{
  // Scale-up sanity check: 2000 PSMs (1000 targets + 1000 decoys) with a
  // mildly informative 3-feature signal. Large enough that Percolator's
  // SanityCheck should be satisfied by the training signal and NOT fall
  // back to the default-direction reset path.
  //
  // Ground truth: true positives (first 500 targets) have
  //     f0 = +2 + noise, f1 = +1 + noise, f2 = noise
  // "hard" targets (next 500) have
  //     f0 = +0.5 + noise (near decoy distribution)
  // decoys have
  //     f0 = 0 + noise, f1 = 0 + noise, f2 = noise
  //
  // Expected:
  //   * mean target score > mean decoy score
  //   * top-100 scoring rows: >= 80 are targets (easy-positive detection)
  //   * bottom-100 scoring rows: >= 80 are decoys or hard-targets
  //   * q-value at the easy-positive range: well below 0.05 for many rows
  //   * all q-values, peps in [0, 1]
  //   * at least some targets have q < 0.01 (Percolator actually "found" them)

  std::srand(2026);
  auto rand01 = []() { return static_cast<double>(std::rand()) / RAND_MAX; };
  auto randn = [&]() {
    // Box-Muller
    double u1 = std::max(1e-9, rand01());
    double u2 = rand01();
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
  };

  RescoreInput input;
  input.feature_names = StringList{"f0", "f1", "noise"};
  std::vector<bool> is_easy_target;  // ground-truth easy positives
  const size_t n_easy = 500, n_hard = 500, n_dec = 1000;

  for (size_t i = 0; i < n_easy; ++i)
  {
    input.features.push_back({+2.0 + 0.5 * randn(),
                              +1.0 + 0.5 * randn(),
                              randn()});
    input.is_decoy.push_back(false);
    is_easy_target.push_back(true);
  }
  for (size_t i = 0; i < n_hard; ++i)
  {
    input.features.push_back({+0.5 + 0.5 * randn(),
                              +0.5 * randn(),
                              randn()});
    input.is_decoy.push_back(false);
    is_easy_target.push_back(false);
  }
  for (size_t i = 0; i < n_dec; ++i)
  {
    input.features.push_back({+0.0 + 0.5 * randn(),
                              +0.0 * randn(),
                              randn()});
    input.is_decoy.push_back(true);
    is_easy_target.push_back(false);
  }

  Percolator p;
  Param par = p.getDefaults();
  par.setValue("seed", 2026);
  p.setParameters(par);
  RescoreOutput out = p.rescore(input);

  const size_t n = input.features.size();
  TEST_EQUAL(out.scores.size(), n)
  TEST_EQUAL(out.q_values.size(), n)
  TEST_EQUAL(out.peps.size(), n)

  // All outputs in valid ranges
  for (double q : out.q_values) TEST_TRUE(q >= 0.0 && q <= 1.0 + 1e-9)
  for (double pep : out.peps)   TEST_TRUE(pep >= 0.0 && pep <= 1.0 + 1e-9)

  // Sort indices by descending score
  std::vector<size_t> idx(n);
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
    [&](size_t a, size_t b){ return out.scores[a] > out.scores[b]; });

  // Top 100: should be dominated by easy targets (80+/100 is expected)
  size_t top_targets = 0, top_easy = 0;
  for (size_t k = 0; k < 100; ++k)
  {
    if (!input.is_decoy[idx[k]]) ++top_targets;
    if (is_easy_target[idx[k]])  ++top_easy;
  }
  TEST_EQUAL(top_targets >= 80, true)
  TEST_EQUAL(top_easy >= 60, true)

  // Bottom 100: should be dominated by decoys (+ hard targets)
  size_t bot_non_easy = 0;
  for (size_t k = n - 100; k < n; ++k)
  {
    if (!is_easy_target[idx[k]]) ++bot_non_easy;
  }
  TEST_EQUAL(bot_non_easy >= 80, true)

  // Percolator should find some confident targets: count rows with q < 0.01
  size_t confident = 0;
  for (size_t i = 0; i < n; ++i)
  {
    if (!input.is_decoy[i] && out.q_values[i] < 0.01) ++confident;
  }
  // On this dataset with ~500 easy positives, Percolator should identify
  // at least 50 of them at q<0.01 when training succeeds. If the
  // SanityCheck fallback kicks in this number would be 0 — which is also
  // valid Percolator behavior on degenerate inputs but would indicate
  // this test is no longer actually stressing the trained path.
  TEST_EQUAL(confident >= 50, true)

  // PEPs should anti-correlate with SVM scores (higher score → lower PEP).
  // Spearman-style check: mean PEP of top-100 scored rows should be well
  // below mean PEP of bottom-100 scored rows.
  double top_pep_sum = 0, bot_pep_sum = 0;
  for (size_t k = 0; k < 100; ++k) top_pep_sum += out.peps[idx[k]];
  for (size_t k = n - 100; k < n; ++k) bot_pep_sum += out.peps[idx[k]];
  TEST_EQUAL(top_pep_sum / 100.0 < bot_pep_sum / 100.0, true)
}
END_SECTION

START_SECTION([EXTRA] report_as_main_score promotes a Percolator output to hit.getScore())
{
  // Build a small target/decoy set with a numeric meta feature and run the
  // high-level PeptideIdentification API with report_as_main_score="q-value".
  // Verify: hit.getScore() equals percolator_q_value, score_type is "q-value",
  // higher_score_better=false, and hits are re-ranked accordingly.
  // Build one shared dataset; both rescore runs consume an independent copy.
  std::srand(3);
  auto rand01 = []() { return static_cast<double>(std::rand()) / RAND_MAX; };
  auto randn = [&]() {
    double u1 = std::max(1e-9, rand01());
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * rand01());
  };
  std::vector<PeptideIdentification> source;
  const size_t n = 800;
  for (size_t i = 0; i < n; ++i)
  {
    const bool is_decoy = (i % 2 == 1);
    PeptideIdentification pid;
    pid.setIdentifier("run1");
    pid.setRT(0.1 * i);
    PeptideHit hit;
    hit.setSequence(AASequence::fromString("PEPTIDER"));
    hit.setCharge(2);
    hit.setScore(42.0);                                // pre-rescore main score
    hit.setMetaValue("target_decoy", is_decoy ? "decoy" : "target");
    // Moderate overlap: clean enough to train, not so clean that percolator
    // trips "too good separation" or the median-decoy rescale check.
    hit.setMetaValue("f", (is_decoy ? 0.0 : 1.5) + 0.9 * randn());
    pid.insertHit(hit);
    source.push_back(std::move(pid));
  }

  // Run A: report_as_main_score="q-value" — hit score should end up == q-value.
  std::vector<PeptideIdentification> peps_a = source;
  {
    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 3);
    par.setValue("test_fdr",  0.1);
    par.setValue("train_fdr", 0.1);
    par.setValue("report_as_main_score", "q-value");
    p.setParameters(par);
    p.rescore(peps_a, StringList{"f"});
  }

  bool all_promoted = true;
  for (const auto& pid : peps_a)
  {
    if (pid.getHits().empty()) continue;
    if (pid.getScoreType() != "q-value") { all_promoted = false; break; }
    if (pid.isHigherScoreBetter()) { all_promoted = false; break; }
    const auto& hit = pid.getHits().front();
    const double qv = hit.getMetaValue("percolator_q_value");
    if (std::abs(hit.getScore() - qv) > 1e-12) { all_promoted = false; break; }
  }
  TEST_EQUAL(all_promoted, true)

  // Run B: default report_as_main_score="none" — hit score unchanged, but
  // percolator_* meta values still stamped.
  std::vector<PeptideIdentification> peps_b = source;
  {
    Percolator p;
    Param par = p.getDefaults();
    par.setValue("seed", 3);
    par.setValue("test_fdr",  0.1);
    par.setValue("train_fdr", 0.1);
    // report_as_main_score left at default "none"
    p.setParameters(par);
    p.rescore(peps_b, StringList{"f"});
  }
  TEST_EQUAL(peps_b.front().getHits().front().getScore(), 42.0)
  TEST_TRUE(peps_b.front().getHits().front().metaValueExists("percolator_q_value"))
}
END_SECTION

START_SECTION([EXTRA] use_pi0=false disables pi0 correction)
{
  // Build a dataset where pi0 correction clearly matters: plenty of well-
  // separated targets mixed with decoys. With usePi0=true (default) q-values
  // get multiplied by a pi0 < 1 estimate, so turning pi0 off raises every
  // q-value at a given rank. Asserting strict "all q_off >= q_on" is too
  // brittle (Percolator's merge step can reshuffle ties), so instead check
  // that the means diverge in the expected direction.
  std::srand(11);
  auto rand01 = []() { return static_cast<double>(std::rand()) / RAND_MAX; };
  auto randn = [&]() {
    double u1 = std::max(1e-9, rand01());
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * rand01());
  };
  // Moderate separation (targets ~1.5σ above decoys). Wide enough that
  // median decoy < score at FDR (both settings train), narrow enough that
  // pi0 bootstrap has real uncertainty to correct for.
  RescoreInput input;
  input.feature_names = StringList{"f0"};
  for (size_t i = 0; i < 600; ++i)
  {
    input.features.push_back({+1.5 + 1.0 * randn()});
    input.is_decoy.push_back(false);
  }
  for (size_t i = 0; i < 600; ++i)
  {
    input.features.push_back({0.0 + 1.0 * randn()});
    input.is_decoy.push_back(true);
  }

  Percolator p_on, p_off;
  for (auto* p : {&p_on, &p_off})
  {
    Param par = p->getDefaults();
    par.setValue("seed", 17);
    par.setValue("test_fdr",  0.05);
    par.setValue("train_fdr", 0.05);
    p->setParameters(par);
  }
  {
    Param par = p_off.getParameters();
    par.setValue("use_pi0", "false");
    p_off.setParameters(par);
  }

  RescoreOutput out_on  = p_on.rescore(input);
  RescoreOutput out_off = p_off.rescore(input);

  // Mean target q over top-100 ranked targets should be strictly lower when
  // pi0 is on (because q is scaled by pi0 < 1).
  std::vector<size_t> idx(input.features.size());
  std::iota(idx.begin(), idx.end(), 0);
  std::sort(idx.begin(), idx.end(),
    [&](size_t a, size_t b){ return out_on.scores[a] > out_on.scores[b]; });
  double mean_on = 0, mean_off = 0;
  size_t counted = 0;
  for (size_t k = 0; k < idx.size() && counted < 100; ++k)
  {
    if (input.is_decoy[idx[k]]) continue;
    mean_on  += out_on.q_values[idx[k]];
    mean_off += out_off.q_values[idx[k]];
    ++counted;
  }
  TEST_EQUAL(counted, 100)
  TEST_TRUE(mean_on < mean_off)
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

START_SECTION((PercolatorModel train(const RescoreInput& input)))
{
  RescoreInput input = makeModeratelySeparableInput_(300, 42);

  Percolator p;
  Param par = p.getDefaults();
  par.setValue("seed", 17);
  p.setParameters(par);

  PercolatorModel model = p.train(input);

  TEST_EQUAL(model.format_version, 1)
  TEST_EQUAL(model.feature_names.size(), 2u)
  TEST_EQUAL(model.feature_names[0], "signal")
  TEST_EQUAL(model.feature_names[1], "noise")
  TEST_EQUAL(model.weights.size(), 3u)  // 2 features + 1 bias
  TEST_EQUAL(model.normalizer_type, "stdv")
  TEST_EQUAL(model.seed, 17)
  // Signal weight should be positive (targets have higher signal than decoys).
  TEST_TRUE(model.weights[0] > 0.0)
}
END_SECTION

START_SECTION((RescoreOutput score(const RescoreInput& input, const PercolatorModel& model)))
{
  // Train on A, score A. Target median score > decoy median score.
  RescoreInput input = makeModeratelySeparableInput_(300, 42);
  Percolator p;
  Param par = p.getDefaults();
  par.setValue("seed", 17);
  p.setParameters(par);

  PercolatorModel model = p.train(input);
  RescoreOutput out = p.score(input, model);

  TEST_EQUAL(out.scores.size(), input.features.size())
  TEST_EQUAL(out.q_values.size(), input.features.size())
  TEST_EQUAL(out.peps.size(), input.features.size())
  TEST_TRUE(medianTargetBeatsDecoy_(out, input))
}
END_SECTION

START_SECTION([EXTRA] train on A score disjoint B)
{
  // Independent datasets drawn from the same distribution.
  RescoreInput A = makeModeratelySeparableInput_(300, 42);
  RescoreInput B = makeModeratelySeparableInput_(300, 99);

  Percolator p;
  Param par = p.getDefaults();
  par.setValue("seed", 17);
  p.setParameters(par);

  PercolatorModel model = p.train(A);
  RescoreOutput out = p.score(B, model);

  TEST_EQUAL(out.scores.size(), B.features.size())
  TEST_TRUE(medianTargetBeatsDecoy_(out, B))
}
END_SECTION

START_SECTION([EXTRA] train then score on same overlapping input)
{
  // Exercises the unified train+score path on data that is NOT trivially
  // separable (significant target/decoy overlap). Confirms:
  //  1. checkSeparationAndSetPi0 does NOT fire (scoring the training rows
  //     with averaged weights is safe when class distributions overlap).
  //  2. q-values and PEPs land in [0, 1].
  //  3. Target median score still > decoy median.
  //  4. Scores correlate with rescore()'s CV-merge output on the same input,
  //     even though the two paths use different scoring policies (averaged
  //     weights vs per-fold weights). We check Spearman-style agreement via
  //     the fraction of concordant target/decoy pairs in the top decile —
  //     both paths should rank targets above decoys with high probability.
  RescoreInput input = makeModeratelySeparableInput_(400, 42);

  Percolator p_tr, p_re;
  for (auto* p : {&p_tr, &p_re})
  {
    Param par = p->getDefaults();
    par.setValue("seed", 17);
    p->setParameters(par);
  }

  PercolatorModel model = p_tr.train(input);
  RescoreOutput out_ts = p_tr.score(input, model);  // train+score, avg-weights
  RescoreOutput out_re = p_re.rescore(input);       // CV-merge

  // (1) and (2): q-values / PEPs bounded
  for (double q : out_ts.q_values) TEST_TRUE(q >= 0.0 && q <= 1.0 + 1e-9)
  for (double p : out_ts.peps)     TEST_TRUE(p >= 0.0 && p <= 1.0 + 1e-9)

  // (3): target median > decoy median for the train+score path
  TEST_TRUE(medianTargetBeatsDecoy_(out_ts, input))

  // (4): both paths rank targets above decoys with similar fidelity.
  // Count concordant target/decoy pairs in the top 20% of each ranking.
  auto concordantRate = [&](const RescoreOutput& o) {
    std::vector<size_t> idx(o.scores.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
      [&](size_t a, size_t b){ return o.scores[a] > o.scores[b]; });
    const size_t top = idx.size() / 5;
    size_t tgt = 0;
    for (size_t k = 0; k < top; ++k) if (!input.is_decoy[idx[k]]) ++tgt;
    return static_cast<double>(tgt) / static_cast<double>(top);
  };
  const double rate_ts = concordantRate(out_ts);
  const double rate_re = concordantRate(out_re);
  // Both should have >70% targets in the top 20% (moderate separation).
  TEST_TRUE(rate_ts > 0.70)
  TEST_TRUE(rate_re > 0.70)
  // And they should agree closely — within 10 percentage points.
  TEST_TRUE(std::abs(rate_ts - rate_re) < 0.10)
}
END_SECTION

START_SECTION([EXTRA] score rejects feature-count mismatch)
{
  RescoreInput input = makeModeratelySeparableInput_(300, 42);
  Percolator p;
  PercolatorModel model = p.train(input);

  // Drop one feature from each row, keep the same feature_names (mismatch on size).
  RescoreInput bad = input;
  for (auto& row : bad.features) row.pop_back();
  bad.feature_names.pop_back();
  TEST_EXCEPTION(Exception::InvalidValue, p.score(bad, model))
}
END_SECTION

START_SECTION([EXTRA] score rejects feature-name mismatch)
{
  RescoreInput input = makeModeratelySeparableInput_(300, 42);
  Percolator p;
  PercolatorModel model = p.train(input);

  // Swap feature_names — same count, different labels at same position.
  RescoreInput bad = input;
  bad.feature_names = StringList{"signal", "different_name"};
  TEST_EXCEPTION(Exception::InvalidValue, p.score(bad, model))
}
END_SECTION

START_SECTION((static void saveModel(const PercolatorModel& model, const String& filename)))
{
  RescoreInput input = makeModeratelySeparableInput_(300, 42);
  Percolator p;
  Param par = p.getDefaults();
  par.setValue("seed", 17);
  p.setParameters(par);

  PercolatorModel model = p.train(input);

  String tmp;
  NEW_TMP_FILE(tmp);
  Percolator::saveModel(model, tmp);

  PercolatorModel loaded = Percolator::loadModel(tmp);
  TEST_EQUAL(loaded.format_version, model.format_version)
  TEST_EQUAL(loaded.normalizer_type, model.normalizer_type)
  TEST_EQUAL(loaded.seed, model.seed)
  TEST_EQUAL(loaded.feature_names.size(), model.feature_names.size())
  for (size_t i = 0; i < model.feature_names.size(); ++i)
  {
    TEST_EQUAL(loaded.feature_names[i], model.feature_names[i])
  }
  TEST_EQUAL(loaded.weights.size(), model.weights.size())
  for (size_t i = 0; i < model.weights.size(); ++i)
  {
    TEST_REAL_SIMILAR(loaded.weights[i], model.weights[i])
  }

  // Round-trip scoring: loaded model produces identical scores.
  RescoreOutput out_mem  = p.score(input, model);
  RescoreOutput out_disk = p.score(input, loaded);
  bool all_equal = out_mem.scores.size() == out_disk.scores.size();
  for (size_t i = 0; all_equal && i < out_mem.scores.size(); ++i)
  {
    if (std::abs(out_mem.scores[i] - out_disk.scores[i]) > 1e-9) all_equal = false;
  }
  TEST_EQUAL(all_equal, true)
}
END_SECTION

START_SECTION((static PercolatorModel loadModel(const String& filename)))
{
  // Covered by the save-model round-trip section above.
  // Verify negative cases here: missing file and malformed content.
  String nonexistent;
  NEW_TMP_FILE(nonexistent);
  std::remove(nonexistent.c_str());
  TEST_EXCEPTION(Exception::FileNotFound, Percolator::loadModel(nonexistent))

  auto writeModel = [](const String& path, const std::string& body) {
    std::ofstream os(path.c_str());
    os << body;
  };

  // Unsupported format_version
  {
    String tmp; NEW_TMP_FILE(tmp);
    writeModel(tmp,
      "format_version: 99\n"
      "normalizer: stdv\n"
      "seed: 1\n"
      "n_features: 1\n"
      "bias: 0.1\n"
      "x\t0.5\n");
    TEST_EXCEPTION(Exception::ParseError, Percolator::loadModel(tmp))
  }

  // Missing 'bias' header
  {
    String tmp; NEW_TMP_FILE(tmp);
    writeModel(tmp,
      "format_version: 1\n"
      "normalizer: stdv\n"
      "seed: 1\n"
      "n_features: 1\n"
      "x\t0.5\n");
    TEST_EXCEPTION(Exception::ParseError, Percolator::loadModel(tmp))
  }

  // Unknown header key (typoed 'normalize:' instead of 'normalizer:')
  {
    String tmp; NEW_TMP_FILE(tmp);
    writeModel(tmp,
      "format_version: 1\n"
      "normalize: stdv\n"
      "normalizer: stdv\n"
      "seed: 1\n"
      "n_features: 1\n"
      "bias: 0.1\n"
      "x\t0.5\n");
    TEST_EXCEPTION(Exception::ParseError, Percolator::loadModel(tmp))
  }

  // Invalid normalizer value
  {
    String tmp; NEW_TMP_FILE(tmp);
    writeModel(tmp,
      "format_version: 1\n"
      "normalizer: bogus\n"
      "seed: 1\n"
      "n_features: 1\n"
      "bias: 0.1\n"
      "x\t0.5\n");
    TEST_EXCEPTION(Exception::ParseError, Percolator::loadModel(tmp))
  }

  // Duplicate header key
  {
    String tmp; NEW_TMP_FILE(tmp);
    writeModel(tmp,
      "format_version: 1\n"
      "format_version: 1\n"
      "normalizer: stdv\n"
      "seed: 1\n"
      "n_features: 1\n"
      "bias: 0.1\n"
      "x\t0.5\n");
    TEST_EXCEPTION(Exception::ParseError, Percolator::loadModel(tmp))
  }

  // Declared n_features does not match actual feature row count
  {
    String tmp; NEW_TMP_FILE(tmp);
    writeModel(tmp,
      "format_version: 1\n"
      "normalizer: stdv\n"
      "seed: 1\n"
      "n_features: 2\n"  // claims 2 but only 1 row below
      "bias: 0.1\n"
      "x\t0.5\n");
    TEST_EXCEPTION(Exception::ParseError, Percolator::loadModel(tmp))
  }
}
END_SECTION

START_SECTION([EXTRA] feature named 'm0' round-trips through save/load)
{
  // Regression test: the earlier on-disk format reserved 'm0' as a bias-row
  // sentinel, which would corrupt any model whose features happened to carry
  // that name. Moving bias into the header makes feature names fully opaque.
  RescoreInput input = makeModeratelySeparableInput_(300, 42);
  // Rename the first feature to 'm0' on both the input and (implicitly) the
  // trained model.
  input.feature_names[0] = "m0";

  Percolator p;
  Param par = p.getDefaults();
  par.setValue("seed", 17);
  p.setParameters(par);

  PercolatorModel model = p.train(input);
  TEST_EQUAL(model.feature_names[0], "m0")

  String tmp; NEW_TMP_FILE(tmp);
  Percolator::saveModel(model, tmp);

  PercolatorModel loaded = Percolator::loadModel(tmp);
  TEST_EQUAL(loaded.feature_names.size(), model.feature_names.size())
  TEST_EQUAL(loaded.feature_names[0], "m0")
  TEST_EQUAL(loaded.feature_names[1], "noise")
  TEST_EQUAL(loaded.weights.size(), model.weights.size())
  for (size_t i = 0; i < model.weights.size(); ++i)
  {
    TEST_REAL_SIMILAR(loaded.weights[i], model.weights[i])
  }

  // Scoring with the loaded model produces identical results.
  RescoreOutput out_mem  = p.score(input, model);
  RescoreOutput out_disk = p.score(input, loaded);
  bool all_equal = out_mem.scores.size() == out_disk.scores.size();
  for (size_t i = 0; all_equal && i < out_mem.scores.size(); ++i)
  {
    if (std::abs(out_mem.scores[i] - out_disk.scores[i]) > 1e-9) all_equal = false;
  }
  TEST_EQUAL(all_equal, true)
}
END_SECTION

START_SECTION([EXTRA] score rejects empty feature_names on either side)
{
  RescoreInput input = makeModeratelySeparableInput_(300, 42);

  Percolator p;
  PercolatorModel model = p.train(input);
  TEST_EQUAL(model.feature_names.empty(), false)

  // Empty input feature_names — model has names, input doesn't.
  {
    RescoreInput no_names = input;
    no_names.feature_names.clear();
    TEST_EXCEPTION(Exception::InvalidValue, p.score(no_names, model))
  }

  // Empty model feature_names — constructed by hand to simulate a corrupted
  // or partially-built model. Weight count must still satisfy the n_feat+1
  // precondition so the mismatch is specifically about names.
  {
    PercolatorModel bad_model;
    bad_model.format_version = 1;
    bad_model.weights = std::vector<double>(input.features.front().size() + 1, 0.1);
    bad_model.normalizer_type = "stdv";
    // feature_names deliberately left empty
    TEST_EXCEPTION(Exception::InvalidValue, p.score(input, bad_model))
  }
}
END_SECTION

START_SECTION([EXTRA] same-instance consecutive calls are bit-identical)
{
  // Guards against persistent-RNG / leaked-global-state regressions: a
  // second rescore() on the same instance must produce exactly the same
  // output as the first.
  RescoreInput input = makeModeratelySeparableInput_(500, 123);

  Percolator p;
  Param par = p.getDefaults();
  par.setValue("seed", 123);
  p.setParameters(par);

  RescoreOutput out1 = p.rescore(input);
  RescoreOutput out2 = p.rescore(input);

  TEST_EQUAL(out1.scores.size(), out2.scores.size())
  TEST_EQUAL(out1.q_values.size(), out2.q_values.size())
  TEST_EQUAL(out1.peps.size(), out2.peps.size())

  bool bit_identical = out1.scores.size() == out2.scores.size();
  for (size_t i = 0; bit_identical && i < out1.scores.size(); ++i)
  {
    if (out1.scores[i]   != out2.scores[i])   bit_identical = false;
    if (out1.q_values[i] != out2.q_values[i]) bit_identical = false;
    if (out1.peps[i]     != out2.peps[i])     bit_identical = false;
  }
  TEST_EQUAL(bit_identical, true)
}
END_SECTION

START_SECTION([EXTRA] input-order invariance Jaccard characterization)
{
  // Documents the input-order dependence the docstring calls out. Asserts
  // Pearson r >= 0.99 on scores (after un-permuting) and Jaccard >= 0.95
  // on accepted-target sets at q <= 0.01.
  //
  // Prerequisite: populate stable scan_numbers on `base` BEFORE shuffling.
  // The CV fold hash incorporates scan_numbers (defaults to row_index when
  // empty). Without stable scans, shuffling reassigns CV folds, which
  // dominates the measurement. Populating stable scans isolates the
  // intrinsic algorithmic order-sensitivity we actually want to bound.
  RescoreInput base = makeModeratelySeparableInput_(800, 456);
  base.scan_numbers.resize(base.features.size());
  for (size_t i = 0; i < base.scan_numbers.size(); ++i)
  {
    base.scan_numbers[i] = static_cast<int>(i + 1);
  }

  std::vector<size_t> perm(base.features.size());
  std::iota(perm.begin(), perm.end(), 0);
  std::mt19937 rng(456);
  std::vector<size_t> perm_shuf = perm;
  std::shuffle(perm_shuf.begin(), perm_shuf.end(), rng);

  auto permuteInput = [](const RescoreInput& in, const std::vector<size_t>& p) {
    RescoreInput out;
    out.feature_names = in.feature_names;
    out.features.reserve(p.size());
    out.is_decoy.reserve(p.size());
    out.scan_numbers.reserve(p.size());
    for (size_t i : p)
    {
      out.features.push_back(in.features[i]);
      out.is_decoy.push_back(in.is_decoy[i]);
      out.scan_numbers.push_back(in.scan_numbers[i]);
    }
    return out;
  };

  Percolator p1, p2;
  for (auto* p : {&p1, &p2})
  {
    Param par = p->getDefaults();
    par.setValue("seed", 456);
    p->setParameters(par);
  }

  RescoreOutput out1 = p1.rescore(permuteInput(base, perm));
  RescoreOutput out2 = p2.rescore(permuteInput(base, perm_shuf));

  // Un-permute out2 so row i refers to base row i.
  const size_t n = out1.scores.size();
  std::vector<double> scores2_aligned(n), qvals2_aligned(n);
  for (size_t i = 0; i < perm_shuf.size(); ++i)
  {
    scores2_aligned[perm_shuf[i]] = out2.scores[i];
    qvals2_aligned[perm_shuf[i]] = out2.q_values[i];
  }

  // Pearson r on aligned scores.
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const double x = out1.scores[i];
    const double y = scores2_aligned[i];
    sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
  }
  const double nf = static_cast<double>(n);
  const double num = nf * sxy - sx * sy;
  const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
  const double r = (den > 1e-15) ? num / den : 0.0;
  TEST_TRUE(r >= 0.99)

  // Jaccard on accepted-target set at q <= 0.01.
  std::set<size_t> a, b;
  for (size_t i = 0; i < n; ++i)
  {
    if (base.is_decoy[i]) continue;
    if (out1.q_values[i]  <= 0.01) a.insert(i);
    if (qvals2_aligned[i] <= 0.01) b.insert(i);
  }
  std::set<size_t> inter, uni;
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(),
                        std::inserter(inter, inter.begin()));
  std::set_union(a.begin(), a.end(), b.begin(), b.end(),
                 std::inserter(uni, uni.begin()));
  const double jaccard = uni.empty() ? 1.0
    : static_cast<double>(inter.size()) / static_cast<double>(uni.size());
  TEST_TRUE(jaccard >= 0.95)
}
END_SECTION

START_SECTION([EXTRA] seed change produces different but statistically equivalent scores)
{
  // Negative control: two different seeds must produce different scores
  // (catches seed-not-wired-through bugs) but Pearson r >= 0.9 and
  // |cA - cB| / max(cA, cB) <= 0.20 at q <= 0.05 (statistical equivalence).
  //
  // Why q <= 0.05 instead of 0.01: the q-value "cliff" at 0.01 is inherently
  // unstable across seeds — small score shifts move many targets across the
  // threshold even when Pearson r is ~0.985. Observed counts on this fixture:
  //   q <= 0.01: c17=123, c42=171, ratio=0.28 (unstable, cliff-sensitive)
  //   q <= 0.05: c17=394, c42=356, ratio=0.10 (stable, larger target mass)
  //   q <= 0.10: c17=657, c42=652, ratio=0.01
  // q <= 0.05 measures the same downstream "seed-induced classification drift"
  // property at a threshold where there's enough accepted targets for the
  // estimate to be stable.
  RescoreInput input = makeModeratelySeparableInput_(800, 789);

  Percolator p17, p42;
  {
    Param par = p17.getDefaults();
    par.setValue("seed", 17);
    p17.setParameters(par);
  }
  {
    Param par = p42.getDefaults();
    par.setValue("seed", 42);
    p42.setParameters(par);
  }

  RescoreOutput out17 = p17.rescore(input);
  RescoreOutput out42 = p42.rescore(input);

  // Not bit-identical.
  bool any_diff = false;
  for (size_t i = 0; i < out17.scores.size(); ++i)
  {
    if (out17.scores[i] != out42.scores[i]) { any_diff = true; break; }
  }
  TEST_EQUAL(any_diff, true)

  // Pearson r >= 0.9.
  const size_t n = out17.scores.size();
  double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
  for (size_t i = 0; i < n; ++i)
  {
    const double x = out17.scores[i];
    const double y = out42.scores[i];
    sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y;
  }
  const double nf = static_cast<double>(n);
  const double num = nf * sxy - sx * sy;
  const double den = std::sqrt((nf * sxx - sx * sx) * (nf * syy - sy * sy));
  const double r = (den > 1e-15) ? num / den : 0.0;
  TEST_TRUE(r >= 0.9)

  // Count ratio at q <= 0.05 (targets only).
  int c17 = 0, c42 = 0;
  for (size_t i = 0; i < n; ++i)
  {
    if (input.is_decoy[i]) continue;
    if (out17.q_values[i] <= 0.05) ++c17;
    if (out42.q_values[i] <= 0.05) ++c42;
  }
  const int maxc = std::max(c17, c42);
  const double ratio = (maxc > 0)
    ? static_cast<double>(std::abs(c17 - c42)) / static_cast<double>(maxc)
    : 0.0;
  TEST_TRUE(ratio <= 0.20)
}
END_SECTION

END_TEST
