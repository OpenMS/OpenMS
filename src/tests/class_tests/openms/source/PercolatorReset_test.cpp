// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/ID/PercolatorReset.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{
  // deterministic pseudo-random numbers (the same on every platform)
  struct Rng
  {
    uint64_t s;
    double uniform() // [0, 1)
    {
      s = s * 6364136223846793005ULL + 1442695040888963407ULL;
      return double(s >> 11) * (1.0 / 9007199254740992.0);
    }
    double normal() { return (uniform() + uniform() + uniform() + uniform() - 2.0) * 1.7320508; }
  };

  std::string randomSequence(Rng& rng, size_t len)
  {
    static const std::string aa = "ACDEFGHILMNPQSTVWY";
    std::string s;
    for (size_t i = 0; i < len; ++i) s += aa[static_cast<size_t>(rng.uniform() * aa.size())];
    return s;
  }

  // Simulated search: correct target PSMs (features shifted), incorrect targets and decoys (same null distribution).
  // Decoys are the reversed bodies of target stems with K or R as C-terminal residue (protein-reversal-like), each
  // spectrum has one hit; @p correct_truth receives per row whether it is a correct target.
  PercolatorReset::Input simulate(uint64_t seed, size_t correct, size_t incorrect, size_t decoys, std::vector<bool>& correct_truth,
                                  double shift = 4.5)
  {
    Rng rng{seed};
    PercolatorReset::Input in;
    in.feature_names = {"score", "signal", "noise", "enzInt"};
    std::vector<std::string> bodies;
    for (size_t i = 0; i < correct + incorrect; ++i) bodies.push_back(randomSequence(rng, 7 + i % 9));
    size_t spectrum = 0;
    auto add = [&](const std::string& stem, bool decoy, bool good)
    {
      const double signal = rng.normal() + (good ? shift : 0.0);
      const double score = 0.6 * signal + 0.8 * rng.normal();
      in.features.push_back({score, signal, rng.normal(), decoy ? 0.0 : 1.0});
      in.is_decoy.push_back(decoy);
      in.search_scores.push_back(score);
      in.spectrum.push_back(spectrum++);
      in.stems.push_back(stem);
      correct_truth.push_back(good);
    };
    for (size_t i = 0; i < correct; ++i) add(bodies[i] + "K", false, true);
    for (size_t i = correct; i < correct + incorrect; ++i) add(bodies[i] + "K", false, false);
    for (size_t i = 0; i < decoys; ++i)
    {
      std::string body = bodies[i % bodies.size()];
      std::reverse(body.begin(), body.end());
      add(body + (i % 2 == 0 ? "K" : "R"), true, false);

    }
    return in;
  }
}

START_TEST(PercolatorReset, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION(static std::vector<double> seqStepQValues(const std::vector<bool>& is_decoy, double factor, double plus_one))
{
  // T T D T D: fdp = 1/1, 1/2, 2/2, 2/3, 3/3 -> q = 0.5, 0.5, 2/3, 2/3, 1
  const std::vector<double> q = PercolatorReset::seqStepQValues({false, false, true, false, true}, 1.0, 1.0);
  TEST_EQUAL(q.size(), 5)
  TEST_REAL_SIMILAR(q[0], 0.5)
  TEST_REAL_SIMILAR(q[1], 0.5)
  TEST_REAL_SIMILAR(q[2], 2.0 / 3.0)
  TEST_REAL_SIMILAR(q[3], 2.0 / 3.0)
  TEST_REAL_SIMILAR(q[4], 1.0)
  // factor 2 (single-decoy RESET), capped at 1; without the +1
  const std::vector<double> q2 = PercolatorReset::seqStepQValues({false, false, false, false, true}, 2.0, 0.0);
  TEST_REAL_SIMILAR(q2[0], 0.0)
  TEST_REAL_SIMILAR(q2[3], 0.0)
  TEST_REAL_SIMILAR(q2[4], 0.5)
}
END_SECTION

START_SECTION(static bool isTrainingDecoy(std::uint64_t salt, const std::string& stem))
{
  Rng rng{7};
  size_t training = 0, same = 0;
  const size_t n = 20000;
  for (size_t i = 0; i < n; ++i)
  {
    const std::string stem = randomSequence(rng, 10);
    const bool a = PercolatorReset::isTrainingDecoy(42, stem);
    TEST_EQUAL(PercolatorReset::isTrainingDecoy(42, stem), a) // a function of (salt, stem)
    training += a ? 1 : 0;
    same += a == PercolatorReset::isTrainingDecoy(43, stem) ? 1 : 0;
  }
  // Bernoulli(1/2): 20000 coins, 5 standard deviations = 354
  TEST_EQUAL(std::abs(static_cast<double>(training) - n / 2.0) < 354.0, true)
  // another salt gives independent coins
  TEST_EQUAL(std::abs(static_cast<double>(same) - n / 2.0) < 354.0, true)
}
END_SECTION

START_SECTION(static Output run(const Input& input, const Options& options))
{
  std::vector<bool> truth;
  const PercolatorReset::Input in = simulate(11, 1500, 1500, 3000, truth);
  PercolatorReset::Options opt;
  opt.salt = 12345;
  opt.threads = 2;
  opt.iterations = 5;          // default 20
  opt.representative = "model"; // experimental option, exercised by the runs below (incl. the masking test)
  const PercolatorReset::Output out = PercolatorReset::run(in, opt);
  TEST_EQUAL(out.scores.size(), in.features.size())
  TEST_EQUAL(out.stem_q_values.size(), in.features.size())
  TEST_EQUAL(out.target_stems, 3000)
  TEST_EQUAL(out.decoy_stems, 3000)
  TEST_EQUAL(out.training_decoys + out.estimating_decoys, 3000)
  // enzInt is zeroed: not a training feature
  TEST_EQUAL(std::find(out.used_features.begin(), out.used_features.end(), "enzInt") == out.used_features.end(), true)
  TEST_EQUAL(out.rounds.size(), 5)
  // power and (loosely) the false discovery proportion at 1 %
  size_t disc = 0, false_disc = 0;
  for (size_t i = 0; i < in.features.size(); ++i)
  {
    if (!in.is_decoy[i] && out.stem_q_values[i] <= 0.01)
    {
      ++disc;
      if (!truth[i]) ++false_disc;
    }
  }
  TEST_EQUAL(disc, out.discoveries_1pct)
  TEST_EQUAL(disc > 600, true)
  TEST_EQUAL(false_disc <= disc / 20, true)
  // training decoys are not ranked
  for (size_t i = 0; i < in.features.size(); ++i)
  {
    if (out.roles[i] == PercolatorReset::Role::TRAINING_DECOY)
    {
      TEST_REAL_SIMILAR(out.stem_q_values[i], 1.0)
      TEST_REAL_SIMILAR(out.psm_q_values[i], 1.0)
    }
    TEST_EQUAL(out.roles[i] == PercolatorReset::Role::TARGET, !in.is_decoy[i])
    TEST_EQUAL(out.psm_peps[i] >= 0.0 && out.psm_peps[i] <= 1.0, true)
  }
  // q-values are monotone in the final score among ranked stems
  std::vector<std::pair<double, double>> sq;
  for (const PercolatorReset::Stem& s : out.stems)
  {
    if (s.considered && s.role != PercolatorReset::Role::TRAINING_DECOY) sq.emplace_back(s.score, s.q_value);
  }
  std::sort(sq.begin(), sq.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
  bool monotone = true;
  for (size_t k = 1; k < sq.size(); ++k) monotone = monotone && sq[k].second >= sq[k - 1].second;
  TEST_EQUAL(monotone, true)

  // reproducible with the same salt, independent of the number of threads
  PercolatorReset::Options opt1 = opt;
  opt1.threads = 1;
  const PercolatorReset::Output again = PercolatorReset::run(in, opt1);
  TEST_EQUAL(again.scores == out.scores, true)
  TEST_EQUAL(again.stem_q_values == out.stem_q_values, true)

  // masking: the training sees the estimating decoys only as pseudo-targets, so relabelling them as targets leaves
  // every score unchanged (only the final estimate changes)
  PercolatorReset::Input relabelled = in;
  for (size_t i = 0; i < in.features.size(); ++i)
  {
    if (out.roles[i] == PercolatorReset::Role::ESTIMATING_DECOY) relabelled.is_decoy[i] = false;
  }
  const PercolatorReset::Output masked = PercolatorReset::run(relabelled, opt);
  TEST_EQUAL(masked.scores == out.scores, true)
  TEST_EQUAL(masked.estimating_decoys, 0)
  TEST_EQUAL(masked.training_decoys, out.training_decoys)

  // another salt: other coins, similar power
  PercolatorReset::Options opt2 = opt;
  opt2.salt = 99;
  const PercolatorReset::Output other = PercolatorReset::run(in, opt2);
  TEST_EQUAL(other.training_decoys != out.training_decoys || other.scores != out.scores, true)
  TEST_EQUAL(std::abs(static_cast<double>(other.discoveries_1pct) - static_cast<double>(out.discoveries_1pct)) < 200.0, true)

  // audit: a symmetric feature set is not flagged; a feature that is 1 for targets and 0 for decoys is
  opt.audit = true;
  const PercolatorReset::Output audited = PercolatorReset::run(in, opt);
  TEST_EQUAL(audited.audit.size(), 4)
  std::map<std::string, PercolatorReset::AuditFeature> af;
  for (const auto& a : audited.audit) af[a.name] = a;
  TEST_EQUAL(af["noise"].flagged, false)
  TEST_EQUAL(af["enzInt"].flagged, true) // the leaky feature (reported, zeroed by default)
  TEST_EQUAL(audited.audit_multivariate_auc < 0.6, true) // enzInt is not a training feature (the low-score null
                                                        // proxy holds some correct targets: not exactly 0.5)
  TEST_EQUAL(audited.scores == out.scores, true)          // the audit only reports

  // the MLP learner runs and finds discoveries
  opt.learner = "mlp";
  opt.mlp_members = 2;
  opt.mlp_epochs = 5;
  opt.iterations = 2;
  const PercolatorReset::Output mlp = PercolatorReset::run(in, opt);
  TEST_EQUAL(mlp.discoveries_1pct > 600, true)

  // FDR control (weak signal, many nulls): mean false discovery proportion over 20 simulated data sets and salts at
  // q <= 0.05 stays near or below 0.05 (expected FDP of RESET <= 0.05; 0.075 leaves room for the sampling error)
  {
    double fdp_sum = 0.0;
    for (uint64_t rep = 0; rep < 20; ++rep)
    {
      std::vector<bool> t2;
      const PercolatorReset::Input sim = simulate(1000 + rep, 400, 2600, 3000, t2, 2.5);
      PercolatorReset::Options o2;
      o2.salt = 777 + rep;
      o2.iterations = 3;
      const PercolatorReset::Output r2 = PercolatorReset::run(sim, o2);
      size_t d2 = 0, f2 = 0;
      for (size_t i = 0; i < sim.features.size(); ++i)
      {
        if (!sim.is_decoy[i] && r2.stem_q_values[i] <= 0.05) { ++d2; if (!t2[i]) ++f2; }
      }
      fdp_sum += d2 > 0 ? static_cast<double>(f2) / static_cast<double>(d2) : 0.0;
    }
    STATUS("mean FDP at 5 %: " << fdp_sum / 20.0)
    TEST_EQUAL(fdp_sum / 20.0 < 0.075, true)
  }

  opt.learner = "svm";
  opt.iterations = 5;
  // review r1: a target and a decoy hit of the same sequence form ONE stem (label-blind stem construction); its label
  // is the label of its best PSM
  {
    PercolatorReset::Input col = in;
    col.stems[0] = "PEPTIDEK";                  // a target row
    col.stems[in.features.size() - 1] = "PEPTIDEK"; // a decoy row
    const PercolatorReset::Output r = PercolatorReset::run(col, opt);
    Size n = 0;
    for (const PercolatorReset::Stem& st : r.stems) n += st.stem == "PEPTIDEK" ? 1 : 0;
    TEST_EQUAL(n, 1)
  }

  // review r2: a target row never inherits the q-value of a decoy stem of its sequence: a strong decoy PSM and a weak
  // target PSM of one sequence -> the stem is a decoy; the target row gets peptide q-value 1 and role TARGET
  {
    PercolatorReset::Input col = in;
    const size_t t = 0, d = in.features.size() - 1;
    col.stems[t] = "PEPTIDEK";
    col.stems[d] = "PEPTIDEK";
    col.search_scores[t] = -50.0;
    col.features[t][0] = -50.0;
    col.search_scores[d] = 50.0;
    col.features[d] = {50.0, 50.0, 0.0, 0.0};
    const PercolatorReset::Output r = PercolatorReset::run(col, opt);
    TEST_EQUAL(r.roles[t] == PercolatorReset::Role::TARGET, true)
    TEST_REAL_SIMILAR(r.stem_q_values[t], 1.0)
    // the exported target discoveries are exactly the target stems with q <= 0.01
    std::set<std::string> exported, accepted;
    for (size_t i = 0; i < col.features.size(); ++i)
    {
      if (!col.is_decoy[i] && r.stem_q_values[i] <= 0.01) exported.insert(col.stems[i]);
    }
    for (const PercolatorReset::Stem& st : r.stems)
    {
      if (!st.decoy && st.q_value <= 0.01) accepted.insert(st.stem);
    }
    TEST_EQUAL(exported == accepted, true)
  }

  // review r1: without a salt, run() draws one and records it
  {
    PercolatorReset::Options no_salt = opt;
    no_salt.salt.reset();
    const PercolatorReset::Output a = PercolatorReset::run(in, no_salt);
    const PercolatorReset::Output b = PercolatorReset::run(in, no_salt);
    TEST_EQUAL(a.salt != b.salt, true)
    PercolatorReset::Options replay = opt;
    replay.salt = a.salt;
    TEST_EQUAL(PercolatorReset::run(in, replay).scores == a.scores, true)
  }

  // review r1: invalid options and features whose standardization overflows are rejected
  {
    PercolatorReset::Options o = opt;
    o.train_fdr = -std::numeric_limits<double>::infinity();
    TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(in, o))
    o = opt; o.train_fdr = 1.0;
    TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(in, o))
    o = opt; o.c_values = {1.0, 0.0};
    TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(in, o))
    o = opt; o.c_values = {1.0, std::numeric_limits<double>::infinity()};
    TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(in, o))
    o = opt; o.learner = "mlp"; o.mlp_members = 0;
    TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(in, o))
    o = opt; o.iterations = 0;
    TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(in, o))
    PercolatorReset::Input huge = in;
    for (size_t i = 0; i < huge.features.size(); ++i) huge.features[i][2] = (i % 2 == 0 ? 1.7e308 : -1.7e308);
    TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(huge, opt))
  }

  // the paper's representative (best search score) runs as well
  {
    PercolatorReset::Options o = opt;
    o.representative = "search";
    const PercolatorReset::Output rs = PercolatorReset::run(in, o);
    TEST_EQUAL(rs.discoveries_1pct > 600, true)
    // masking also with the default representative
    PercolatorReset::Input rl = in;
    for (size_t i = 0; i < in.features.size(); ++i)
    {
      if (rs.roles[i] == PercolatorReset::Role::ESTIMATING_DECOY) rl.is_decoy[i] = false;
    }
    TEST_EQUAL(PercolatorReset::run(rl, o).scores == rs.scores, true)
    o.representative = "best";
    TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(in, o))
  }

  // malformed input
  PercolatorReset::Input bad = in;
  bad.features[3][1] = std::numeric_limits<double>::quiet_NaN();
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(bad, opt))
  bad = in;
  bad.stems.pop_back();
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(bad, opt))
  bad = in;
  bad.is_decoy.assign(bad.is_decoy.size(), false);
  TEST_EXCEPTION(Exception::InvalidValue, PercolatorReset::run(bad, opt))
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
