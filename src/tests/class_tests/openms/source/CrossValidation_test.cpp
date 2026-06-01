// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ML/CROSSVALIDATION/CrossValidation.h>

#include <vector>
#include <cmath>
#include <limits>

using namespace OpenMS;
using std::vector;

START_TEST(CrossValidation, "$Id$")

// --- makeKFolds --------------------------------------------------------------
START_SECTION(makeKFolds basic and edge cases)
{
  // Basic: n=5, K=2  -> round-robin [0,2,4], [1,3]
  {
    const auto folds = CrossValidation::makeKFolds(5, 2);
    TEST_EQUAL(folds.size(), 2);
    TEST_EQUAL(folds[0].size(), 3);
    TEST_EQUAL(folds[1].size(), 2);
    TEST_EQUAL(folds[0][0], 0);
    TEST_EQUAL(folds[0][1], 2);
    TEST_EQUAL(folds[0][2], 4);
    TEST_EQUAL(folds[1][0], 1);
    TEST_EQUAL(folds[1][1], 3);
  }

  // K > n gets clamped to n (LOO-like singletons)
  {
    const auto folds = CrossValidation::makeKFolds(5, 7);
    TEST_EQUAL(folds.size(), 5);
    for (Size i = 0; i < folds.size(); ++i)
    {
      TEST_EQUAL(folds[i].size(), 1);
      TEST_EQUAL(folds[i][0], i);
    }
  }

  // Exactly one sample, K==1 → single fold holding {0}
  {
    const auto folds = CrossValidation::makeKFolds(1, 1);
    TEST_EQUAL(folds.size(), 1);
    TEST_EQUAL(folds[0].size(), 1);
    TEST_EQUAL(folds[0][0], 0);
  }

  // K > n gets clamped to n (still a single fold with {0})
  {
    const auto folds = CrossValidation::makeKFolds(1, 5);
    TEST_EQUAL(folds.size(), 1);
    TEST_EQUAL(folds[0].size(), 1);
    TEST_EQUAL(folds[0][0], 0);
  }

  // Invalid inputs
  TEST_EXCEPTION(Exception::InvalidValue, CrossValidation::makeKFolds(0, 1)); // n==0
  TEST_EXCEPTION(Exception::InvalidValue, CrossValidation::makeKFolds(5, 0)); // K==0
}
END_SECTION

// --- gridSearch1D: best candidate selection ---------------------------------
START_SECTION(gridSearch1D selects the true best candidate by score)
{
  // Candidates around the optimum 0.5
  const vector<double> cands{0.2, 0.5, 0.8};

  // Any deterministic folds; content does not change the score logic here
  const auto folds = CrossValidation::makeKFolds(6, 3); // [[0,3],[1,4],[2,5]]

  // Train/eval: append |cand - 0.5| per held-out sample (constant error per sample)
  auto train_eval = [](double cand,
                       const vector<vector<Size>>& flds,
                       vector<double>& abs_errs)
  {
    const double e = std::fabs(cand - 0.5);
    for (const auto& fold : flds)
    {
      for (Size idx : fold) { (void)idx; abs_errs.push_back(e); }
    }
  };

  // Score: mean of abs_errs
  auto score = [](const vector<double>& errs) -> double
  {
    double s = 0.0;
    for (double v : errs) s += v;
    return errs.empty() ? std::numeric_limits<double>::infinity() : s / errs.size();
  };

  const auto result = CrossValidation::gridSearch1D(
    cands.begin(), cands.end(),
    folds, train_eval, score);

  TEST_REAL_SIMILAR(result.first, 0.5); // best candidate
  TEST_REAL_SIMILAR(result.second, 0.0); // zero error at optimum
}
END_SECTION

// --- gridSearch1D: tie-breaking policies ------------------------------------
START_SECTION(gridSearch1D tie-breaking (PreferLarger / PreferSmaller / PreferAny))
{
  const vector<double> cands{0.4, 0.6}; // symmetric around 0.5 -> equal scores
  const auto folds = CrossValidation::makeKFolds(4, 2);

  auto train_eval = [](double cand,
                       const vector<vector<Size>>& flds,
                       vector<double>& abs_errs)
  {
    const double e = std::fabs(cand - 0.5); // 0.1 for both
    for (const auto& fold : flds)
      for (Size idx : fold) { (void)idx; abs_errs.push_back(e); }
  };
  auto score = [](const vector<double>& errs)
  {
    double s = 0.0;
    for (double v : errs) s += v;
    return s / errs.size();
  };

  // PreferLarger: picks 0.6 on tie
  {
    const auto [cand, sc] = CrossValidation::gridSearch1D(
      cands.begin(), cands.end(), folds, train_eval, score,
      1e-12, CrossValidation::CandidateTieBreak::PreferLarger);
    TEST_REAL_SIMILAR(cand, 0.6);
    TEST_REAL_SIMILAR(sc, 0.1);
  }

  // PreferSmaller: picks 0.4 on tie
  {
    const auto [cand, sc] = CrossValidation::gridSearch1D(
      cands.begin(), cands.end(), folds, train_eval, score,
      1e-12, CrossValidation::CandidateTieBreak::PreferSmaller);
    TEST_REAL_SIMILAR(cand, 0.4);
    TEST_REAL_SIMILAR(sc, 0.1);
  }

  // PreferAny: keeps first encountered (0.4)
  {
    const auto [cand, sc] = CrossValidation::gridSearch1D(
      cands.begin(), cands.end(), folds, train_eval, score,
      1e-12, CrossValidation::CandidateTieBreak::PreferAny);
    TEST_REAL_SIMILAR(cand, 0.4);
    TEST_REAL_SIMILAR(sc, 0.1);
  }
}
END_SECTION

// --- gridSearch1D: tie tolerance behavior -----------------------------------
START_SECTION(gridSearch1D respects tie tolerance)
{
  // First candidate is within tie tolerance to perfect optimum; second is exact optimum.
  const double near = 0.5 - 5e-13; // |error| = 5e-13
  const vector<double> cands{near, 0.5};
  const auto folds = CrossValidation::makeKFolds(3, 3);

  auto train_eval = [](double cand,
                       const vector<vector<Size>>& flds,
                       vector<double>& abs_errs)
  {
    const double e = std::fabs(cand - 0.5);
    for (const auto& fold : flds)
      for (Size idx : fold) { (void)idx; abs_errs.push_back(e); }
  };
  auto score = [](const vector<double>& errs)
  {
    double s = 0.0;
    for (double v : errs) s += v;
    return s / errs.size();
  };

  // With default tie_tol=1e-12, (5e-13 vs 0) are considered a tie → PreferLarger picks 0.5
  {
    const auto [cand, sc] = CrossValidation::gridSearch1D(
      cands.begin(), cands.end(), folds, train_eval, score);
    TEST_REAL_SIMILAR(cand, 0.5);
    TEST_REAL_SIMILAR(sc, 0.0);
  }

  // If we tighten tie_tol below 5e-13, the smaller score (0.0) must win regardless of tie-break policy
  {
    const auto [cand, sc] = CrossValidation::gridSearch1D(
      cands.begin(), cands.end(), folds, train_eval, score,
      1e-13, CrossValidation::CandidateTieBreak::PreferSmaller);
    TEST_REAL_SIMILAR(cand, 0.5);
    TEST_REAL_SIMILAR(sc, 0.0);
  }
}
END_SECTION

// --- gridSearch1D: works with single-sample LOO (n==1) -----------------------
START_SECTION(gridSearch1D works with single-sample LOO (n==1))
{
  const std::vector<double> cands{0.2, 0.5};
  const auto folds = CrossValidation::makeKFolds(1, 1); // [[0]]

  // Append one error per held-out idx (here exactly one)
  auto train_eval = [](double cand,
                       const std::vector<std::vector<Size>>& flds,
                       std::vector<double>& abs_errs)
  {
    const double e = std::fabs(cand - 0.3); // optimum at cand=0.3
    for (const auto& fold : flds)
      for (Size idx : fold) { (void)idx; abs_errs.push_back(e); }
  };

  // Mean absolute error
  auto score = [](const std::vector<double>& errs) -> double
  {
    double s = 0.0;
    for (double v : errs) s += v;
    return errs.empty() ? std::numeric_limits<double>::infinity() : s / errs.size();
  };

  const auto [best_cand, best_score] =
    CrossValidation::gridSearch1D(cands.begin(), cands.end(),
                                  folds, train_eval, score);

  TEST_REAL_SIMILAR(best_cand, 0.2);                 // closer to 0.3 than 0.5
  TEST_REAL_SIMILAR(best_score, std::fabs(0.2 - 0.3));
}
END_SECTION

// --- gridSearch1D: empty candidate range throws -----------------------------
START_SECTION(gridSearch1D throws on empty candidate range)
{
  const vector<double> empty;
  const auto folds = CrossValidation::makeKFolds(3, 3);

  auto train_eval = [](double /*cand*/,
                       const vector<vector<Size>>& /*flds*/,
                       vector<double>& /*abs_errs*/) {};
  auto score = [](const vector<double>& errs)
  {
    return errs.empty() ? std::numeric_limits<double>::infinity() : errs.front();
  };

  TEST_EXCEPTION(Exception::InvalidRange,
                 CrossValidation::gridSearch1D(empty.begin(), empty.end(),
                                               folds, train_eval, score));
}
END_SECTION

END_TEST
