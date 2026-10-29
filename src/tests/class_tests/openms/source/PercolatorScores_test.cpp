// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

// Vendored Percolator (OpenMS::Internal::Percolator): this test links the
// static OpenMS_Percolator library directly, its symbols are hidden in libOpenMS.
#include "DataSet.h"
#include "FeatureNames.h"
#include "PSMDescription.h"
#include "ScoreHolder.h"
#include "Scores.h"
#include "ssl.h"

#include <algorithm>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;
using OpenMS::Internal::Percolator::AlgIn;
using OpenMS::Internal::Percolator::DataSet;
using OpenMS::Internal::Percolator::LabelType;
using OpenMS::Internal::Percolator::OrderScanHash;
using OpenMS::Internal::Percolator::PSMDescription;
using OpenMS::Internal::Percolator::ScoreHolder;
using OpenMS::Internal::Percolator::Scores;

namespace
{
  // The PSMs of a test, scored by their single feature. PSMDescription does not own its features.
  struct PSMStore
  {
    vector<unique_ptr<PSMDescription>> psms;
    vector<unique_ptr<double[]>> features;

    void add(Scores& scores, unsigned file, unsigned scan, const string& id, double score, LabelType label)
    {
      features.emplace_back(new double[1]{score});
      auto psm = make_unique<PSMDescription>();
      psm->setId(id);
      psm->specFileNr = file;
      psm->scan = scan;
      psm->expMass = 1000.0 + scan;
      psm->features = features.back().get();
      scores.addScoreHolder(ScoreHolder(0.0, label, psm.get()));
      psms.push_back(std::move(psm));
    }
  };

  // Scores and q-values as CrossValidation::doStep computes them before it selects the positives.
  void scoreAndQvals(Scores& scores, double fdr)
  {
    scores.recalculateSizes();
    vector<double> w = {1.0, 0.0}; // weight of the feature, bias
    scores.calcScoresAndQvals(w, fdr, true);
  }

  // The ids of the positive training examples, in their order.
  vector<string> positiveIds(Scores& scores, double fdr, bool train_best_positive)
  {
    AlgIn data(scores.size(), 2);
    scores.generateNegativeTrainingSet(data, 1.0);
    scores.generatePositiveTrainingSet(data, fdr, 1.0, train_best_positive);
    map<const double*, string> id_of;
    for (const ScoreHolder& sh : scores) id_of[sh.pPSM->features] = sh.pPSM->getId();
    vector<string> ids;
    for (int i = data.negatives; i < data.m; ++i) ids.push_back(id_of[data.vals[i]]);
    return ids;
  }

  vector<const PSMDescription*> psmOrder(const Scores& scores)
  {
    vector<const PSMDescription*> order;
    for (const ScoreHolder& sh : scores) order.push_back(sh.pPSM);
    return order;
  }
}

START_TEST(PercolatorScores, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DataSet::resetFeatureNames();
DataSet::getFeatureNames().insertFeature("score");
DataSet::getFeatureNames().initFeatures();

const double fdr = 0.01;

START_SECTION([EXTRA] generatePositiveTrainingSet with trainBestPositive: the best-scoring target per spectrum)
{
  // Two files with the same scan numbers. Each spectrum has three targets above all decoys (q = 0),
  // one target below them (q > fdr) and two decoys, added with the best target last.
  // Issue #10720: the PSM kept per spectrum was arbitrary, and a spectrum whose kept PSM failed the
  // FDR was missing from the positives.
  Scores scores(false);
  PSMStore store;
  vector<string> best_ids; // by decreasing score
  for (unsigned file = 0; file < 2; ++file)
  {
    for (unsigned scan = 1; scan <= 40; ++scan)
    {
      const double k = 0.1 * (file * 40 + scan);
      const string spectrum = "f" + to_string(file) + "_s" + to_string(scan);
      store.add(scores, file, scan, spectrum + "_low", -100.0 - k, LabelType::TARGET);
      store.add(scores, file, scan, spectrum + "_d1", -10.0 - k, LabelType::DECOY);
      store.add(scores, file, scan, spectrum + "_third", 20.0 + k, LabelType::TARGET);
      store.add(scores, file, scan, spectrum + "_second", 50.0 + k, LabelType::TARGET);
      store.add(scores, file, scan, spectrum + "_d2", -20.0 - k, LabelType::DECOY);
      store.add(scores, file, scan, spectrum + "_best", 100.0 + k, LabelType::TARGET);
      best_ids.insert(best_ids.begin(), spectrum + "_best");
    }
  }
  scoreAndQvals(scores, fdr);

  // input sorted by score, as in CrossValidation::doStep
  const vector<const PSMDescription*> before = psmOrder(scores);
  TEST_TRUE(positiveIds(scores, fdr, true) == best_ids)
  // scores_ is not changed (std::unique overwrote PSMs in it)
  TEST_TRUE(psmOrder(scores) == before)

  // input in spectrum hash order, as CrossValidation::mergeCpCnPairs gets it (nested-xval-bins > 1)
  std::sort(scores.begin(), scores.end(), OrderScanHash());
  const vector<const PSMDescription*> hashed = psmOrder(scores);
  TEST_TRUE(positiveIds(scores, fdr, true) == best_ids)
  TEST_TRUE(psmOrder(scores) == hashed)

  // without trainBestPositive: every target that passes the FDR, in the order of scores_
  vector<string> passing;
  for (const ScoreHolder& sh : scores)
  {
    if (sh.isTarget() && sh.q <= fdr) passing.push_back(sh.pPSM->getId());
  }
  TEST_EQUAL(passing.size(), 240u)
  TEST_TRUE(positiveIds(scores, fdr, false) == passing)
}
END_SECTION

START_SECTION([EXTRA] generatePositiveTrainingSet with trainBestPositive: equal scores)
{
  // Two targets of a spectrum with the same, best score: the one with the smaller PSM id, in
  // either order of the input, so that every standard library selects the same PSM.
  for (bool a_first : {true, false})
  {
    Scores scores(false);
    PSMStore store;
    for (unsigned scan = 1; scan <= 20; ++scan)
    {
      const string spectrum = "s" + to_string(scan);
      if (a_first) store.add(scores, 0, scan, spectrum + "_a", 10.0 + scan, LabelType::TARGET);
      store.add(scores, 0, scan, spectrum + "_b", 10.0 + scan, LabelType::TARGET);
      if (!a_first) store.add(scores, 0, scan, spectrum + "_a", 10.0 + scan, LabelType::TARGET);
      store.add(scores, 0, scan, spectrum + "_c", 5.0 + scan, LabelType::TARGET);
      store.add(scores, 0, scan, spectrum + "_d", -5.0 - scan, LabelType::DECOY);
    }
    scoreAndQvals(scores, fdr);
    vector<string> expected;
    for (unsigned scan = 20; scan >= 1; --scan) expected.push_back("s" + to_string(scan) + "_a");
    TEST_TRUE(positiveIds(scores, fdr, true) == expected)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
