// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

// Vendored Percolator (OpenMS::Internal::Percolator): this test links the
// static OpenMS_Percolator library directly, its symbols are hidden in libOpenMS.
#include "DataSet.h"
#include "FeatureNames.h"
#include "Scores.h"
#include "SetHandler.h"
#include "ssl.h"

#include <cstddef>
#include <map>
#include <set>
#include <vector>

using namespace OpenMS::Internal::Percolator;

namespace
{
  // n_spectra spectra with hits_per_spectrum target hits each (one feature: the score) and as many decoys
  // (one per spectrum, scoring below every target); the target hit h of spectrum s scores
  // 1000 + 100 * s - h * 3.7 (h = 0 the best); the hits are registered from the worst to the best.
  struct Input
  {
    SetHandler handler{0};
    std::map<const double*, std::pair<int, int>> hit_of;  // feature row -> (spectrum, hit), targets
  };

  void fill(Input& in, int n_spectra, int hits_per_spectrum)
  {
    FeatureMemoryPool& pool = in.handler.getFeaturePool();
    pool.createPool(1);
    DataSet::resetFeatureNames();
    FeatureNames& names = DataSet::getFeatureNames();
    names.insertFeature("score");
    names.initFeatures();
    auto* targets = new DataSet();
    auto* decoys = new DataSet();
    targets->setLabel(LabelType::TARGET);
    decoys->setLabel(LabelType::DECOY);
    for (int s = 0; s < n_spectra; ++s)
    {
      for (int h = hits_per_spectrum - 1; h >= 0; --h)
      {
        auto* psm = new PSMDescription();
        psm->features = pool.allocate();
        psm->features[0] = 1000.0 + 100.0 * s - h * 3.7;
        psm->scan = static_cast<unsigned int>(s + 1);
        psm->expMass = 500.0;
        targets->registerPsm(psm);
        in.hit_of[psm->features] = {s, h};
      }
      auto* decoy = new PSMDescription();
      decoy->features = pool.allocate();
      decoy->features[0] = -1.0 - s;
      decoy->scan = static_cast<unsigned int>(s + 1);
      decoy->expMass = 500.0;
      decoys->registerPsm(decoy);
    }
    in.handler.push_back_dataset(targets);
    in.handler.push_back_dataset(decoys);
  }
}

START_TEST(PercolatorBestPositive, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] Scores::generatePositiveTrainingSet / generateTrainingSetsOfSortedSet with trainBestPositive: bestByScore))
{
  // ProSE's rescoring trains on one PSM per spectrum (train_best_positive). Percolator keeps the first PSM of each
  // spectrum that an unstable sort by spectrum leaves, which with many hits per spectrum is mostly not the best one;
  // bestByScore keeps the best-scoring one.
  const int n_spectra = 40, hits = 10;
  Input in;
  fill(in, n_spectra, hits);
  Scores scores(false);
  scores.populateWithPSMs(in.handler);
  std::vector<double> w = {1.0, 0.0};  // score = the feature
  scores.calcScoresAndQvals(w, 0.5);

  for (bool best_by_score : {true, false})
  {
    for (bool sorted_path : {false, true})
    {
      Scores copy = scores;
      AlgIn data(static_cast<unsigned int>(n_spectra * (hits + 1)), 1);
      if (sorted_path)
      {
        TEST_TRUE(copy.generateTrainingSetsOfSortedSet(data, 0.5, true, best_by_score))
      }
      else
      {
        copy.generateNegativeTrainingSet(data, 1.0);
        copy.generatePositiveTrainingSet(data, 0.5, 1.0, true, best_by_score);
      }
      TEST_EQUAL(data.negatives, n_spectra)
      TEST_EQUAL(data.positives, n_spectra)  // one per spectrum (all targets are above every decoy)
      std::set<int> spectra;
      int best = 0;
      for (int i = data.negatives; i < data.negatives + data.positives; ++i)
      {
        TEST_EQUAL(data.Y[i], 1.0)
        const auto it = in.hit_of.find(data.vals[i]);
        TEST_TRUE(it != in.hit_of.end())
        if (it == in.hit_of.end()) continue;
        spectra.insert(it->second.first);
        best += it->second.second == 0;
      }
      TEST_EQUAL(spectra.size(), static_cast<std::size_t>(n_spectra))
      if (best_by_score)
      {
        TEST_EQUAL(best, n_spectra)  // every positive is the best hit of its spectrum
      }
      else
      {
        TEST_TRUE(best < n_spectra)  // Percolator's choice (documents the behaviour bestByScore changes)
      }
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
