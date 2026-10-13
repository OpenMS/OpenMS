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
#include "Globals.h"
#include "PosteriorEstimator.h"

#include <cmath>
#include <limits>
#include <utility>
#include <vector>

using namespace OpenMS;
using namespace std;
using OpenMS::Internal::Percolator::Globals;
using OpenMS::Internal::Percolator::PosteriorEstimator;

namespace
{
  // PSMs best score first (as Scores::getScoreLabelPairs returns them): for each
  // score, n_targets targets and n_decoys decoys.
  struct Group { double score; int n_targets; int n_decoys; };

  vector<pair<double, bool>> makePSMs(const vector<Group>& groups)
  {
    vector<pair<double, bool>> combined;
    for (const Group& g : groups)
    {
      for (int i = 0; i < g.n_targets; ++i) combined.emplace_back(g.score, true);
      for (int i = 0; i < g.n_decoys; ++i) combined.emplace_back(g.score, false);
    }
    return combined;
  }

  // PEPs must be probabilities and must not decrease with decreasing score.
  bool validPEPs(const vector<pair<double, bool>>& combined, const vector<double>& peps)
  {
    if (peps.size() != combined.size()) return false;
    for (size_t i = 0; i < peps.size(); ++i)
    {
      if (!std::isfinite(peps[i]) || peps[i] < 0.0 || peps[i] > 1.0) return false;
      if (i > 0 && combined[i].first < combined[i - 1].first && peps[i] < peps[i - 1]) return false;
    }
    return true;
  }
}

START_TEST(PercolatorPosteriorEstimator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// The PEP spline (BaseSpline::initiateQR) needs strictly increasing bin medians;
// tied medians asserted in Debug builds (dx > 0) and divided by zero in Release
// builds. Seen with ProSE -rescore (nonparametric PEP, no-terminate flag set).
Globals::getInstance()->setNoTerminate(true);

START_SECTION([EXTRA] estimatePEP with fewer than 4 bins)
{
  // Fewer than 4 distinct scores give fewer than 4 bins, which are padded to 4.
  // Upstream padded with an offset of at most 1e-20, lost to rounding here.
  for (bool use_pi0 : {false, true})
  {
    for (double offset : {0.0, -10.0, 1000.0})
    {
      vector<pair<double, bool>> combined = makePSMs({{offset + 2.0, 20, 0}, {offset + 1.0, 10, 10}});
      vector<double> peps;
      PosteriorEstimator::estimatePEP(combined, use_pi0, 1.0, peps, true);
      TEST_TRUE(validPEPs(combined, peps))
      TEST_TRUE(peps.front() < peps.back())

      combined = makePSMs({{offset + 2.0, 20, 0}, {offset + 1.0, 10, 10}, {offset, 5, 15}});
      PosteriorEstimator::estimatePEP(combined, use_pi0, 1.0, peps, true);
      TEST_TRUE(validPEPs(combined, peps))
      TEST_TRUE(peps.front() < peps.back())

      // a single score: one bin
      combined = makePSMs({{offset + 2.0, 20, 5}});
      PosteriorEstimator::estimatePEP(combined, use_pi0, 1.0, peps, true);
      TEST_TRUE(validPEPs(combined, peps))
    }
  }
}
END_SECTION

START_SECTION([EXTRA] estimatePEP with bin medians closer than the spline resolves)
{
  // Scores 1 ulp apart are separate bins, but the spline stores their distance
  // (~1e-16) in sparse vectors that treat |x| <= 1e-15 as 0.
  const double near_minus_one = std::nextafter(-1.0, 0.0);
  for (bool use_pi0 : {false, true})
  {
    vector<pair<double, bool>> combined = makePSMs({{2.0, 10, 0}, {1.0, 8, 2}, {0.0, 6, 4},
      {near_minus_one, 4, 6}, {-1.0, 4, 6}, {-2.0, 2, 8}});
    vector<double> peps;
    PosteriorEstimator::estimatePEP(combined, use_pi0, 1.0, peps, true);
    TEST_TRUE(validPEPs(combined, peps))
    TEST_TRUE(peps.front() < peps.back())
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
