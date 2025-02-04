// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include "OpenMS/OPENSWATHALGO/ALGO/StatsHelpers.h"

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;


START_TEST(DiaPrescore2, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION ( testscorefunction)
{
  static const double arr1[] = {
    10, 20, 50, 100, 50, 20, 10, // peak at 600
    3, 7, 15, 30, 15, 7, 3,      // peak at 601
    1, 3, 9, 15, 9, 3, 1,        // peak at 602
    3, 9, 3                      // peak at 603
  };
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(double) );
  static const double arr2[] = {
    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
    600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(double) );
  double norm = OpenSwath::norm(mz.begin(),mz.end());
  std::vector<double> normalized;
  OpenSwath::normalize(mz,norm,normalized);
  TEST_REAL_SIMILAR(OpenSwath::norm(normalized.begin(),normalized.end()), 1.);
  double x = dotProd(normalized.begin(),normalized.end(),normalized.begin());
  TEST_REAL_SIMILAR(x, 1.);
  double man = manhattanDist(normalized.begin(),normalized.end(),normalized.begin());
  TEST_REAL_SIMILAR(man, 0.);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

