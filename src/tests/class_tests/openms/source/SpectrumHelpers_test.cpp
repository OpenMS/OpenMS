// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest, Witold Wolski $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SpectrumHelpers.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/MockObjects.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;

START_TEST(DiaPrescore2, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION ( [EXTRA] testscorefunction)
{
  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2(new OpenSwath::BinaryDataArray);

  static const double arr1[] = {

    10, 20, 50, 100, 50, 20, 10, // peak at 499
    3, 7, 15, 30, 15, 7, 3,      // peak at 500
    1, 3, 9, 15, 9, 3, 1,        // peak at 501
    3, 9, 3,                     // peak at 502

    10, 20, 50, 100, 50, 20, 10, // peak at 600
    3, 7, 15, 30, 15, 7, 3,      // peak at 601
    1, 3, 9, 15, 9, 3, 1,        // peak at 602
    3, 9, 3                      // peak at 603
  };
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = {
    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03,
    499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03,
    500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
    501.99, 502.0, 502.01,
    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
    600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  data1->data = mz;
  data2->data = intensity;
  sptr->setMZArray( data1);
  sptr->setIntensityArray( data2);

  std::vector<OpenSwath::SpectrumPtr> sptrArr;
  sptrArr.push_back(sptr);

  double mzres, intensityres, imres;
  // mz range from 499 to 501
  RangeMZ mz_range(500.);
  RangeMobility im_range_empty;
  mz_range.minSpanIfSingular(2.);
  DIAHelpers::integrateWindow(sptrArr, mzres, imres, intensityres, mz_range, im_range_empty);

  TEST_REAL_SIMILAR(mzres, 499.392014652015);
  TEST_REAL_SIMILAR(intensityres, 273 );


  // >> exp = [240, 74, 39, 15, 0] > 121 / 500.338842975207
  // >> theo = [1, 0.325757771553019, 0.0678711748364005, 0.0105918703087134, 0.00134955223787482]
  // >> from scipy.stats.stats import pearsonr
  // >> pearsonr(exp, theo)
  // (0.99463189043051314, 0.00047175434098498532)
  //
  mz_range.setMin(499.6);
  mz_range.setMax(501.4);
  DIAHelpers::integrateWindow(sptrArr, mzres, imres, intensityres, mz_range, im_range_empty);

  TEST_REAL_SIMILAR(mzres, 500.338842975207);
  TEST_REAL_SIMILAR(intensityres,121 );

  std::vector<double> wincenter, mzresv, intresv, imresv;
  wincenter.push_back(300.);
  wincenter.push_back(200.);
  wincenter.push_back(500.);
  wincenter.push_back(600.);
  DIAHelpers::integrateWindows(sptrArr,wincenter,0.5, intresv, mzresv, imresv, im_range_empty);
  TEST_REAL_SIMILAR(mzresv[0], 300);
  TEST_REAL_SIMILAR(intresv[0],0 );
  TEST_REAL_SIMILAR(mzresv[1],200 );
  TEST_REAL_SIMILAR(intresv[1],0 );
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

