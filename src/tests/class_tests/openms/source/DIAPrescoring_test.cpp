// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAPrescoring.h>
#include "OpenMS/OPENSWATHALGO/DATAACCESS/MockObjects.h"
#include <OpenMS/KERNEL/RangeManager.h>

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;


START_TEST(DiaPrescore2, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

const std::string ION_MOBILITY_DESCRIPTION = "Ion Mobility";

DiaPrescore* ptr = nullptr;
DiaPrescore* nullPointer = nullptr;

START_SECTION(DiaPrescore())
{
  ptr = new DiaPrescore();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DiaPrescore())
{
  delete ptr;
}
END_SECTION

START_SECTION ( test score function with perfect first transition and ion mobility filtering )
{
  OpenSwath::LightTransition mock_tr1;
  mock_tr1.product_mz = 500.;
  mock_tr1.fragment_charge = 1;
  mock_tr1.transition_name = "group1";
  mock_tr1.library_intensity = 5.;

  OpenSwath::LightTransition mock_tr2;
  mock_tr2.product_mz = 600.;
  mock_tr2.fragment_charge = 1;
  mock_tr2.transition_name = "group2";
  mock_tr2.library_intensity = 5.;

  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2(new OpenSwath::BinaryDataArray);
  OpenMS::RangeMobility im_range_empty;

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
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(double) );
  static const double arr2[] = {
      499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03,
      500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
      501.97, 501.98, 501.99, 502.0, 502.01, 502.02, 502.03,
      502.99, 503.0, 503.01,

      599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
      600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
      601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
      602.99, 603.0, 603.01,
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(double) );

  data1->data = mz;
  data2->data = intensity;
  sptr->setMZArray( data1);
  sptr->setIntensityArray( data2);

  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr1);
  transitions.push_back(mock_tr2);

  DiaPrescore diaprescore(0.05);
  double manhattan = 0., dotprod = 0.;

  std::vector <OpenSwath::SpectrumPtr> sptrArr;
  sptrArr.push_back(sptr);

  diaprescore.score(sptrArr, transitions , im_range_empty, dotprod, manhattan);
  //std::cout << "dotprod : " << dotprod << std::endl;
  //std::cout << "manhattan : " << manhattan << std::endl;
  // >> exp = [240, 74, 39, 15, 0]
  // >> theo = [1, 0.325757771553019, 0.0678711748364005, 0.0105918703087134, 0.00134955223787482]
  // >> from scipy.stats.stats import pearsonr
  // >> pearsonr(exp, theo)
  // (0.99463189043051314, 0.00047175434098498532)
  //
  TEST_REAL_SIMILAR(dotprod, 0.916131286812994)
  TEST_REAL_SIMILAR(manhattan, 0.23670593984202)
}
END_SECTION

START_SECTION ( test score function missing first transition )
{
  OpenSwath::LightTransition mock_tr1;
  mock_tr1.product_mz = 500.;
  mock_tr1.fragment_charge = 1;
  mock_tr1.transition_name = "group1";
  mock_tr1.library_intensity = 5.;

  OpenSwath::LightTransition mock_tr2;
  mock_tr2.product_mz = 600.;
  mock_tr2.fragment_charge = 1;
  mock_tr2.transition_name = "group2";
  mock_tr2.library_intensity = 5.;

  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);

  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2(new OpenSwath::BinaryDataArray);

  static const double arr1[] = {
    /*
    10, 20, 50, 100, 50, 20, 10, // peak at 499
    3, 7, 15, 30, 15, 7, 3,      // peak at 500
    1, 3, 9, 15, 9, 3, 1,        // peak at 501
    3, 9, 3,                     // peak at 502
    */

    10, 20, 50, 100, 50, 20, 10, // peak at 600
    3, 7, 15, 30, 15, 7, 3,      // peak at 601
    1, 3, 9, 15, 9, 3, 1,        // peak at 602
    3, 9, 3                      // peak at 603
  };
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(double) );
  static const double arr2[] = {
    /*
    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03,
    499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03,
    500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
    501.99, 502.0, 502.01,
    */
    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
    600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(double) );
  data1->data = mz;
  data2->data = intensity;
  sptr->setMZArray( data1);
  sptr->setIntensityArray( data2);

  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr1);
  transitions.push_back(mock_tr2);

  DiaPrescore diaprescore(0.05);
  OpenMS::RangeMobility im_range_empty;
  double manhattan = 0., dotprod = 0.;

  std::vector <OpenSwath::SpectrumPtr> sptrArr;
  sptrArr.push_back(sptr);
  diaprescore.score(sptrArr, transitions, im_range_empty, dotprod, manhattan);
  //std::cout << "dotprod : " << dotprod << std::endl;
  //std::cout << "manhattan : " << manhattan << std::endl;
  // >> exp = [240, 74, 39, 15, 0]
  // >> theo = [1, 0.325757771553019, 0.0678711748364005, 0.0105918703087134, 0.00134955223787482]
  // >> from scipy.stats.stats import pearsonr
  // >> pearsonr(exp, theo)
  // (0.99463189043051314, 0.00047175434098498532)
  //
  TEST_REAL_SIMILAR(dotprod, 0.627263258948172)
  TEST_REAL_SIMILAR(manhattan, 0.984211129641047)
}
END_SECTION

START_SECTION ( test score function with shifted first transition )
{
  OpenSwath::LightTransition mock_tr1;
  mock_tr1.product_mz = 500.;
  mock_tr1.fragment_charge = 1;
  mock_tr1.transition_name = "group1";
  mock_tr1.library_intensity = 5.;

  OpenSwath::LightTransition mock_tr2;
  mock_tr2.product_mz = 600.;
  mock_tr2.fragment_charge = 1;
  mock_tr2.transition_name = "group2";
  mock_tr2.library_intensity = 5.;

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
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(double) );
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
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(double) );
  data1->data = mz;
  data2->data = intensity;
  sptr->setMZArray( data1);
  sptr->setIntensityArray( data2);

  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr1);
  transitions.push_back(mock_tr2);

  DiaPrescore diaprescore(0.05);
  OpenMS::RangeMobility im_range_empty;
  double manhattan = 0., dotprod = 0.;

  std::vector <OpenSwath::SpectrumPtr> sptrArr;
  sptrArr.push_back(sptr);

  diaprescore.score(sptrArr, transitions, im_range_empty, dotprod, manhattan);
  //std::cout << "dotprod : " << dotprod << std::endl;
  //std::cout << "manhattan : " << manhattan << std::endl;
  // >> exp = [240, 74, 39, 15, 0]
  // >> theo = [1, 0.325757771553019, 0.0678711748364005, 0.0105918703087134, 0.00134955223787482]
  // >> from scipy.stats.stats import pearsonr
  // >> pearsonr(exp, theo)
  // (0.99463189043051314, 0.00047175434098498532)
  //
  TEST_REAL_SIMILAR(dotprod, 0.43738312515644)
  TEST_REAL_SIMILAR(manhattan, 0.557433222328531)
}
END_SECTION

START_SECTION ( test score function missing first transition due to different ion mobility )
{
  OpenSwath::LightTransition mock_tr1;
  mock_tr1.product_mz = 500.;
  mock_tr1.fragment_charge = 1;
  mock_tr1.transition_name = "group1";
  mock_tr1.library_intensity = 5.;

  OpenSwath::LightTransition mock_tr2;
  mock_tr2.product_mz = 600.;
  mock_tr2.fragment_charge = 1;
  mock_tr2.transition_name = "group2";
  mock_tr2.library_intensity = 5.;

  double PRECURSOR_ION_MOBILITY = 7;
  double ION_MOBILITY_WIDTH = 2;

  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data3(new OpenSwath::BinaryDataArray);

  std::vector<double> intensity {

    10, 20, 50, 100, 50, 20, 10, // peak at 499
    3, 7, 15, 30, 15, 7, 3,      // peak at 500
    1, 3, 9, 15, 9, 3, 1,        // peak at 501
    3, 9, 3,                     // peak at 502

    10, 20, 50, 100, 50, 20, 10, // peak at 600
    3, 7, 15, 30, 15, 7, 3,      // peak at 601
    1, 3, 9, 15, 9, 3, 1,        // peak at 602
    3, 9, 3                      // peak at 603
  };
  std::vector<double> mz {

    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03,
    499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03,
    500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
    501.99, 502.0, 502.01,

    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
    600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
  };

  std::vector<double> im {
      1, 1, 3, 1, 1, 1, 1,      // peak at 499
      2, 2, 3, 1, 2, 1, 2,      // peak at 500
      1, 2, 1, 2, 3, 2, 2,      // peak at 501
      2, 2, 2,                  // peak at 502

      7, 6, 7, 6, 8, 6, 8,      // peak at 600
      7, 7, 7, 8, 7, 7, 8,      // peak at 601
      8, 6, 8, 8, 6, 6, 6,      // peak at 602
      6, 8, 6                   // peak at 603
  };

  data1->data = mz;
  data2->data = intensity;
  data3->data = im;
  sptr->setMZArray(data1);
  sptr->setIntensityArray(data2);
  data3->description = ION_MOBILITY_DESCRIPTION;
  sptr->getDataArrays().push_back(data3);

  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr1);
  transitions.push_back(mock_tr2);

  DiaPrescore diaprescore(0.05);
  OpenMS::RangeMobility im_range(PRECURSOR_ION_MOBILITY);
  im_range.minSpanIfSingular(ION_MOBILITY_WIDTH);
  double manhattan = 0., dotprod = 0.;

  std::vector <OpenSwath::SpectrumPtr> sptrArr;
  sptrArr.push_back(sptr);
  diaprescore.score(sptrArr, transitions, im_range, dotprod, manhattan);
  //std::cout << "dotprod : " << dotprod << std::endl;
  //std::cout << "manhattan : " << manhattan << std::endl;
  // >> exp = [240, 74, 39, 15, 0]
  // >> theo = [1, 0.325757771553019, 0.0678711748364005, 0.0105918703087134, 0.00134955223787482]
  // >> from scipy.stats.stats import pearsonr
  // >> pearsonr(exp, theo)
  // (0.99463189043051314, 0.00047175434098498532)
  //
  TEST_REAL_SIMILAR(dotprod, 0.627263258948172)
  TEST_REAL_SIMILAR(manhattan, 0.984211129641047)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
