// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmUnlabeled.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(FeatureGroupingAlgorithmUnlabeled, "$Id FeatureGroupingAlgorithmUnlabeled_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmUnlabeled* ptr = nullptr;
FeatureGroupingAlgorithmUnlabeled* nullPointer = nullptr;
START_SECTION((FeatureGroupingAlgorithmUnlabeled()))
	ptr = new FeatureGroupingAlgorithmUnlabeled();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmUnlabeled()))
	delete ptr;
END_SECTION

START_SECTION((static FeatureGroupingAlgorithm* create()))
  FeatureGroupingAlgorithm* ptr2 = nullptr;
  FeatureGroupingAlgorithm* base_NullPointer = nullptr;
	ptr2 = FeatureGroupingAlgorithmUnlabeled::create();
  TEST_NOT_EQUAL(ptr2, base_NullPointer)
  delete ptr2;
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(FeatureGroupingAlgorithmUnlabeled::getProductName(),"unlabeled")
END_SECTION

START_SECTION((virtual void group(const std::vector< FeatureMap > &maps, ConsensusMap &out)))
	// This is tested extensively in TEST/TOPP
	NOT_TESTABLE;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



