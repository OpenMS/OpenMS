// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(FeatureGroupingAlgorithmQT, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmQT* ptr = nullptr;
FeatureGroupingAlgorithmQT* nullPointer = nullptr;
START_SECTION((FeatureGroupingAlgorithmQT()))
	ptr = new FeatureGroupingAlgorithmQT();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmQT()))
	delete ptr;
END_SECTION

START_SECTION((static FeatureGroupingAlgorithm* create()))
	FeatureGroupingAlgorithm* ptr2 = nullptr;
  FeatureGroupingAlgorithm* base_NullPointer = nullptr;
	ptr2 = FeatureGroupingAlgorithmQT::create();
  TEST_NOT_EQUAL(ptr2, base_NullPointer)
  delete ptr2;
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(FeatureGroupingAlgorithmQT::getProductName(), "unlabeled_qt")
END_SECTION

START_SECTION((virtual void group(const std::vector< FeatureMap >& maps, ConsensusMap& out)))
	// This is tested extensively in TEST/TOPP
	NOT_TESTABLE;
END_SECTION

START_SECTION((virtual void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out)))
	// This is tested extensively in TEST/TOPP
	NOT_TESTABLE;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



