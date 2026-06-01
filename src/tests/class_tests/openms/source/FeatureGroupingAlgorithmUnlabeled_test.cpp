// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

START_SECTION((virtual void group(const std::vector< FeatureMap > &maps, ConsensusMap &out)))
	// This is tested extensively in TEST/TOPP
	NOT_TESTABLE;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



