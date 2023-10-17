// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureGroupingAlgorithmKD, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmKD* ptr = nullptr;
FeatureGroupingAlgorithmKD* nullPointer = nullptr;
START_SECTION((FeatureGroupingAlgorithmKD()))
  ptr = new FeatureGroupingAlgorithmKD();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmKD()))
  delete ptr;
END_SECTION

START_SECTION((static FeatureGroupingAlgorithm* create()))
  FeatureGroupingAlgorithm* ptr2 = nullptr;
  FeatureGroupingAlgorithm* base_NullPointer = nullptr;
  ptr2 = FeatureGroupingAlgorithmKD::create();
  TEST_NOT_EQUAL(ptr2, base_NullPointer)
  delete ptr2;
END_SECTION

START_SECTION((static String getProductName()))
  TEST_EQUAL(FeatureGroupingAlgorithmKD::getProductName(), "unlabeled_kd")
END_SECTION

START_SECTION((virtual void group(const std::vector<FeatureMap>& maps, ConsensusMap& out)))
  // This is tested in the tool
  NOT_TESTABLE;
END_SECTION

START_SECTION((virtual void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out)))
  // This is tested in the tool
  NOT_TESTABLE;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
