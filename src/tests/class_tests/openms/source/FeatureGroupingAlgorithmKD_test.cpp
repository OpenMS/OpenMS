// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

START_SECTION((virtual void group(const std::vector<FeatureMap>& maps, ConsensusMap& out)))
  // This is tested in the tool
  NOT_TESTABLE;
END_SECTION

START_SECTION((virtual void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out)))
  // This is tested in the tool
  NOT_TESTABLE;
END_SECTION

START_SECTION((group() with two empty FeatureMaps does not crash))
{
  // Regression test: two input maps with no features yield an empty mass range.
  // Previously group_() dereferenced the empty mass range (front()/back()/size()-1)
  // -> crash/UB. After the fix an early guard runs postprocess_() and returns.
  FeatureGroupingAlgorithmKD algo;
  std::vector<FeatureMap> maps(2); // two feature-empty maps
  // metadata carried by the (feature-empty) maps must NOT be silently dropped by the
  // empty-input early return: postprocess_() still transfers protein / unassigned IDs.
  ProteinIdentification prot;
  prot.setIdentifier("run0");
  maps[0].getProteinIdentifications().push_back(prot);
  PeptideIdentification upep;
  upep.setIdentifier("run0");
  maps[0].getUnassignedPeptideIdentifications().push_back(upep);

  ConsensusMap out;
  algo.group(maps, out);
  TEST_EQUAL(out.size(), 0)
  TEST_EQUAL(out.getProteinIdentifications().size(), 1)
  TEST_EQUAL(out.getUnassignedPeptideIdentifications().size(), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
