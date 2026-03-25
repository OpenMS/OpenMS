// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Michal Startek $
// $Authors: Michal Startek $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmWNet.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureGroupingAlgorithmWNet, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmWNet* ptr = nullptr;
FeatureGroupingAlgorithmWNet* nullPointer = nullptr;

START_SECTION((FeatureGroupingAlgorithmWNet()))
  ptr = new FeatureGroupingAlgorithmWNet();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmWNet()))
  delete ptr;
END_SECTION

START_SECTION((virtual void group(const std::vector<FeatureMap>& maps, ConsensusMap& out)))
{
  FeatureGroupingAlgorithmWNet algo;

  // test exception on too few maps
  vector<FeatureMap> empty_maps;
  ConsensusMap out;
  TEST_EXCEPTION(Exception::IllegalArgument, algo.group(empty_maps, out));

  vector<FeatureMap> one_map(1);
  TEST_EXCEPTION(Exception::IllegalArgument, algo.group(one_map, out));

  // create two simple maps with known matching features
  // Map 1: features at (mz=500, rt=100), (mz=600, rt=200), (mz=700, rt=300)
  // Map 2: features at (mz=500.1, rt=101), (mz=600.1, rt=201), (mz=900, rt=900)
  // Expected: first two pairs should match, third features unmatched
  vector<FeatureMap> maps(2);

  Feature f;
  f.setUniqueId(1);
  f.setMZ(500.0);
  f.setRT(100.0);
  f.setIntensity(1000.0f);
  maps[0].push_back(f);

  f.setUniqueId(2);
  f.setMZ(600.0);
  f.setRT(200.0);
  f.setIntensity(2000.0f);
  maps[0].push_back(f);

  f.setUniqueId(3);
  f.setMZ(700.0);
  f.setRT(300.0);
  f.setIntensity(1500.0f);
  maps[0].push_back(f);

  f.setUniqueId(4);
  f.setMZ(500.1);
  f.setRT(101.0);
  f.setIntensity(1100.0f);
  maps[1].push_back(f);

  f.setUniqueId(5);
  f.setMZ(600.1);
  f.setRT(201.0);
  f.setIntensity(1900.0f);
  maps[1].push_back(f);

  f.setUniqueId(6);
  f.setMZ(900.0);
  f.setRT(900.0);
  f.setIntensity(500.0f);
  maps[1].push_back(f);

  for (auto& m : maps) m.updateRanges();

  // default params: max_distance=100, trash_cost=100, LINF
  // auto mz_scale ≈ 800/400 = 2, so close features (mz diff 0.1*2=0.2, rt diff 1)
  // are well within threshold; distant feature (mz=900, rt=900) is unmatched

  ConsensusMap result;
  algo.group(maps, result);

  // We should get exactly 2 consensus features (the two close pairs)
  TEST_EQUAL(result.size(), 2)

  // Each consensus feature should have exactly 2 handles (one from each map)
  for (Size i = 0; i < result.size(); ++i)
  {
    TEST_EQUAL(result[i].size(), 2)
  }
}
END_SECTION

START_SECTION((virtual void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out)))
  // ConsensusMap grouping follows the same code path via template
  NOT_TESTABLE;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
