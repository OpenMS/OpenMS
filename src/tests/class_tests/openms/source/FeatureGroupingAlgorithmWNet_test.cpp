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

  // Use Da mode: mz diff 0.1 Da is within 0.3 Da threshold,
  // rt diff 1 s is within 100 s threshold.
  // Distant feature (mz=900, rt=900) exceeds both thresholds.
  Param p = algo.getParameters();
  p.setValue("mz_unit", "Da");
  algo.setParameters(p);

  ConsensusMap result;
  algo.group(maps, result);

  // We should get 4 consensus features: 2 matched pairs + 2 singletons
  // (mz=700/rt=300 from map0 and mz=900/rt=900 from map1 are unmatched)
  TEST_EQUAL(result.size(), 4)

  // Count matched pairs (size 2) and singletons (size 1)
  Size pairs = 0, singletons = 0;
  for (Size i = 0; i < result.size(); ++i)
  {
    if (result[i].size() == 2) ++pairs;
    else if (result[i].size() == 1) ++singletons;
  }
  TEST_EQUAL(pairs, 2)
  TEST_EQUAL(singletons, 2)

  // ppm mode: same maps, same expected outcome
  {
    FeatureGroupingAlgorithmWNet algo_ppm;
    Param pp = algo_ppm.getParameters();
    pp.setValue("mz_unit", "ppm");
    pp.setValue("max_mz_shift_ppm", 500.0); // 500 ppm @ 500 Da ≈ 0.25 Da
    algo_ppm.setParameters(pp);

    ConsensusMap result_ppm;
    algo_ppm.group(maps, result_ppm);

    TEST_EQUAL(result_ppm.size(), 4)

    Size pairs_ppm = 0, singletons_ppm = 0;
    for (Size i = 0; i < result_ppm.size(); ++i)
    {
      if (result_ppm[i].size() == 2) ++pairs_ppm;
      else if (result_ppm[i].size() == 1) ++singletons_ppm;
    }
    TEST_EQUAL(pairs_ppm, 2)
    TEST_EQUAL(singletons_ppm, 2)
  }

  // decharge_mz mode: features at z=2 should be projected to [M+H]+ before matching
  {
    vector<FeatureMap> dc_maps(2);

    Feature fa;
    fa.setUniqueId(11);
    fa.setMZ(250.5); // z=2: neutral = 2*250.5 - 2*1.00728 ≈ 498.985; [M+H]+ ≈ 499.993
    fa.setRT(100.0);
    fa.setIntensity(1000.0f);
    fa.setCharge(2);
    dc_maps[0].push_back(fa);

    Feature fb;
    fb.setUniqueId(12);
    fb.setMZ(500.0); // z=1: already singly charged
    fb.setRT(100.0);
    fb.setIntensity(900.0f);
    fb.setCharge(1);
    dc_maps[1].push_back(fb);

    for (auto& m : dc_maps) m.updateRanges();

    FeatureGroupingAlgorithmWNet algo_dc;
    Param pdc = algo_dc.getParameters();
    pdc.setValue("mz_unit", "Da");
    pdc.setValue("decharge_mz", "true");
    algo_dc.setParameters(pdc);

    ConsensusMap result_dc;
    algo_dc.group(dc_maps, result_dc);

    // After decharging, fa maps to ~499.993 Da and fb is 500.0 Da,
    // difference ~0.007 Da, well within 0.3 Da threshold.
    Size pairs_dc = 0;
    for (Size i = 0; i < result_dc.size(); ++i)
    {
      if (result_dc[i].size() == 2) ++pairs_dc;
    }
    TEST_EQUAL(pairs_dc, 1)
  }
}
END_SECTION

START_SECTION((virtual void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out)))
{
  FeatureGroupingAlgorithmWNet algo;
  Param p = algo.getParameters();
  p.setValue("mz_unit", "Da");
  algo.setParameters(p);

  // Build two ConsensusMap inputs with matching ConsensusFeatures
  std::vector<ConsensusMap> maps(2);

  ConsensusFeature cf1;
  cf1.setUniqueId(101);
  cf1.setMZ(400.0);
  cf1.setRT(50.0);
  cf1.setIntensity(500.0f);
  maps[0].push_back(cf1);

  ConsensusFeature cf2;
  cf2.setUniqueId(102);
  cf2.setMZ(500.0);
  cf2.setRT(150.0);
  cf2.setIntensity(800.0f);
  maps[0].push_back(cf2);

  ConsensusFeature cf3;
  cf3.setUniqueId(103);
  cf3.setMZ(400.1);
  cf3.setRT(51.0);
  cf3.setIntensity(520.0f);
  maps[1].push_back(cf3);

  ConsensusFeature cf4;
  cf4.setUniqueId(104);
  cf4.setMZ(900.0);
  cf4.setRT(900.0);
  cf4.setIntensity(100.0f);
  maps[1].push_back(cf4);

  for (auto& m : maps) m.updateRanges();

  ConsensusMap result;
  algo.group(maps, result);

  // cf1/cf3 should match (Δmz=0.1 Da within 0.3 Da, Δrt=1 s within 100 s)
  // cf2 and cf4 should remain as singletons
  TEST_EQUAL(result.size(), 3)

  Size pairs = 0, singletons = 0;
  for (Size i = 0; i < result.size(); ++i)
  {
    if (result[i].size() == 2) ++pairs;
    else if (result[i].size() == 1) ++singletons;
  }
  TEST_EQUAL(pairs, 1)
  TEST_EQUAL(singletons, 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
