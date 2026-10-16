// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Veit $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmKD.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmKD.h>
#include <OpenMS/KERNEL/Feature.h>

using namespace OpenMS;
using namespace std;

class MapAlignmentAlgorithmKDTest : public MapAlignmentAlgorithmKD
{
public:
  using MapAlignmentAlgorithmKD::MapAlignmentAlgorithmKD;
  using MapAlignmentAlgorithmKD::filterCCs_;
};

START_TEST(MapAlignmentAlgorithmKD, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MapAlignmentAlgorithmKD* ptr = nullptr;
MapAlignmentAlgorithmKD* nullPointer = nullptr;

START_SECTION((MapAlignmentAlgorithmKD(Size num_maps, const Param& param)))
  ptr = new MapAlignmentAlgorithmKD(42, Param());
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MapAlignmentAlgorithmKD()))
  delete ptr;
END_SECTION

START_SECTION((void addRTFitData(const KDTreeFeatureMaps& kd_data)))
  NOT_TESTABLE;
END_SECTION

START_SECTION((void fitLOWESS()))
  NOT_TESTABLE;
END_SECTION

START_SECTION((void transform(KDTreeFeatureMaps& kd_data) const))
  NOT_TESTABLE;
END_SECTION

START_SECTION([EXTRA] filterCCs_ rejects connected components with mixed charge states)
{
  // two input maps: min CC size is max(2, 0.5 * 2) = 2, no conflicts allowed
  MapAlignmentAlgorithmKDTest aligner(2, FeatureGroupingAlgorithmKD().getParameters());

  vector<Feature> features(6);
  // CC 0: same charge in both maps -> kept
  features[0].setCharge(2);
  features[1].setCharge(2);
  // CC 1: charges 2 and 3 -> rejected
  features[2].setCharge(2);
  features[3].setCharge(3);
  // CC 2: an unknown charge (0) does not count as a conflict -> kept
  features[4].setCharge(2);
  features[5].setCharge(0);

  KDTreeFeatureMaps kd_data;
  for (Size i = 0; i < features.size(); ++i)
  {
    features[i].setRT(100.0 * i);
    features[i].setMZ(500.0 + i);
    kd_data.addFeature(i % 2, &features[i]); // alternate maps, so no CC has a map conflict
  }

  map<Size, vector<Size>> ccs;
  ccs[0] = {0, 1};
  ccs[1] = {2, 3};
  ccs[2] = {4, 5};

  map<Size, vector<Size>> filtered_ccs;
  aligner.filterCCs_(kd_data, ccs, filtered_ccs);

  TEST_EQUAL(filtered_ccs.size(), 2)
  TEST_EQUAL(filtered_ccs.count(0), 1)
  TEST_EQUAL(filtered_ccs.count(1), 0)
  TEST_EQUAL(filtered_ccs.count(2), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
