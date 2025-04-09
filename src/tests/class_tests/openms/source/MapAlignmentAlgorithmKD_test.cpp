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

using namespace OpenMS;
using namespace std;

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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
