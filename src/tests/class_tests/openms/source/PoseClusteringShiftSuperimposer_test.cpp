// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h>
///////////////////////////

#include <OpenMS/KERNEL/Feature.h>

using namespace OpenMS;
using namespace std;

typedef DPosition <2> PositionType;

START_TEST(PoseClusteringShiftSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PoseClusteringShiftSuperimposer* ptr = nullptr;
PoseClusteringShiftSuperimposer* nullPointer = nullptr;

START_SECTION((PoseClusteringShiftSuperimposer()))
	ptr = new PoseClusteringShiftSuperimposer();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PoseClusteringShiftSuperimposer()))
	delete ptr;
END_SECTION

START_SECTION((virtual void run(const ConsensusMap& map_model, const ConsensusMap& map_scene, TransformationDescription& transformation)))

  std::vector<ConsensusMap> input(2);

  Feature feat1;
  Feature feat2;
  PositionType pos1(1,1);
  PositionType pos2(5,5);
  feat1.setPosition(pos1);
  feat1.setIntensity(100.0f);
  feat2.setPosition(pos2);
  feat2.setIntensity(100.0f);
  input[0].push_back(ConsensusFeature(feat1));
  input[0].push_back(ConsensusFeature(feat2));

  Feature feat3;
  Feature feat4;
  PositionType pos3(21.4,1.02);
  PositionType pos4(25.4,5.02);
  feat3.setPosition(pos3);
  feat3.setIntensity(100.0f);
  feat4.setPosition(pos4);
  feat4.setIntensity(100.0f);
  input[1].push_back(ConsensusFeature(feat3));
  input[1].push_back(ConsensusFeature(feat4));

  TransformationDescription transformation;
  PoseClusteringShiftSuperimposer pcat;
	Param params;
#if 0 // switch this on for debugging
  params.setValue("dump_buckets","tmp_PoseClusteringShiftSuperimposer_buckets");
  params.setValue("dump_pairs","tmp_PoseClusteringShiftSuperimposer_pairs");
  pcat.setParameters(params);
#endif
	pcat.run(input[0], input[1], transformation);

  TEST_STRING_EQUAL(transformation.getModelType(), "linear")
	params = transformation.getModelParameters();
	TEST_EQUAL(params.size(), 2)
  TEST_REAL_SIMILAR(params.getValue("slope"), 1.0)
  TEST_REAL_SIMILAR(params.getValue("intercept"), -20.4)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



