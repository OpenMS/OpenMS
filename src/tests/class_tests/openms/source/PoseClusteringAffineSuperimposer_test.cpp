// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
///////////////////////////

#include <OpenMS/KERNEL/Feature.h>

using namespace OpenMS;
using namespace std;

typedef DPosition <2> PositionType;

START_TEST(PoseClusteringAffineSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PoseClusteringAffineSuperimposer* ptr = nullptr;
PoseClusteringAffineSuperimposer* nullPointer = nullptr;

START_SECTION((PoseClusteringAffineSuperimposer()))
{
  ptr = new PoseClusteringAffineSuperimposer();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((virtual ~PoseClusteringAffineSuperimposer()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual void run(const ConsensusMap& map_model, const ConsensusMap& map_scene, TransformationDescription& transformation)))
{
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
  PositionType pos3(1.4,1.02);
  PositionType pos4(5.4,5.02);
  feat3.setPosition(pos3);
  feat3.setIntensity(100.0f);
  feat4.setPosition(pos4);
  feat4.setIntensity(100.0f);
  input[1].push_back(ConsensusFeature(feat3));
  input[1].push_back(ConsensusFeature(feat4));

  Param parameters;
  parameters.setValue(String("scaling_bucket_size"), 0.01);
  parameters.setValue(String("shift_bucket_size"), 0.1);

  // If hashing goes wrong, get debug output with the following:
  //  parameters.setValue(String("dump_buckets"),"pcast_buckets");
  //  parameters.setValue(String("dump_pairs"),"pcast_pairs");

  TransformationDescription transformation;
  PoseClusteringAffineSuperimposer pcat;
  pcat.setParameters(parameters);

  // That's a precondition for run()!  Now even documented :-)
  input[0].updateRanges();
  input[1].updateRanges();

  pcat.run(input[0], input[1], transformation);

  TEST_STRING_EQUAL(transformation.getModelType(), "linear")
  parameters = transformation.getModelParameters();
  TEST_EQUAL(parameters.size(), 2)
  TEST_REAL_SIMILAR(parameters.getValue("slope"), 1.0)
  TEST_REAL_SIMILAR(parameters.getValue("intercept"), -0.4)
}
END_SECTION

START_SECTION((virtual void run(const std::vector<Peak2D> & map_model, const std::vector<Peak2D> & map_scene, TransformationDescription& transformation)))
{
  std::vector<Peak2D> map_model, map_scene;

  Peak2D p1;
  p1.setRT(1);
  p1.setMZ(1);
  p1.setIntensity(100.0f);
  Peak2D p2;
  p2.setRT(5);
  p2.setMZ(5);
  p2.setIntensity(100.0f);
  map_model.push_back(p1);
  map_model.push_back(p2);

  Peak2D p3;
  p3.setRT(1.4);
  p3.setMZ(1.02);
  p3.setIntensity(100.0f);
  Peak2D p4;
  p4.setRT(5.4);
  p4.setMZ(5.02);
  p4.setIntensity(100.0f);
  map_scene.push_back(p3);
  map_scene.push_back(p4);

  Param parameters;
  parameters.setValue(String("scaling_bucket_size"), 0.01);
  parameters.setValue(String("shift_bucket_size"), 0.1);

  // If hashing goes wrong, get debug output with the following:
  //  parameters.setValue(String("dump_buckets"),"pcast_buckets");
  //  parameters.setValue(String("dump_pairs"),"pcast_pairs");

  TransformationDescription transformation;
  PoseClusteringAffineSuperimposer pcat;
  pcat.setParameters(parameters);

  pcat.run(map_model, map_scene, transformation);

  TEST_STRING_EQUAL(transformation.getModelType(), "linear")
  parameters = transformation.getModelParameters();
  TEST_EQUAL(parameters.size(), 2)
  TEST_REAL_SIMILAR(parameters.getValue("slope"), 1.0)
  TEST_REAL_SIMILAR(parameters.getValue("intercept"), -0.4)
}
END_SECTION

START_SECTION(([EXTRA]virtual void run(const std::vector<Peak2D> & map_model, const std::vector<Peak2D> & map_scene, TransformationDescription& transformation)))
{
  std::vector<Peak2D> map_model, map_scene;

  // map1_rt = 
  double map1_rt[] = {1.0, 5.0};
  double map2_rt[] = {1.4, 5.4};

  double map1_mz[] = {1.0 , 5.0 };
  double map2_mz[] = {1.02, 5.02};

  double map1_int[] = {100, 100};
  double map2_int[] = {100, 100};

  for (Size i = 0; i < 2; i++)
  {
    Peak2D p;
    p.setRT(map1_rt[i]);
    p.setMZ(map1_mz[i]);
    p.setIntensity(map1_int[i]);
    map_model.push_back(p);
  }
  for (Size i = 0; i < 2; i++)
  {
    Peak2D p;
    p.setRT(map2_rt[i]);
    p.setMZ(map2_mz[i]);
    p.setIntensity(map2_int[i]);
    map_scene.push_back(p);
  }


  Param parameters;
  parameters.setValue(String("scaling_bucket_size"), 0.01);
  parameters.setValue(String("shift_bucket_size"), 0.1);

  // If hashing goes wrong, get debug output with the following:
  //  parameters.setValue(String("dump_buckets"),"pcast_buckets");
  //  parameters.setValue(String("dump_pairs"),"pcast_pairs");

  TransformationDescription transformation;
  PoseClusteringAffineSuperimposer pcat;
  pcat.setParameters(parameters);

  pcat.run(map_model, map_scene, transformation);

  TEST_STRING_EQUAL(transformation.getModelType(), "linear")
  parameters = transformation.getModelParameters();
  TEST_EQUAL(parameters.size(), 2)
  TEST_REAL_SIMILAR(parameters.getValue("slope"), 1.0)
  TEST_REAL_SIMILAR(parameters.getValue("intercept"), -0.4)
}
END_SECTION

START_SECTION(([EXTRA]virtual void run(const std::vector<Peak2D> & map_model, const std::vector<Peak2D> & map_scene, TransformationDescription& transformation)))
{
  std::vector<Peak2D> map_model, map_scene;

  // add another point at 5.2 -> 5.8 RT (and add some chaff in the middle)
  double map1_rt[] = {1.0, 5.0, 1.3, 2.2, 5.2};
  double map2_rt[] = {1.4, 5.4, 4.4, 4.4, 5.8};

  double map1_mz[] = {1.0 , 5.0 , 800, 900, 5.0 };
  double map2_mz[] = {1.02, 5.02, 800, 900, 5.02};

  double map1_int[] = {100, 100, 41, 20, 50};
  double map2_int[] = {100, 100, 40, 20, 50};

  for (Size i = 0; i < 5; i++)
  {
    Peak2D p;
    p.setRT(map1_rt[i]);
    p.setMZ(map1_mz[i]);
    p.setIntensity(map1_int[i]);
    map_model.push_back(p);
  }
  for (Size i = 0; i < 5; i++)
  {
    Peak2D p;
    p.setRT(map2_rt[i]);
    p.setMZ(map2_mz[i]);
    p.setIntensity(map2_int[i]);
    map_scene.push_back(p);
  }

  // make sure vector is not really sorted
  std::reverse(map_model.begin(), map_model.end() );
  std::reverse(map_scene.begin(), map_scene.end() );

  // using 2 points
  {
    Param parameters;
    parameters.setValue(String("scaling_bucket_size"), 0.01);
    parameters.setValue(String("shift_bucket_size"), 0.1);
    parameters.setValue(String("num_used_points"), 2); // only use first two points -> same results as before expected

    TransformationDescription transformation;
    PoseClusteringAffineSuperimposer pcat;
    pcat.setParameters(parameters);

    pcat.run(map_model, map_scene, transformation);

    TEST_STRING_EQUAL(transformation.getModelType(), "linear")
    parameters = transformation.getModelParameters();
    TEST_EQUAL(parameters.size(), 2)
    TEST_REAL_SIMILAR(parameters.getValue("slope"), 1.0)
    TEST_REAL_SIMILAR(parameters.getValue("intercept"), -0.4)
  }

  // using 3 points
  {
    Param parameters;
    parameters.setValue(String("scaling_bucket_size"), 0.01);
    parameters.setValue(String("shift_bucket_size"), 0.1);
    parameters.setValue(String("num_used_points"), 3); // only use first three points -> different results as before expected

    TransformationDescription transformation;
    PoseClusteringAffineSuperimposer pcat;
    pcat.setParameters(parameters);

    pcat.run(map_model, map_scene, transformation);

    TEST_STRING_EQUAL(transformation.getModelType(), "linear")
    parameters = transformation.getModelParameters();
    TEST_EQUAL(parameters.size(), 2)
    TEST_REAL_SIMILAR(parameters.getValue("slope"), 0.977273) // slope should be less than before
    TEST_REAL_SIMILAR(parameters.getValue("intercept"), -0.368182) // intercept should be higher than before
  }

  // what happens if we set the wrong parameters?
  {
    Param parameters;
    parameters.setValue(String("scaling_bucket_size"), 0.01);
    parameters.setValue(String("shift_bucket_size"), 0.1);
    parameters.setValue(String("num_used_points"), 3); // only use first three points -> different results as before expected
    parameters.setValue(String("max_shift"), 0.2);
    parameters.setValue(String("max_scaling"), 1.001);

    TransformationDescription transformation;
    PoseClusteringAffineSuperimposer pcat;
    pcat.setParameters(parameters);

    pcat.run(map_model, map_scene, transformation);

    // quite easy: we get the wrong results!
    // TODO: don't let this happen, so easy to prevent!
    TEST_STRING_EQUAL(transformation.getModelType(), "linear")
    parameters = transformation.getModelParameters();
    TEST_EQUAL(parameters.size(), 2)
    TEST_REAL_SIMILAR(parameters.getValue("slope"), 1.0) // TODO this is completely wrong
    TEST_REAL_SIMILAR(parameters.getValue("intercept"), -0.4) // TODO this is completely wrong
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


