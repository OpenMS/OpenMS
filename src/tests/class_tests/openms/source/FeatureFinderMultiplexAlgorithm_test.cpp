// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFinderMultiplexAlgorithm, "$Id$")

FeatureFinderMultiplexAlgorithm* ptr = 0;
FeatureFinderMultiplexAlgorithm* null_ptr = 0;
START_SECTION(FeatureFinderMultiplexAlgorithm())
{
  ptr = new FeatureFinderMultiplexAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION(~FeatureFinderMultiplexAlgorithm())
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual void run()))
{
  MzMLFile mzml_file;
  MSExperiment exp;
  ConsensusMap result;
  
  mzml_file.getOptions().addMSLevel(1);
  mzml_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderMultiplex_1_input.mzML"), exp);
  exp.updateRanges(1);
  
  Param param;
  ParamXMLFile paramFile;
  paramFile.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderMultiplex_1_parameters.ini"), param);
  param = param.copy("FeatureFinderMultiplex:1:",true);
  param.remove("in");
  param.remove("out");
  param.remove("out_multiplets");
  param.remove("log");
  param.remove("debug");
  param.remove("threads");
  param.remove("no_progress");
  param.remove("force");
  param.remove("test");
  
  FeatureFinderMultiplexAlgorithm algorithm;
  algorithm.setParameters(param);
  algorithm.run(exp, true);
  result = algorithm.getConsensusMap();
  
  TEST_EQUAL(result.size(), 2);
  
  double L = result[0].getFeatures().begin()->getIntensity();
  double H = (++(result[0].getFeatures().begin()))->getIntensity();

  // Check that the HEAVY:LIGHT ratio is close to the expected 3:1 ratio
  TOLERANCE_ABSOLUTE(0.2);
  TEST_REAL_SIMILAR(H/L, 3.0);
}
END_SECTION

END_TEST
