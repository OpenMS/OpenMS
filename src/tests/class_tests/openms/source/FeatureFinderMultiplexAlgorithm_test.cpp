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
  exp.updateRanges();
  
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

  // Check that the per-mass-trace meta values (XICs of the individual isotopes) are attached to the features.
  FeatureMap feature_map = algorithm.getFeatureMap();
  TEST_NOT_EQUAL(feature_map.size(), 0);
  for (const Feature& f : feature_map)
  {
    TEST_EQUAL(f.metaValueExists("masstrace_intensity"), true);
    TEST_EQUAL(f.metaValueExists("masstrace_centroid_rt"), true);
    TEST_EQUAL(f.metaValueExists("masstrace_centroid_mz"), true);
    TEST_EQUAL(f.metaValueExists("num_of_masstraces"), true);

    std::vector<double> intensities = f.getMetaValue("masstrace_intensity");
    std::vector<double> centroid_rt = f.getMetaValue("masstrace_centroid_rt");
    std::vector<double> centroid_mz = f.getMetaValue("masstrace_centroid_mz");
    Size num_traces = (Size)f.getMetaValue("num_of_masstraces");

    // one entry per isotope mass trace, sizes consistent across the three lists
    TEST_NOT_EQUAL(num_traces, 0);
    TEST_EQUAL(intensities.size(), num_traces);
    TEST_EQUAL(centroid_rt.size(), num_traces);
    TEST_EQUAL(centroid_mz.size(), num_traces);

    // mass trace intensities are non-negative integrated areas
    for (double intensity : intensities)
    {
      TEST_EQUAL(intensity >= 0.0, true);
    }
  }
}
END_SECTION

END_TEST
