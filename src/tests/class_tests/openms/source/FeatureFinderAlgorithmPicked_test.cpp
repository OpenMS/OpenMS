// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/Constants.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>
///////////////////////////

#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <utility>

START_TEST(FeatureFinderAlgorithmPicked, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

typedef FeatureFinderAlgorithmPicked FFPP;

FFPP* ptr = nullptr;
FFPP* nullPointer = nullptr;

START_SECTION((FeatureFinderAlgorithmPicked()))
  ptr = new FFPP;
  TEST_NOT_EQUAL(ptr,nullPointer)
END_SECTION

START_SECTION((~FeatureFinderAlgorithmPicked()))
  delete ptr;
END_SECTION

START_SECTION((virtual void run()))
  //input and output
  PeakMap input;
  MzMLFile mzml_file;
  mzml_file.getOptions().addMSLevel(1);
  mzml_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.mzML"),input);
  input.updateRanges();
  FeatureMap output;

  //parameters
  Param param;
  ParamXMLFile paramFile;
  paramFile.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.ini"), param);
  param = param.copy("FeatureFinder:1:algorithm:", true);

  FFPP ffpp;
  ffpp.run(std::move(input), output, param, FeatureMap());

  TEST_EQUAL(output.size(), 8);

  // test some of the metavalue number_of_datapoints
  TEST_EQUAL(output[0].getMetaValue(Constants::UserParam::NUM_OF_DATAPOINTS), 88);
  TEST_EQUAL(output[3].getMetaValue(Constants::UserParam::NUM_OF_DATAPOINTS), 71);
  TEST_EQUAL(output[7].getMetaValue(Constants::UserParam::NUM_OF_DATAPOINTS), 47);

  TOLERANCE_ABSOLUTE(0.001);
  TEST_REAL_SIMILAR(output[0].getOverallQuality(), 0.8826);
  TEST_REAL_SIMILAR(output[1].getOverallQuality(), 0.8680);
  TEST_REAL_SIMILAR(output[2].getOverallQuality(), 0.9077);
  TEST_REAL_SIMILAR(output[3].getOverallQuality(), 0.9270);
  TEST_REAL_SIMILAR(output[4].getOverallQuality(), 0.9398);
  TEST_REAL_SIMILAR(output[5].getOverallQuality(), 0.9098);
  TEST_REAL_SIMILAR(output[6].getOverallQuality(), 0.9403);
  TEST_REAL_SIMILAR(output[7].getOverallQuality(), 0.9245);

  TOLERANCE_ABSOLUTE(20.0);
  TEST_REAL_SIMILAR(output[0].getIntensity(), 51366.2);
  TEST_REAL_SIMILAR(output[1].getIntensity(), 44767.6);
  TEST_REAL_SIMILAR(output[2].getIntensity(), 34731.1);
  TEST_REAL_SIMILAR(output[3].getIntensity(), 19494.2);
  TEST_REAL_SIMILAR(output[4].getIntensity(), 12570.2);
  TEST_REAL_SIMILAR(output[5].getIntensity(), 8532.26);
  TEST_REAL_SIMILAR(output[6].getIntensity(), 7318.62);
  TEST_REAL_SIMILAR(output[7].getIntensity(), 5038.81);

END_SECTION

START_SECTION(([EXTRA] isotopic_pattern:mz_tolerance and mass_trace:mz_tolerance are not interchangeable (#9247)))
{
  // PR #9247 fixed a swap in updateMembers_(): pattern_tolerance_ <- isotopic_pattern:mz_tolerance and
  // trace_tolerance_ <- mass_trace:mz_tolerance. Pin both asymmetric configurations directionally, so a
  // repeated assignment swap exchanges the outcomes and fails the test.
  MzMLFile mzml_file;
  mzml_file.getOptions().addMSLevel(1);
  PeakMap input_template;
  mzml_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.mzML"), input_template);
  input_template.updateRanges();

  Param base;
  ParamXMLFile paramFile;
  paramFile.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.ini"), base);
  base = base.copy("FeatureFinder:1:algorithm:", true);

  // tight isotope-pattern tolerance, loose mass-trace tolerance
  Param p1 = base;
  p1.setValue("isotopic_pattern:mz_tolerance", 0.005);
  p1.setValue("mass_trace:mz_tolerance", 0.5);
  FeatureMap out1;
  {
    PeakMap in = input_template;
    FFPP ff;
    ff.run(std::move(in), out1, p1, FeatureMap());
  }

  // swapped assignment: loose isotope-pattern tolerance, tight mass-trace tolerance
  Param p2 = base;
  p2.setValue("isotopic_pattern:mz_tolerance", 0.5);
  p2.setValue("mass_trace:mz_tolerance", 0.005);
  FeatureMap out2;
  {
    PeakMap in = input_template;
    FFPP ff;
    ff.run(std::move(in), out2, p2, FeatureMap());
  }

  TEST_EQUAL(out1.size(), 1)
  TEST_EQUAL(out2.size(), 0)

  if (out1.size() == 1)
  {
    TEST_EQUAL(out1[0].getMetaValue(Constants::UserParam::NUM_OF_DATAPOINTS), 33)

    TOLERANCE_ABSOLUTE(0.001);
    TEST_REAL_SIMILAR(out1[0].getRT(), 4278.1601)
    TEST_REAL_SIMILAR(out1[0].getMZ(), 653.7722)
    TEST_REAL_SIMILAR(out1[0].getOverallQuality(), 0.9609)

    TOLERANCE_ABSOLUTE(20.0);
    TEST_REAL_SIMILAR(out1[0].getIntensity(), 18467.8)
  }
}
END_SECTION

START_SECTION(([EXTRA] isotopic_pattern:charge_low above charge_high is rejected))
{
  // charge_high - charge_low + 1 was stored in a UInt and wrapped, so the score arrays were written out of bounds.
  FFPP ffpp;
  Param p = ffpp.getParameters();
  p.setValue("isotopic_pattern:charge_low", 4);
  p.setValue("isotopic_pattern:charge_high", 2);
  TEST_EXCEPTION(Exception::InvalidParameter, ffpp.setParameters(p))
  p.setValue("isotopic_pattern:charge_low", 3);
  TEST_EXCEPTION(Exception::InvalidParameter, ffpp.setParameters(p))
  p.setValue("isotopic_pattern:charge_low", 2);
  ffpp.setParameters(p); // a single charge is fine
  TEST_EQUAL(int(ffpp.getParameters().getValue("isotopic_pattern:charge_low")), 2)
}
END_SECTION

START_SECTION(([EXTRA] feature:min_isotope_fit 0 aborts seeds without an isotope pattern))
{
  // findBestIsotopeFit_ returns 0 exactly when it finds no placement (it only stores a pattern with a score > 0),
  // and that empty pattern was then extended (out-of-bounds read). Such seeds must be aborted, so a bound of 0 has to
  // give the same features as a bound just above 0.
  MzMLFile mzml_file;
  mzml_file.getOptions().addMSLevel(1);
  PeakMap input_template;
  mzml_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.mzML"), input_template);
  input_template.updateRanges();
  Param base;
  ParamXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.ini"), base);
  base = base.copy("FeatureFinder:1:algorithm:", true);
  base.setValue("feature:reported_mz", "average");
  base.setValue("feature:min_trace_score", 0.0);
  base.setValue("feature:min_score", 0.0);
  base.setValue("seed:min_score", 0.0);

  FeatureMap out_zero, out_tiny;
  {
    Param p = base;
    p.setValue("feature:min_isotope_fit", 0.0);
    PeakMap in = input_template;
    FFPP ff;
    ff.run(std::move(in), out_zero, p, FeatureMap());
  }
  {
    Param p = base;
    p.setValue("feature:min_isotope_fit", 1e-300);
    PeakMap in = input_template;
    FFPP ff;
    ff.run(std::move(in), out_tiny, p, FeatureMap());
  }
  TEST_NOT_EQUAL(out_tiny.size(), 0)
  TEST_EQUAL(out_zero.size(), out_tiny.size())
  for (Size i = 0; i < std::min(out_zero.size(), out_tiny.size()); ++i)
  {
    TEST_REAL_SIMILAR(out_zero[i].getRT(), out_tiny[i].getRT())
    TEST_REAL_SIMILAR(out_zero[i].getMZ(), out_tiny[i].getMZ())
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
