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
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

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
  MzDataFile mzdata_file;
  mzdata_file.getOptions().addMSLevel(1);
  mzdata_file.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.mzData"),input);
  input.updateRanges(1);
  FeatureMap output;

  //parameters
  Param param;
  ParamXMLFile paramFile;
  paramFile.load(OPENMS_GET_TEST_DATA_PATH("FeatureFinderAlgorithmPicked.ini"), param);
  param = param.copy("FeatureFinder:1:algorithm:", true);

  FFPP ffpp;
  ffpp.run(input, output, param, FeatureMap());

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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
