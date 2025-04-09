// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/SYSTEM/File.h>

///////////////////////////

#include <OpenMS/FORMAT/AbsoluteQuantitationMethodFile.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

using namespace OpenMS;
using namespace std;

class AbsoluteQuantitationMethodFile_facade : AbsoluteQuantitationMethodFile
{
  public:
    void parseLine_(StringList & line, std::map<String,Size> & headers, AbsoluteQuantitationMethod & aqm)
    {
      AbsoluteQuantitationMethodFile::parseLine_(line, headers, aqm);
    }
};

///////////////////////////

START_TEST(AbsoluteQuantitationMethodFile, "$Id$")

/////////////////////////////////////////////////////////////

AbsoluteQuantitationMethodFile* ptr = nullptr;
AbsoluteQuantitationMethodFile* nullPointer = nullptr;
const String in_file_1 = OPENMS_GET_TEST_DATA_PATH("AbsoluteQuantitationMethodFile_in_1.csv");
const String in_file_2 = OPENMS_GET_TEST_DATA_PATH("AbsoluteQuantitationMethodFile_in_2.csv");
const String out_file = File::getTemporaryFile();

START_SECTION((AbsoluteQuantitationMethodFile()))
	ptr = new AbsoluteQuantitationMethodFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~AbsoluteQuantitationMethodFile()))
	delete ptr;
END_SECTION

START_SECTION(void parseLine_(StringList & line, std::map<String,Size> & headers, AbsoluteQuantitationMethod & aqm) const)
  AbsoluteQuantitationMethodFile_facade aqmf;
  AbsoluteQuantitationMethod aqm;

  // headers
  std::map<String, Size> headers;
  headers["IS_name"] = 0;
  headers["component_name"] = 1;
  headers["feature_name"] = 2;
  headers["concentration_units"] = 3;
  headers["llod"] = 4;
  headers["ulod"] = 5;
  headers["lloq"] = 6;
  headers["uloq"] = 7;
  headers["correlation_coefficient"] = 8;
  headers["n_points"] = 9;
  headers["transformation_model"] = 10;
  headers["transformation_model_param_slope"] = 11;
  headers["transformation_model_param_intercept"] = 12;

  // line test 1
  StringList line1;
  line1.push_back("IS1");
  line1.push_back("component1");
  line1.push_back("feature1");
  line1.push_back("uM");
  line1.push_back("3.0");
  line1.push_back("  "); //test for empty string
  line1.push_back(" 2.0  "); //test for leading and trailing white spaces
  line1.push_back("8.0");
  line1.push_back("0.99");
  line1.push_back("5");
  line1.push_back("TransformationModelLinear");
  line1.push_back("2.0");
  line1.push_back("1.0");

  aqmf.parseLine_(line1, headers, aqm);

  TEST_EQUAL(aqm.getISName(), "IS1");
  TEST_EQUAL(aqm.getComponentName(), "component1");
  TEST_EQUAL(aqm.getFeatureName(), "feature1");
  TEST_EQUAL(aqm.getConcentrationUnits(), "uM");
  TEST_REAL_SIMILAR(aqm.getLLOD(), 3.0);
  TEST_REAL_SIMILAR(aqm.getULOD(), 0.0);
  TEST_REAL_SIMILAR(aqm.getLLOQ(), 2.0);
  TEST_REAL_SIMILAR(aqm.getULOQ(), 8.0);
  TEST_REAL_SIMILAR(aqm.getCorrelationCoefficient(), 0.99);
  TEST_EQUAL(aqm.getNPoints(), 5);
  TEST_EQUAL(aqm.getTransformationModel(), "TransformationModelLinear");
  const Param transformation_model_params = aqm.getTransformationModelParams();
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"), 2.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"), 1.0);
END_SECTION

START_SECTION(void load(const String & filename, std::vector<AbsoluteQuantitationMethod> & aqm_list))
  AbsoluteQuantitationMethodFile aqmf;
  std::vector<AbsoluteQuantitationMethod> aqm_list;

  aqmf.load(in_file_1, aqm_list);
  TEST_EQUAL(aqm_list[0].getComponentName(), "component1");
  TEST_EQUAL(aqm_list[0].getISName(), "IS1");
  TEST_EQUAL(aqm_list[0].getFeatureName(), "feature1");
  TEST_EQUAL(aqm_list[0].getConcentrationUnits(), "uM");
  TEST_REAL_SIMILAR(aqm_list[0].getLLOD(), 0.0);
  TEST_REAL_SIMILAR(aqm_list[0].getULOD(), 10.0);
  TEST_REAL_SIMILAR(aqm_list[0].getLLOQ(), 2.0);
  TEST_REAL_SIMILAR(aqm_list[0].getULOQ(), 8.0);
  TEST_REAL_SIMILAR(aqm_list[0].getCorrelationCoefficient(), 0.99);
  TEST_EQUAL(aqm_list[0].getNPoints(), 5);
  TEST_EQUAL(aqm_list[0].getTransformationModel(), "TransformationModelLinear");
  Param transformation_model_params;
  transformation_model_params = aqm_list[0].getTransformationModelParams();
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"), 2.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"), 1.0);

  TEST_EQUAL(aqm_list[1].getComponentName(), "component2");
  TEST_EQUAL(aqm_list[1].getISName(), "IS2");
  TEST_EQUAL(aqm_list[1].getFeatureName(), "feature2");
  TEST_EQUAL(aqm_list[1].getConcentrationUnits(), "uM");
  TEST_REAL_SIMILAR(aqm_list[1].getLLOD(), 1.0);
  TEST_REAL_SIMILAR(aqm_list[1].getULOD(), 9.0);
  TEST_REAL_SIMILAR(aqm_list[1].getLLOQ(), 3.0);
  TEST_REAL_SIMILAR(aqm_list[1].getULOQ(), 7.0);
  TEST_REAL_SIMILAR(aqm_list[1].getCorrelationCoefficient(), 0.98);
  TEST_EQUAL(aqm_list[1].getNPoints(), 4);
  TEST_EQUAL(aqm_list[1].getTransformationModel(), "TransformationModelLinear");
  transformation_model_params = aqm_list[1].getTransformationModelParams();
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"), 2.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"), 2.0);

  TEST_EQUAL(aqm_list[2].getComponentName(), "component3");
  TEST_EQUAL(aqm_list[2].getISName(), "IS3");
  TEST_EQUAL(aqm_list[2].getFeatureName(), "feature3");
  TEST_EQUAL(aqm_list[2].getConcentrationUnits(), "uM");
  TEST_REAL_SIMILAR(aqm_list[2].getLLOD(), 2.0);
  TEST_REAL_SIMILAR(aqm_list[2].getULOD(), 8.0);
  TEST_REAL_SIMILAR(aqm_list[2].getLLOQ(), 4.0);
  TEST_REAL_SIMILAR(aqm_list[2].getULOQ(), 6.0);
  TEST_REAL_SIMILAR(aqm_list[2].getCorrelationCoefficient(), 0.97);
  TEST_EQUAL(aqm_list[2].getNPoints(), 3);
  TEST_EQUAL(aqm_list[2].getTransformationModel(), "TransformationModelLinear");
  transformation_model_params = aqm_list[2].getTransformationModelParams();
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"), 1.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"), 2.0);

  TEST_EQUAL(aqm_list[6].getComponentName(), "component 7"); // Checking for space within the name
  TEST_EQUAL(aqm_list[6].getISName(), ""); // empty cell, default value is used.
  TEST_EQUAL(aqm_list[6].getFeatureName(), "feature 7");
  TEST_EQUAL(aqm_list[6].getConcentrationUnits(), "");
  TEST_REAL_SIMILAR(aqm_list[6].getLLOD(), 0.0); // empty cell, default value is used.
  TEST_REAL_SIMILAR(aqm_list[6].getULOD(), 0.0);
  TEST_REAL_SIMILAR(aqm_list[6].getLLOQ(), 0.0);
  TEST_REAL_SIMILAR(aqm_list[6].getULOQ(), 0.0);
  TEST_REAL_SIMILAR(aqm_list[6].getCorrelationCoefficient(), 0.0);
  TEST_EQUAL(aqm_list[6].getNPoints(), 0);
  TEST_EQUAL(aqm_list[6].getTransformationModel(), "");
  transformation_model_params = aqm_list[6].getTransformationModelParams();
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"), 0.0); // empty cell, default value is used.
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"), 2.0);

  TEST_EQUAL(aqm_list[7].getComponentName(), "component8");
  TEST_EQUAL(aqm_list[7].getISName(), "IS8");
  TEST_EQUAL(aqm_list[7].getFeatureName(), "feature8");
  TEST_EQUAL(aqm_list[7].getConcentrationUnits(), "uM");
  TEST_REAL_SIMILAR(aqm_list[7].getLLOD(), 7.0);
  TEST_REAL_SIMILAR(aqm_list[7].getULOD(), 3.0);
  TEST_REAL_SIMILAR(aqm_list[7].getLLOQ(), 0.0);
  TEST_REAL_SIMILAR(aqm_list[7].getULOQ(), 1.0);
  TEST_REAL_SIMILAR(aqm_list[7].getCorrelationCoefficient(), 0.92);
  TEST_EQUAL(aqm_list[7].getNPoints(), 1);
  TEST_EQUAL(aqm_list[7].getTransformationModel(), "TransformationModelLinear");
  transformation_model_params = aqm_list[7].getTransformationModelParams();
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"), 1.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"), 2.0);

  // The following input file doesn't have the headers: component_name, llod
  // Note that a default value of "" and 0 is given for these missing columns.
  aqmf.load(in_file_2, aqm_list);
  TEST_EQUAL(aqm_list[0].getComponentName(), ""); // A component name with a default value.
  TEST_EQUAL(aqm_list[0].getISName(), "IS1");
  TEST_EQUAL(aqm_list[0].getFeatureName(), "feature1");
  TEST_EQUAL(aqm_list[0].getConcentrationUnits(), "uM");
  TEST_REAL_SIMILAR(aqm_list[0].getLLOD(), 0.0); // A LLOD with a default value.
  TEST_REAL_SIMILAR(aqm_list[0].getULOD(), 10.0);
  TEST_REAL_SIMILAR(aqm_list[0].getLLOQ(), 2.0);
  TEST_REAL_SIMILAR(aqm_list[0].getULOQ(), 8.0);
  TEST_REAL_SIMILAR(aqm_list[0].getCorrelationCoefficient(), 0.99);
  TEST_EQUAL(aqm_list[0].getNPoints(), 5);
  TEST_EQUAL(aqm_list[0].getTransformationModel(), "TransformationModelLinear");
  transformation_model_params = aqm_list[0].getTransformationModelParams();
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"), 2.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"), 1.0);

  TEST_EQUAL(aqm_list[1].getComponentName(), ""); // A component name with a default value.
  TEST_EQUAL(aqm_list[1].getISName(), ""); // empty cell, default value is used.
  TEST_EQUAL(aqm_list[1].getFeatureName(), "feature 7");
  TEST_EQUAL(aqm_list[1].getConcentrationUnits(), "");
  TEST_REAL_SIMILAR(aqm_list[1].getLLOD(), 0.0); // empty cell, default value is used.
  TEST_REAL_SIMILAR(aqm_list[1].getULOD(), 0.0);
  TEST_REAL_SIMILAR(aqm_list[1].getLLOQ(), 0.0);
  TEST_REAL_SIMILAR(aqm_list[1].getULOQ(), 0.0);
  TEST_REAL_SIMILAR(aqm_list[1].getCorrelationCoefficient(), 0.0);
  TEST_EQUAL(aqm_list[1].getNPoints(), 0);
  TEST_EQUAL(aqm_list[1].getTransformationModel(), "");
  transformation_model_params = aqm_list[1].getTransformationModelParams();
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"), 0.0); // empty cell, default value is used.
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"), 2.0);

  TEST_EQUAL(aqm_list[2].getComponentName(), ""); // A component name with a default value.
  TEST_EQUAL(aqm_list[2].getISName(), "IS8");
  TEST_EQUAL(aqm_list[2].getFeatureName(), "feature8");
  TEST_EQUAL(aqm_list[2].getConcentrationUnits(), "uM");
  TEST_REAL_SIMILAR(aqm_list[2].getLLOD(), 0.0); // A LLOD with a default value.
  TEST_REAL_SIMILAR(aqm_list[2].getULOD(), 3.0);
  TEST_REAL_SIMILAR(aqm_list[2].getLLOQ(), 0.0);
  TEST_REAL_SIMILAR(aqm_list[2].getULOQ(), 1.0);
  TEST_REAL_SIMILAR(aqm_list[2].getCorrelationCoefficient(), 0.92);
  TEST_EQUAL(aqm_list[2].getNPoints(), 1);
  TEST_EQUAL(aqm_list[2].getTransformationModel(), "TransformationModelLinear");
  transformation_model_params = aqm_list[2].getTransformationModelParams();
  TEST_REAL_SIMILAR(transformation_model_params.getValue("slope"), 1.0);
  TEST_REAL_SIMILAR(transformation_model_params.getValue("intercept"), 2.0);
END_SECTION

START_SECTION(void store(const String & filename, const std::vector<AbsoluteQuantitationMethod> & aqm_list) const)
  AbsoluteQuantitationMethodFile aqmf;
  vector<AbsoluteQuantitationMethod> aqm_list1, aqm_list2;
  aqmf.load(in_file_1, aqm_list1);
  aqmf.store(out_file, aqm_list1);
  aqmf.load(out_file, aqm_list2);
  TEST_EQUAL(aqm_list1.size(), aqm_list2.size())
  for (Size i = 0; i < aqm_list1.size(); ++i)
  {
    TEST_EQUAL(aqm_list1[i] == aqm_list2[i], true)
  }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
