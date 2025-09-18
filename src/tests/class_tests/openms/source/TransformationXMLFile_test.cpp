// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/TransformationXMLFile.h>

///////////////////////////

START_TEST(FASTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationXMLFile* ptr = nullptr;
TransformationXMLFile* nullPointer = nullptr;

START_SECTION((TransformationXMLFile()))
{
  ptr = new TransformationXMLFile();
  TEST_NOT_EQUAL(ptr, nullPointer);
  delete ptr;
}
END_SECTION

START_SECTION([EXTRA] static bool isValid(const String& filename))
{
  TransformationXMLFile f;
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML"), std::cerr), true);
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_2.trafoXML"), std::cerr), true);
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_3.trafoXML"), std::cerr), false);
  TEST_EQUAL(f.isValid(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_4.trafoXML"), std::cerr), true);
}
END_SECTION


START_SECTION(void load(const String & filename, TransformationDescription & transformation, bool fit_model=true))
{
  TransformationDescription trafo;
  TransformationXMLFile trafo_xml;
  Param params;

  trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_1.trafoXML"), trafo);
  TEST_STRING_EQUAL(trafo.getModelType(), "none");
  params = trafo.getModelParameters();
  TEST_EQUAL(params.empty(), true);

  trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_2.trafoXML"), trafo);
  TEST_STRING_EQUAL(trafo.getModelType(), "linear");
  params = trafo.getModelParameters();
  TEST_EQUAL(params.size(), 2);
  TEST_REAL_SIMILAR(params.getValue("slope"), 3.141592653589793238);
  TEST_REAL_SIMILAR(params.getValue("intercept"), 2.718281828459045235);

  trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_4.trafoXML"), trafo);
  TEST_STRING_EQUAL(trafo.getModelType(), "interpolated");
  params = trafo.getModelParameters();
  TEST_EQUAL(params.getValue("interpolation_type"), "linear");
  TEST_EQUAL(trafo.getDataPoints().size(), 3);
  TEST_REAL_SIMILAR(trafo.getDataPoints()[0].first, 1.2);
  TEST_REAL_SIMILAR(trafo.getDataPoints()[1].first, 2.2);
  TEST_REAL_SIMILAR(trafo.getDataPoints()[2].first, 3.2);
  TEST_REAL_SIMILAR(trafo.getDataPoints()[0].second, 5.2);
  TEST_REAL_SIMILAR(trafo.getDataPoints()[1].second, 6.25);
  TEST_REAL_SIMILAR(trafo.getDataPoints()[2].second, 7.3);

  // also test the option of not performing the actual model fit
  trafo_xml.load(OPENMS_GET_TEST_DATA_PATH("TransformationXMLFile_2.trafoXML"), trafo, false);
  TEST_STRING_EQUAL(trafo.getModelType(), "none");
  params = trafo.getModelParameters();
  TEST_EQUAL(params.empty(), true);
}
END_SECTION

START_SECTION(void store(String filename, const TransformationDescription& transformation))
{
  TransformationDescription trafo, trafo2;
  TransformationXMLFile trafo_xml;

  String tmp_file_none;
  Param params;
  trafo.fitModel("none", params);
  NEW_TMP_FILE(tmp_file_none);
  trafo_xml.store(tmp_file_none, trafo);
  trafo_xml.load(tmp_file_none, trafo2);
  TEST_STRING_EQUAL(trafo2.getModelType(), "none");
  params = trafo2.getModelParameters();
  TEST_EQUAL(params.empty(), true);
  {
    double pre_image = 234255132.43212;
    double image = trafo.apply(pre_image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }
  String tmp_file_linear;
  NEW_TMP_FILE(tmp_file_linear);
  params.setValue("slope", 3.141592653589793238);
  params.setValue("intercept", 2.718281828459045235);
  trafo.fitModel("linear", params);
  trafo_xml.store(tmp_file_linear, trafo);
  trafo_xml.load(tmp_file_linear, trafo2);
  TEST_STRING_EQUAL(trafo.getModelType(), "linear");
  params.clear();
  params = trafo2.getModelParameters();
  TEST_EQUAL(params.size(), 2);
  TEST_REAL_SIMILAR(params.getValue("slope"), 3.141592653589793238);
  TEST_REAL_SIMILAR(params.getValue("intercept"), 2.718281828459045235);
  {
    double pre_image = 234255132.43212;
    double image = trafo.apply(pre_image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }

  String tmp_file_pairs;
  NEW_TMP_FILE(tmp_file_pairs);
  TransformationDescription::DataPoints pairs;
  pairs.push_back(make_pair(1.2, 5.2));
  pairs.push_back(make_pair(2.2, 6.25));
  pairs.push_back(make_pair(3.2, 7.3));
  trafo.setDataPoints(pairs);
  params.clear();
  params.setValue("interpolation_type", "linear");
  trafo.fitModel("interpolated", params);
  trafo_xml.store(tmp_file_pairs, trafo);
  trafo_xml.load(tmp_file_pairs, trafo2);
  TEST_STRING_EQUAL(trafo2.getModelType(), "interpolated");
  params = trafo2.getModelParameters();
  TEST_EQUAL(params.size(), 2);
  TEST_STRING_EQUAL(params.getValue("interpolation_type"), "linear");
  TEST_STRING_EQUAL(params.getValue("extrapolation_type"), "two-point-linear");
  TEST_EQUAL(trafo2.getDataPoints().size(), 3);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].first, 1.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[1].first, 2.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[2].first, 3.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].second, 5.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[1].second, 6.25);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[2].second, 7.3);
  {
    double pre_image = 234255132.43212;
    double image = trafo.apply(pre_image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }

  TEST_EXCEPTION(Exception::IllegalArgument, trafo.fitModel("mumble_pfrwoarpfz"));
#if 0
  String tmp_file_bspline;
  NEW_TMP_FILE(tmp_file_bspline);
  pairs.clear();
  pairs.push_back(make_pair(1.2, 5.2));
  pairs.push_back(make_pair(3.2, 7.3));
  pairs.push_back(make_pair(2.2, 6.25));
  pairs.push_back(make_pair(2.2, 3.1));
  pairs.push_back(make_pair(2.2, 7.25));
  pairs.push_back(make_pair(3.0, 8.5));
  pairs.push_back(make_pair(3.1, 4.7));
  pairs.push_back(make_pair(1.7, 6.0));
  pairs.push_back(make_pair(2.9, 4.7));
  pairs.push_back(make_pair(4.2, 5.0));
  pairs.push_back(make_pair(3.7, -2.4));
  trafo.setDataPoints(pairs);
  params.clear();
  params.setValue("num_breakpoints", 4);
  params.setValue("break_positions", "uniform");
  trafo.fitModel("b_spline", params);
  trafo_xml.store(tmp_file_pairs, trafo);
  trafo_xml.load(tmp_file_pairs, trafo2);
  TEST_STRING_EQUAL(trafo2.getModelType(), "b_spline");
  params.clear();
  params = trafo2.getModelParameters();
  TEST_EQUAL(params.getValue("num_breakpoints"), 4);
  TEST_EQUAL(params.getValue("break_positions"), "uniform");
  TEST_EQUAL(params.size(), 8);
  TEST_EQUAL(trafo2.getDataPoints().size(), 11);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].first, 1.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[0].second, 5.2);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[10].first, 3.7);
  TEST_REAL_SIMILAR(trafo2.getDataPoints()[10].second, -2.4);
  for (Int breaks = 0; breaks < 10; ++breaks)
  {
    if (breaks == 1) continue;
    params.setValue("num_breakpoints", breaks);
    trafo.fitModel("b_spline", params);
    STATUS("breaks: " << breaks);
    double pre_image = 234255132.43212;
    double image = trafo.apply(pre_image);
    STATUS("Here is an invocation of trafo.apply():   pre_image: " << pre_image << "  image: " << image);
  }
#endif
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
