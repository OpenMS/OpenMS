// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>

///////////////////////////

START_TEST(TransformationModelLinear, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationModelLinear* ptr = nullptr;
TransformationModelLinear* nullPointer = nullptr;

TransformationModel::DataPoints data, empty;
TransformationModel::DataPoint point;
point.first = 0.0;
point.second = 1.0;
data.push_back(point);
point.first = 1.0;
point.second = 2.0;
data.push_back(point);
point.first = 1.0;
point.second = 4.0;
data.push_back(point);

START_SECTION((TransformationModelLinear(const DataPoints &, const Param &)))
{
  TEST_EXCEPTION(Exception::IllegalArgument, TransformationModelLinear lm(empty, Param())); // need data
  ptr = new TransformationModelLinear(data, Param());
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~TransformationModelLinear()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual double evaluate(double value) const))
{
  ptr = new TransformationModelLinear(data, Param());

  TEST_REAL_SIMILAR(ptr->evaluate(-0.5), 0.0);
  TEST_REAL_SIMILAR(ptr->evaluate(0.0), 1.0);
  TEST_REAL_SIMILAR(ptr->evaluate(0.5), 2.0);
  TEST_REAL_SIMILAR(ptr->evaluate(1.0), 3.0);
  TEST_REAL_SIMILAR(ptr->evaluate(1.5), 4.0);

  delete ptr;
}
END_SECTION

START_SECTION((void getParameters(Param & params) const))
{  

  TransformationModel::DataPoint point;
  point.first = 2.0;
  point.second = 2.0;
  data.push_back(point);
  Param p_in;
  //test weightings
  p_in.setValue("symmetric_regression", "true");
  p_in.setValue("x_weight", "ln(x)");
  p_in.setValue("y_weight", "ln(y)");
  p_in.setValue("x_datum_min", 10e-5);
  p_in.setValue("x_datum_max", 1e15);
  p_in.setValue("y_datum_min", 10e-8);
  p_in.setValue("y_datum_max", 1e15);
  TransformationModelLinear lm0(data, p_in);
  Param p_out = p_in;
  p_out.setValue("slope", 0.095036911971605034);
  p_out.setValue("intercept", 0.89550911545438994);
  TEST_EQUAL(lm0.getParameters(), p_out);

  //add additional data and test without weightings
  p_in.setValue("x_weight", "x");
  p_in.setValue("y_weight", "y");
  p_in.setValue("x_datum_min", 10e-5);
  p_in.setValue("x_datum_max", 1e15);
  p_in.setValue("y_datum_min", 10e-8);
  p_in.setValue("y_datum_max", 1e15);
  TransformationModelLinear lm(data, p_in);
  p_out = p_in;
  p_out.setValue("slope", 0.5);
  p_out.setValue("intercept", 1.75);
  TEST_EQUAL(lm.getParameters(), p_out);

  //test with empty data
  p_in.clear();
  p_in.setValue("slope", 12.3);
  p_in.setValue("intercept", -45.6);
  p_in.setValue("x_weight", "x");
  p_in.setValue("y_weight", "y");
  p_in.setValue("x_datum_min", 10e-5);
  p_in.setValue("x_datum_max", 1e15);
  p_in.setValue("y_datum_min", 10e-8);
  p_in.setValue("y_datum_max", 1e15);
  TransformationModelLinear lm2(empty, p_in);
  TEST_EQUAL(lm2.getParameters(), p_in);
}
END_SECTION

START_SECTION(([EXTRA] void getParameters(double&, double&, String&, String&, double&, double&, double&, double&)))
{
  Param param;
  param.setValue("slope", 12.3);
  param.setValue("intercept", -45.6);  
  String x_weight_test, y_weight_test;
  x_weight_test = "x";
  y_weight_test = "ln(y)";
  param.setValue("x_weight", x_weight_test);
  param.setValue("y_weight", y_weight_test);
  param.setValue("x_datum_min", 1e-15);
  param.setValue("x_datum_max", 1e15);
  param.setValue("y_datum_min", 1e-15);
  param.setValue("y_datum_max", 1e15);
  TransformationModelLinear lm(empty, param);
  double slope, intercept, x_datum_min, x_datum_max, y_datum_min, y_datum_max;
  String x_weight, y_weight;
  lm.getParameters(slope, intercept, x_weight, y_weight, x_datum_min, x_datum_max, y_datum_min, y_datum_max);
  TEST_REAL_SIMILAR(param.getValue("slope"), slope);
  TEST_REAL_SIMILAR(param.getValue("intercept"), intercept);
  TEST_EQUAL(param.getValue("x_weight"), x_weight);
  TEST_EQUAL(param.getValue("y_weight"), y_weight);
  TEST_REAL_SIMILAR(param.getValue("x_datum_min"), x_datum_min);
  TEST_REAL_SIMILAR(param.getValue("x_datum_max"), x_datum_max);
  TEST_REAL_SIMILAR(param.getValue("y_datum_min"), y_datum_min);
  TEST_REAL_SIMILAR(param.getValue("y_datum_max"), y_datum_max);
}
END_SECTION

START_SECTION((TransformationModelLinear(const DataPoints &, const Param &)))
{
  // weighting/unweighting test 1
  // set-up the parameters
  Param param; 
  String x_weight_test, y_weight_test;
  x_weight_test = "ln(x)";
  y_weight_test = "ln(y)";
  param.setValue("x_weight", x_weight_test);
  param.setValue("y_weight", y_weight_test);
  param.setValue("x_datum_min", 1e-15);
  param.setValue("x_datum_max", 1e8);
  param.setValue("y_datum_min", 1e-8);
  param.setValue("y_datum_max", 1e15);

  // set-up the data and test
  TransformationModel::DataPoints data1;
  data1.clear();
  TransformationModel::DataPoint point;
  point.first = 1.0;
  point.second = 2.0;
  data1.push_back(point);
  point.first = 2.0;
  point.second = 4.0;
  data1.push_back(point);
  point.first = 4.0;
  point.second = 8.0;
  data1.push_back(point);

  // test evaluate
  TransformationModelLinear lm(data1, param);
  TEST_REAL_SIMILAR(lm.evaluate(2), 4);

  // test evaluate using the inverted model
  lm.invert();
  TEST_REAL_SIMILAR(lm.evaluate(4), 2);

  // weighting/unweighting test 2
  // set-up the parameters
  x_weight_test = "1/x";
  y_weight_test = "y";
  param.setValue("x_weight", x_weight_test);
  param.setValue("y_weight", y_weight_test);

  // test evaluate
  TransformationModelLinear lm1(data1, param);
  TEST_REAL_SIMILAR(lm1.evaluate(2),5.285714286);

  // test evaluate using the inverted model
  lm1.invert();
  TEST_REAL_SIMILAR(lm1.evaluate(5.285714286),2);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
