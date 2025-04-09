// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: Hendrik Weisser, Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>

///////////////////////////

START_TEST(TransformationModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationModel* ptr = nullptr;
TransformationModel* nullPointer = nullptr;

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
START_SECTION((TransformationModel()))
{
  ptr = new TransformationModel();
  TEST_NOT_EQUAL(ptr, nullPointer)
  delete ptr;
}
END_SECTION

START_SECTION((TransformationModel(const DataPoints&, const Param&)))
{
  ptr = new TransformationModel(TransformationModel::DataPoints(), Param());
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~TransformationModel()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual double evaluate(double value) const))
{
  // null model (identity):
  TransformationModel tm;
  TEST_REAL_SIMILAR(tm.evaluate(-3.14159), -3.14159);
  TEST_REAL_SIMILAR(tm.evaluate(0.0), 0.0);
  TEST_REAL_SIMILAR(tm.evaluate(12345678.9), 12345678.9);
}
END_SECTION

START_SECTION((void getParameters(Param& params) const))
{
  TransformationModel tm;
  Param p = tm.getParameters();
  TEST_EQUAL(p.empty(), true)
}
END_SECTION

START_SECTION(([EXTRA] static void getDefaultParameters(Param& params)))
{
  Param param;
  param.setValue("some-value", 12.3);
  TransformationModel::getDefaultParameters(param);
  TEST_EQUAL(param.empty(), true)
}
END_SECTION

START_SECTION((bool checkValidWeight(const string& weight, const vector<string>& valid_weights) const))
{
  Param param;
  TransformationModel dw(data, param);
  string test;
  test = "ln(x)";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), true);
  test = "1/y";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidYWeights()), true);
  test = "1/x2";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), true);
  test = "x";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), true);
  test = "y";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidYWeights()), true);
  test = "none";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), false);
  test = "x2";
  TEST_EQUAL(dw.checkValidWeight(test,dw.getValidXWeights()), false);
}
END_SECTION

START_SECTION((double weightDatum(double& datum, const string& weight) const))
{
  Param param;
  TransformationModel dw(data, param);
  string test;
  test = "x";
  double inf = std::numeric_limits<double>::infinity();
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 2.0);
  test = "none";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 2.0);
  test = "ln(x)";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), -inf);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), std::log(2.0));
  test = "1/x";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), inf);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 1/std::abs(2.0));
  test = "1/x2";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), inf);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 1/std::abs(std::pow(2.0,2)));
  test = "ln(y)";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), -inf);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), std::log(2.0));
  test = "1/y";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), inf);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 1/std::abs(2.0));
  test = "1/y2";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), inf);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 1/std::abs(std::pow(2.0,2)));
}
END_SECTION

START_SECTION((double checkDatumRange(const double& datum, const double& datum_min, const double& datum_max))) //new
{
  Param param;
  TransformationModel dw(data, param);
  double dmin = 10e-6;
  double dmax = 10e9;
  TEST_REAL_SIMILAR(dw.checkDatumRange(10e-7, dmin, dmax), dmin);
  TEST_REAL_SIMILAR(dw.checkDatumRange(10e12, dmin, dmax), dmax);
  TEST_REAL_SIMILAR(dw.checkDatumRange(100, dmin, dmax), 100);
}
END_SECTION

START_SECTION((double weightDatum(double& datum, const string& weight) const))
{
  Param param;
  TransformationModel dw(data, param);
  string test;
  test = "x";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 2.0);
  TEST_REAL_SIMILAR(dw.weightDatum(10e13,test), 10e13);
  double xmin = 10e-5;
  double xmax = 10e12;
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(0.0, xmin, xmax),test), xmin);
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(2.0, xmin, xmax),test), 2.0);
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(10e13, xmin, xmax),test), xmax);
  test = "none";
  TEST_REAL_SIMILAR(dw.weightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.weightDatum(2.0,test), 2.0);
  test = "ln(x)";
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(0.0, xmin, xmax),test), std::log(xmin));
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(2.0, xmin, xmax),test), std::log(2.0));
  test = "1/x";
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(0.0, xmin, xmax),test), 1/xmin);
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(2.0, xmin, xmax),test), 1/std::abs(2.0));
  test = "1/x2";
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(0.0, xmin, xmax),test), 1/std::pow(xmin,2));
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(2.0, xmin, xmax),test), 1/std::abs(std::pow(2.0,2)));
  test = "ln(y)";
  double ymin = 10e-8;
  double ymax = 10e12;
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(0.0, ymin, ymax),test), std::log(ymin));
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(2.0, ymin, ymax),test), std::log(2.0));
  test = "1/y";
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(0.0, ymin, ymax),test), 1/ymin);
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(2.0, ymin, ymax),test), 1/std::abs(2.0));
  test = "1/y2";
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(0.0, ymin, ymax),test), 1/std::pow(ymin,2));
  TEST_REAL_SIMILAR(dw.weightDatum(dw.checkDatumRange(2.0, ymin, ymax),test), 1/std::abs(std::pow(2.0,2)));
}
END_SECTION

START_SECTION((virtual void weightData(DataPoints& data, const Param& params)))
{
  TransformationModel::DataPoints data1;
  TransformationModel::DataPoints test1;
  TransformationModel::DataPoint point;
  Param param;
  TransformationModel::getDefaultParameters(param);

  {
    double xmin = 10e-5;
    double xmax = 10e12; 
    double ymin = 10e-8;
    double ymax = 10e12;

    param.setValue("x_weight", "ln(x)");
    param.setValue("y_weight", "y");
    TransformationModel dw(data, param);
    test1.clear();
    point.first = std::log(xmin);
    point.second = 1.0;
    test1.push_back(point);
    point.first = std::log(1.0);
    point.second = 2.0;
    test1.push_back(point);
    point.first = std::log(2.0);
    point.second = 4.0;
    test1.push_back(point); 
    data1.clear();
    point.first = dw.checkDatumRange(0.0, xmin, xmax);
    point.second = dw.checkDatumRange(1.0, ymin, ymax);
    data1.push_back(point);
    point.first = dw.checkDatumRange(1.0, xmin, xmax);
    point.second = dw.checkDatumRange(2.0, ymin, ymax);
    data1.push_back(point);
    point.first = dw.checkDatumRange(2.0, xmin, xmax);
    point.second = dw.checkDatumRange(4.0, ymin, ymax);
    data1.push_back(point);
    dw.weightData(data1);
    for (size_t i = 0; i < data1.size(); ++i)
    {
      TEST_REAL_SIMILAR(data1[i].first,test1[i].first);
      TEST_REAL_SIMILAR(data1[i].second,test1[i].second);
    }
  }

  {
    param.setValue("x_weight", "x");
    param.setValue("y_weight", "ln(y)");
    TransformationModel dw(data, param);
    test1.clear();
    point.first = 0.0;
    point.second = std::log(1.0);
    test1.push_back(point);
    point.first = 1.0;
    point.second = std::log(2.0);
    test1.push_back(point);
    point.first = 2.0;
    point.second = std::log(4.0);
    test1.push_back(point);  
    data1.clear();
    point.first = 0.0;
    point.second = 1.0;
    data1.push_back(point);
    point.first = 1.0;
    point.second = 2.0;
    data1.push_back(point);
    point.first = 2.0;
    point.second = 4.0;
    data1.push_back(point);
    dw.weightData(data1);
    for (size_t i = 0; i < data1.size(); ++i)
    {
      TEST_REAL_SIMILAR(data1[i].first,test1[i].first);
      TEST_REAL_SIMILAR(data1[i].second,test1[i].second);
    }
  }
}
END_SECTION

START_SECTION((double unWeightDatum(double& datum, const string& weight) const))
{
  Param param;
  TransformationModel dw(data, param);
  string test;
  test = "x";
  TEST_REAL_SIMILAR(dw.unWeightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), 2.0);
  test = "none";
  TEST_REAL_SIMILAR(dw.unWeightDatum(0.0,test), 0.0);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), 2.0);
  test = "ln(x)";
  TEST_REAL_SIMILAR(dw.unWeightDatum(std::log(11.0e5),test), 11e5);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), std::exp(2.0));
  test = "1/x";
  TEST_REAL_SIMILAR(dw.unWeightDatum(1/std::abs(9.0e-5),test), 9.0e-5);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), 1/std::abs(2.0));
  test = "1/x2";
  TEST_REAL_SIMILAR(dw.unWeightDatum(1/std::pow(9.0e-5,2),test), 9.0e-5);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), std::sqrt(1/std::abs(2.0)));
  test = "ln(y)";
  TEST_REAL_SIMILAR(dw.unWeightDatum(std::log(11.0e8),test), 11.0e8);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), std::abs(std::exp(2.0)));
  test = "1/y";
  TEST_REAL_SIMILAR(dw.unWeightDatum(1/std::abs(9.0e-8),test), 9e-8);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), 1/std::abs(2.0));
  test = "1/y2";
  TEST_REAL_SIMILAR(dw.unWeightDatum(1/std::pow(9.0e-8,2),test), 9e-8);
  TEST_REAL_SIMILAR(dw.unWeightDatum(2.0,test), std::sqrt(1/std::abs(2.0)));
}
END_SECTION

START_SECTION((virtual void unWeightData(DataPoints& data, const Param& params)))
{

  TransformationModel::DataPoints data1;
  TransformationModel::DataPoints test1;
  TransformationModel::DataPoint point;

  {
    Param param;
    TransformationModel::getDefaultParameters(param);
    param.setValue("x_weight", "ln(x)");
    param.setValue("y_weight", "y");
    TransformationModel dw(data, param);
    test1.clear();
    point.first = std::exp(0.0);
    point.second = 1.0;
    test1.push_back(point);
    point.first = std::exp(1.0);
    point.second = 2.0;
    test1.push_back(point);
    point.first = std::exp(2.0);
    point.second = 4.0;
    test1.push_back(point);   
    data1.clear();
    point.first = 0.0;
    point.second = 1.0;
    data1.push_back(point);
    point.first = 1.0;
    point.second = 2.0;
    data1.push_back(point);
    point.first = 2.0;
    point.second = 4.0;
    data1.push_back(point);
    dw.unWeightData(data1);
    for (size_t i = 0; i < data1.size(); ++i)
    {
      TEST_REAL_SIMILAR(data1[i].first,test1[i].first);
      TEST_REAL_SIMILAR(data1[i].second,test1[i].second);
    }
  }

  {
    Param param;
    TransformationModel::getDefaultParameters(param);
    param.setValue("x_weight", "x");
    param.setValue("y_weight", "ln(y)");
    TransformationModel dw(data, param);
    test1.clear();
    point.first = 0.0;
    point.second = std::exp(1.0);
    test1.push_back(point);
    point.first = 1.0;
    point.second = std::exp(2.0);
    test1.push_back(point);
    point.first = 2.0;
    point.second = std::exp(4.0);
    test1.push_back(point);   
    data1.clear();
    point.first = 0.0;
    point.second = 1.0;
    data1.push_back(point);
    point.first = 1.0;
    point.second = 2.0;
    data1.push_back(point);
    point.first = 2.0;
    point.second = 4.0;
    data1.push_back(point);
    dw.unWeightData(data1);
    for (size_t i = 0; i < data1.size(); ++i)
    {
      TEST_REAL_SIMILAR(data1[i].first,test1[i].first);
      TEST_REAL_SIMILAR(data1[i].second,test1[i].second);
    }
  }
}
END_SECTION

START_SECTION(([EXTRA] DataPoint::DataPoint(double, double, const String&)))
{
  NOT_TESTABLE // tested below
}
END_SECTION

START_SECTION(([EXTRA] DataPoint::DataPoint(const pair<double, double>&)))
{
  NOT_TESTABLE // tested below
}
END_SECTION

START_SECTION(([EXTRA] bool DataPoint::operator<(const DataPoint& other) const))
{
  TransformationModel::DataPoint p1(1.0, 2.0, "abc");
  TransformationModel::DataPoint p2(make_pair(1.0, 2.0));
  TEST_EQUAL(p1 < p2, false);
  TEST_EQUAL(p2 < p1, true);
  p2.note = "def";
  TEST_EQUAL(p1 < p2, true);
  TEST_EQUAL(p2 < p1, false);
  p1.first = 1.5;
  TEST_EQUAL(p1 < p2, false);
  TEST_EQUAL(p2 < p1, true);
}
END_SECTION

START_SECTION(([EXTRA] bool DataPoint::operator==(const DataPoint& other) const))
{
  TransformationModel::DataPoint p1(1.0, 2.0, "abc");
  TransformationModel::DataPoint p2(make_pair(1.0, 2.0));
  TEST_EQUAL(p1 == p2, false);
  p2.note = "abc";
  TEST_TRUE(p1 == p2);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
