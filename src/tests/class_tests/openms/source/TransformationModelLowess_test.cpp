// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>

///////////////////////////

START_TEST(TransformationModelLowess, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TransformationModelLowess* ptr = nullptr;

TransformationModel::DataPoints data, empty;
data.push_back(make_pair(1.2, 5.2));
data.push_back(make_pair(3.2, 7.3));
data.push_back(make_pair(2.2, 6.25));
data.push_back(make_pair(2.2, 3.1));
data.push_back(make_pair(2.2, 7.25));
data.push_back(make_pair(3.0, 8.5));
data.push_back(make_pair(3.1, 4.7));
data.push_back(make_pair(1.7, 6.0));
data.push_back(make_pair(2.9, 4.7));
data.push_back(make_pair(4.2, 5.0));
data.push_back(make_pair(3.7, -2.4));

Param params_default;
TransformationModelLowess::getDefaultParameters(params_default);

START_SECTION((TransformationModelLowess(const DataPoints&, const Param&)))
{
  TEST_EXCEPTION(Exception::IllegalArgument, TransformationModelLowess tm(empty, params_default)); // need data
  ptr = new TransformationModelLowess(data, params_default);
  TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((~TransformationModelLowess()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual double evaluate(double value) const))
{
  // test data: sine function with added noise
  double x[] = {-0.547062107104045, -2.14564213748743, -3.07082880304281, 0.470273389368586, 1.79367651606654, 0.595846950617167, 1.58738829599701, -3.11534942614546, -2.55761408378404, -0.996199010293142, -0.553164304142189, 3.11858532047631, 0.74970539948485, 0.276411185223925, 1.85962696821902, 0.960234253336655, -1.62536120645258, -2.72457034250236, 1.67812366716942, -0.838775352531627, -0.654629712755158, 1.8220799029759, -1.8653140724926, -0.235789436296459, -0.29890807257244, 0.405216494893513, 0.233453956340058, -2.82471832316488, -3.08393846252989, -1.41524590344969, -0.199886448130033};
  double y[] = {-0.584809756448807, -0.866407723341462, -0.0471640435125096, 0.435337754412529, 0.861949333280581, 0.616243288851563, 1.1228424073836, -0.0483419751019981, -0.532873307735754, -0.917205998701872, -0.301045308942404, 0.0120964875551685, 0.758584328691163, 0.405241179450931, 1.00118722437611, 0.765459021914008, -1.03191739643009, -0.477999500942485, 0.872168291767237, -0.770691257861706, -0.496027498267174, 0.743777383059081, -0.982264617804229, -0.398462173815226, -0.40498973770553, 0.348305878579121, 0.0755855659375029, -0.457381746018402, 0.245483195014945, -1.07618910469392, -0.0880708165561682};
  // results validated by visual inspection:
  // double pred_spl[] = {0.846137, 0.689856, 0.5094,   0.31183, 0.10421,    -0.106399, -0.312921, -0.508271, -0.685362, -0.837111, -0.95643, -1.03623, -1.06944, -1.05016, -0.981868, -0.872412, -0.729666, -0.561505, -0.375803, -0.180434, 0.016728, 0.20827, 0.38867, 0.55289, 0.695895, 0.812645, 0.898104, 0.947234, 0.955013, 0.919484, 0.845545, 0.739022, 0.60574, 0.451526, 0.282206, 0.103606, -0.0784482, -0.258084, -0.429374, -0.586387, -0.723191};
  double pred_low[] = {0.81549, 0.62377, 0.432051, 0.240331, 0.0486111, -0.142406, -0.326079, -0.505898, -0.656341, -0.774182, -0.879421, -0.948104, -0.973173, -0.977378, -0.919574, -0.829349, -0.716313, -0.546735, -0.398708, -0.236083, -0.0312524, 0.171936, 0.376142, 0.553262, 0.711474, 0.825752, 0.868956, 0.873851, 0.881359, 0.905248, 0.888206, 0.814595, 0.694699, 0.538594, 0.356359, 0.158072, -0.0192127, -0.157386, -0.295559, -0.433733, -0.571906};

  /* R code:
   *
   * Note that compared to the spline, the lowess contains a linear model and
   * is thus less suited for highly non-linear data such as a sine. However, it
   * performs quite well on the data.
   *

    x = c(-0.547062107104045, -2.14564213748743, -3.07082880304281, 0.470273389368586, 1.79367651606654, 0.595846950617167, 1.58738829599701, -3.11534942614546, -2.55761408378404, -0.996199010293142, -0.553164304142189, 3.11858532047631, 0.74970539948485, 0.276411185223925, 1.85962696821902, 0.960234253336655, -1.62536120645258, -2.72457034250236, 1.67812366716942, -0.838775352531627, -0.654629712755158, 1.8220799029759, -1.8653140724926, -0.235789436296459, -0.29890807257244, 0.405216494893513, 0.233453956340058, -2.82471832316488, -3.08393846252989, -1.41524590344969, -0.199886448130033)
    y = c(-0.584809756448807, -0.866407723341462, -0.0471640435125096, 0.435337754412529, 0.861949333280581, 0.616243288851563, 1.1228424073836, -0.0483419751019981, -0.532873307735754, -0.917205998701872, -0.301045308942404, 0.0120964875551685, 0.758584328691163, 0.405241179450931, 1.00118722437611, 0.765459021914008, -1.03191739643009, -0.477999500942485, 0.872168291767237, -0.770691257861706, -0.496027498267174, 0.743777383059081, -0.982264617804229, -0.398462173815226, -0.40498973770553, 0.348305878579121, 0.0755855659375029, -0.457381746018402, 0.245483195014945, -1.07618910469392, -0.0880708165561682)

    pred_spl = c(0.846137, 0.689856, 0.5094,   0.31183, 0.10421,    -0.106399, -0.312921, -0.508271, -0.685362, -0.837111, -0.95643, -1.03623, -1.06944, -1.05016, -0.981868, -0.872412, -0.729666, -0.561505, -0.375803, -0.180434, 0.016728, 0.20827, 0.38867, 0.55289, 0.695895, 0.812645, 0.898104, 0.947234, 0.955013, 0.919484, 0.845545, 0.739022, 0.60574, 0.451526, 0.282206, 0.103606, -0.0784482, -0.258084, -0.429374, -0.586387, -0.723191)
    pred_low = c(0.81549, 0.62377, 0.432051, 0.240331, 0.0486111, -0.142406, -0.326079, -0.505898, -0.656341, -0.774182, -0.879421, -0.948104, -0.973173, -0.977378, -0.919574, -0.829349, -0.716313, -0.546735, -0.398708, -0.236083, -0.0312524, 0.171936, 0.376142, 0.553262, 0.711474, 0.825752, 0.868956, 0.873851, 0.881359, 0.905248, 0.888206, 0.814595, 0.694699, 0.538594, 0.356359, 0.158072, -0.0192127, -0.157386, -0.295559, -0.433733, -0.571906)


    plot(x, y, xlim=c(-4,4))
    lines( seq(-4, 4.1, 0.2), pred_spl, col="blue")
    lines( seq(-4, 4.1, 0.2), pred_low, col="red")

  */

  data.resize(31);
  for (Size i = 0; i < 31; ++i)
  {
    data[i] = make_pair(x[i], y[i]);
  }

  Param params;
  params.setValue("span", 0.3); // use a very low span value for non-linear data (only use 30 % of all points at a time ...)
  params.setValue("num_iterations", 3);
  params.setValue("delta", -1.0);
  params.setValue("interpolation_type", "cspline");
  params.setValue("extrapolation_type", "four-point-linear");
  TransformationModelLowess tm(data, params);
  
  vector<double> results;
  Size index = 0;
  for (double v = -4; v < 4.1; v += 0.2, index++)
  {
    TEST_REAL_SIMILAR(tm.evaluate(v), pred_low[index]);
  }

  // test extrapolation:
  params.setValue("extrapolation_type", "four-point-linear");
  TransformationModelLowess tm_lin(data, params);
  TEST_REAL_SIMILAR(tm_lin.evaluate(-4.0), 0.815490292172986);
  TEST_REAL_SIMILAR(tm_lin.evaluate(4.0), -0.571905836956494);

  params.setValue("extrapolation_type", "two-point-linear");
  TransformationModelLowess tm_const(data, params);
  TEST_REAL_SIMILAR(tm_const.evaluate(-4.0), -0.04240732863);
  TEST_REAL_SIMILAR(tm_const.evaluate(4.0), 0.046870277);

  params.setValue("extrapolation_type", "global-linear");
  TransformationModelLowess tm_global(data, params);
  TEST_REAL_SIMILAR(tm_global.evaluate(-4.0), -0.9501004);
  TEST_REAL_SIMILAR(tm_global.evaluate(4.0), 1.08486397);
}
END_SECTION

START_SECTION((void getParameters(Param& params) const))
{
  Param p_in;

  p_in.setValue("span", 0.3);
  p_in.setValue("num_iterations", 8);
  p_in.setValue("delta", 1.0);
  p_in.setValue("extrapolate", "b_spline");
  p_in.setValue("interpolation_type", "cspline");
  p_in.setValue("extrapolation_type", "four-point-linear");

  TransformationModelLowess tm(data, p_in);
  TEST_EQUAL(tm.getParameters().getValue("num_iterations"),
             p_in.getValue("num_iterations"));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
