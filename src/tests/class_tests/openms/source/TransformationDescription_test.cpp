// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <vector>
#include <sstream>

///////////////////////////

START_TEST(TransformationDescription, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;


TransformationDescription* ptr = nullptr;
TransformationDescription* nullPointer = nullptr;
START_SECTION((TransformationDescription()))
	ptr = new TransformationDescription;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~TransformationDescription()))
	delete ptr;
END_SECTION

TransformationDescription::DataPoints data;
data.push_back(make_pair(0.0, 1.0));
data.push_back(make_pair(0.25, 1.5));
data.push_back(make_pair(0.5, 2.0));
data.push_back(make_pair(1.0, 3.0));

TransformationDescription::DataPoints data_nonlinear;
data_nonlinear.push_back(make_pair(0.0, 1.0));
data_nonlinear.push_back(make_pair(0.25, 1.0625));
data_nonlinear.push_back(make_pair(0.5, 1.25));
data_nonlinear.push_back(make_pair(1.0, 2.0));


START_SECTION((TransformationDescription(const DataPoints& data)))
{
	ptr = new TransformationDescription(data);
  TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->getDataPoints() == data, true);
	delete ptr;
}
END_SECTION

START_SECTION((const DataPoints& getDataPoints() const))
{
	TransformationDescription td;
	TEST_EQUAL(td.getDataPoints().empty(), true);
}
END_SECTION

START_SECTION((void setDataPoints(const DataPoints& data)))
{
	TransformationDescription td;
	td.fitModel("identity", Param());
	TEST_EQUAL(td.getModelType(), "identity");
	td.setDataPoints(data);
	// setting data points clears the model:
	TEST_EQUAL(td.getModelType(), "none");
	TEST_EQUAL(td.getDataPoints().size(), 4)
	TEST_EQUAL(td.getDataPoints() == data, true);

	TransformationDescription::DataPoints empty;
  empty.clear();
	td.setDataPoints(empty);
	TEST_EQUAL(td.getDataPoints().empty(), true);
}
END_SECTION

START_SECTION((void setDataPoints(const vector<pair<double, double> >& data)))
{
	TransformationDescription td;
	td.fitModel("identity", Param());
	TEST_EQUAL(td.getModelType(), "identity");
  vector<pair<double, double> > pairs;
  pairs.push_back(make_pair(0.0, 1.0));
  pairs.push_back(make_pair(0.25, 1.5));
  pairs.push_back(make_pair(0.5, 2.0));
  pairs.push_back(make_pair(1.0, 3.0));
	td.setDataPoints(pairs);
	// setting data points clears the model:
	TEST_EQUAL(td.getModelType(), "none");
	TEST_EQUAL(td.getDataPoints().size(), 4)
	TEST_EQUAL(td.getDataPoints() == data, true);

  pairs.clear();
	td.setDataPoints(pairs);
	TEST_EQUAL(td.getDataPoints().empty(), true);
}
END_SECTION

START_SECTION((double apply(double value) const))
{
	TransformationDescription td;
	TEST_EQUAL(td.apply(-0.5), -0.5);
	TEST_EQUAL(td.apply(1000), 1000);
	// tested further together with "fitModel"
}
END_SECTION

START_SECTION((const String& getModelType() const))
{
	TransformationDescription td;
	TEST_EQUAL(td.getModelType(), "none");
}
END_SECTION

START_SECTION((static void getModelTypes(StringList& result)))
{
	StringList result;
	TransformationDescription::getModelTypes(result);
	TEST_EQUAL(result.size(), 4);
	TEST_EQUAL(result.at(0), "linear");
	TEST_EQUAL(result.at(1), "b_spline");
	TEST_EQUAL(result.at(2), "interpolated");
	TEST_EQUAL(result.at(3), "lowess");
}
END_SECTION

START_SECTION((void fitModel(const String& model_type, const Param& params=Param())))
{
	TransformationDescription td(data);
	Param params;
	td.fitModel("linear", params);
	TEST_EQUAL(td.getModelType(), "linear");
	TEST_REAL_SIMILAR(td.apply(0.0), 1.0);
	TEST_REAL_SIMILAR(td.apply(0.5), 2.0);
	TEST_REAL_SIMILAR(td.apply(1.0), 3.0);

  // non-linear model (b spline)
	td.fitModel("b_spline", params);
	TEST_EQUAL(td.getModelType(), "b_spline");
	TEST_REAL_SIMILAR(td.apply(0.0), 1.064201730);
	TEST_REAL_SIMILAR(td.apply(0.5), 1.957836652);
	TEST_REAL_SIMILAR(td.apply(1.0), 2.927541901);

  // non-linear model (lowess)
	td.fitModel("lowess", params);
	TEST_EQUAL(td.getModelType(), "lowess");
	TEST_REAL_SIMILAR(td.apply(0.0), 1.0);
	TEST_REAL_SIMILAR(td.apply(0.5), 2.0);
	TEST_REAL_SIMILAR(td.apply(1.0), 3.0);

  // special model type for reference files:
	td.fitModel("identity", Param());
	TEST_EQUAL(td.getModelType(), "identity");
	TEST_REAL_SIMILAR(td.apply(0.0), 0.0);
	TEST_REAL_SIMILAR(td.apply(0.5), 0.5);
	TEST_REAL_SIMILAR(td.apply(1.0), 1.0);

	// can't fit a different model to an "identity" transformation:
	td.fitModel("linear", params);
	TEST_EQUAL(td.getModelType(), "identity");

  {
    // non-linear model (b spline)
    TransformationDescription td_nl(data_nonlinear);
    td_nl.fitModel("b_spline", params);
    TEST_EQUAL(td_nl.getModelType(), "b_spline");
    TEST_REAL_SIMILAR(td_nl.apply(0.0), 1.01084556836969);
    TEST_REAL_SIMILAR(td_nl.apply(0.5), 1.26289804387079);
    TEST_REAL_SIMILAR(td_nl.apply(0.75), 1.53463130131214);
    TEST_REAL_SIMILAR(td_nl.apply(1.0), 1.94984504419826);

    // non-linear model (lowess)
    td_nl.fitModel("lowess", params);
    TEST_EQUAL(td_nl.getModelType(), "lowess");
    TEST_REAL_SIMILAR(td_nl.apply(0.0), 1.0);
    TEST_REAL_SIMILAR(td_nl.apply(0.5), 1.25);
    TEST_REAL_SIMILAR(td_nl.apply(0.75), 1.58423913043478);
    TEST_REAL_SIMILAR(td_nl.apply(1.0), 2.0);
  }
}
END_SECTION

START_SECTION(([EXTRA]void fitModel(const String& model_type, const Param& params=Param())))
{
  // Check whether we can change the parameters and get different behavior
	TransformationDescription td(data);
	Param params;

  // for lowess
  params.setValue("interpolation_type", "linear");

  // for b spline
  params.setValue("extrapolate", "b_spline");

  // non-linear model (b spline)
  TransformationDescription td_nl(data_nonlinear);
  td_nl.fitModel("b_spline", params);
  TEST_EQUAL(td_nl.getModelType(), "b_spline");
  TEST_REAL_SIMILAR(td_nl.apply(0.0), 1.01084556836969);
  TEST_REAL_SIMILAR(td_nl.apply(0.5), 1.26289804387079);
  TEST_REAL_SIMILAR(td_nl.apply(0.75), 1.53463130131214);
  TEST_REAL_SIMILAR(td_nl.apply(1.0), 1.94984504419826);
  TEST_REAL_SIMILAR(td_nl.apply(2.0), 1.328125); // b-spline extrapolation

  // non-linear model (lowess)
  td_nl.fitModel("lowess", params);
  TEST_EQUAL(td_nl.getModelType(), "lowess");
  TEST_REAL_SIMILAR(td_nl.apply(0.0), 1.0);
  TEST_REAL_SIMILAR(td_nl.apply(0.5), 1.25);
  TEST_REAL_SIMILAR(td_nl.apply(0.75), 1.625); // linear interpolation between points
  TEST_REAL_SIMILAR(td_nl.apply(1.0), 2.0);
  TEST_REAL_SIMILAR(td_nl.apply(2.0), 3.5);
}
END_SECTION

START_SECTION((void getModelParameters(Param& params) const))
{
	TransformationDescription td;
	Param params = td.getModelParameters();
	TEST_EQUAL(params, Param());
	params.setValue("slope", 2.5);
	params.setValue("intercept", -100.0);
  params.setValue("x_weight", "x");
  params.setValue("y_weight", "y");
  params.setValue("x_datum_min", 1e-15);
  params.setValue("y_datum_min", 1e-15);
  params.setValue("x_datum_max", 1e15);
  params.setValue("y_datum_max", 1e15);
	const Param const_params = params;
	td.fitModel("linear", const_params);
	params = td.getModelParameters();
	TEST_EQUAL(params, const_params);
}
END_SECTION

START_SECTION((TransformationDescription(const TransformationDescription& rhs)))
{
	TransformationDescription td(data);
	td.fitModel("linear", Param());
	TransformationDescription td2 = td;
	TEST_EQUAL(td.getModelType(), td2.getModelType());
	TEST_EQUAL(td.getDataPoints() == td2.getDataPoints(), true);
	Param params = td.getModelParameters();
	Param params2 = td2.getModelParameters();
	TEST_EQUAL(params, params2);
}
END_SECTION

START_SECTION((TransformationDescription& operator=(const TransformationDescription& rhs)))
{
	TransformationDescription td(data);
	td.fitModel("linear", Param());
	TransformationDescription td2;
	td2 = td;
	TEST_EQUAL(td.getModelType(), td2.getModelType());
	TEST_EQUAL(td.getDataPoints() == td2.getDataPoints(), true);
  Param params = td.getModelParameters();
  Param params2 = td2.getModelParameters();
	TEST_EQUAL(params, params2);
}
END_SECTION

START_SECTION((void invert()))
{
	// test null transformation:
	TransformationDescription td;
	td.fitModel("none", Param());
	td.invert();
	TEST_EQUAL(td.getModelType(), "none");

  // test inversion of data points:
	TransformationDescription td1;
  td1.setDataPoints(data);
  td1.invert();
  td1.invert();
  TEST_EQUAL(td1.getDataPoints() == data, true);

	// test linear transformation:
  TransformationDescription td2;
  td2.setDataPoints(data);
  Param params2;
	params2.setValue("slope", 2.0);
	params2.setValue("intercept", 47.12);
	td2.fitModel("linear", params2);
	TEST_REAL_SIMILAR(params2.getValue("slope"), 2.0);
  TEST_REAL_SIMILAR(params2.getValue("intercept"), 47.12);

  TransformationDescription td3;
  td3.setDataPoints(data);
  Param params3;
  td3.fitModel("linear", params3);
	td3.invert();
	TEST_EQUAL(td3.getModelType(), "linear");
	TEST_REAL_SIMILAR(td3.apply(1.0), 0.0);//control input values
	TEST_REAL_SIMILAR(td3.apply(2.0), 0.5);//control input values
	TEST_REAL_SIMILAR(td3.apply(3.0), 1.0);//control input values
	TEST_REAL_SIMILAR(td3.apply(4.0), 1.5);//control interpolation values
	TEST_REAL_SIMILAR(td3.apply(5.0), 2.0);//control interpolation values

	// test interpolated-linear transformation:
	TransformationDescription td4;
	td4.setDataPoints(data);
  Param params4;
	params4.setValue("interpolation_type", "cspline");
	td4.fitModel("interpolated", params4);
	td4.invert();
	TEST_EQUAL(td4.getModelType(), "interpolated");
	// pairs have changed...
	TEST_EQUAL(td4.getDataPoints() != data, true);
	td4.invert();
	// ... now they're back to the original:
	TEST_EQUAL(td4.getDataPoints() == data, true);
}
END_SECTION

START_SECTION((void getDeviations(std::vector<double>& diffs, bool do_apply = false, bool do_sort = true) const))
{
  vector<double> diffs;
	TransformationDescription td(data_nonlinear);
  td.fitModel("linear");
  td.getDeviations(diffs);
  TEST_EQUAL(diffs.size(), 4);
  TEST_REAL_SIMILAR(diffs[0], 0.75);
  TEST_REAL_SIMILAR(diffs[1], 0.8125);
  TEST_REAL_SIMILAR(diffs[2], 1.0);
  TEST_REAL_SIMILAR(diffs[3], 1.0);

  td.getDeviations(diffs, true, false);
  TEST_EQUAL(diffs.size(), 4);
  TEST_REAL_SIMILAR(diffs[0], 0.125);
  TEST_REAL_SIMILAR(diffs[1], 0.0714286);
  TEST_REAL_SIMILAR(diffs[2], 0.142857);
  TEST_REAL_SIMILAR(diffs[3], 0.0892857);
}
END_SECTION

START_SECTION((void printSummary(std::ostream& os) const))
{
  stringstream ss;
	TransformationDescription td(data_nonlinear);
  td.printSummary(ss);
  string expected =
    "Number of data points (x/y pairs): 4\n"
    "Data range (x): 0 to 1\n"
    "Data range (y): 1 to 2\n"
    "Summary of x/y deviations:\n"
    "- 100% of data points within (+/-)1\n"
    "-  99% of data points within (+/-)1\n"
    "-  95% of data points within (+/-)1\n"
    "-  90% of data points within (+/-)1\n"
    "-  75% of data points within (+/-)1\n"
    "-  50% of data points within (+/-)0.8125\n"
    "-  25% of data points within (+/-)0.75\n\n";
  TEST_STRING_EQUAL(ss.str(), expected);

  ss.str("");
	td.fitModel("linear");
  td.printSummary(ss);
  expected =
    "Number of data points (x/y pairs): 4\n"
    "Data range (x): 0 to 1\n"
    "Data range (y): 1 to 2\n"
    "Summary of x/y deviations before transformation:\n"
    "- 100% of data points within (+/-)1\n"
    "-  99% of data points within (+/-)1\n"
    "-  95% of data points within (+/-)1\n"
    "-  90% of data points within (+/-)1\n"
    "-  75% of data points within (+/-)1\n"
    "-  50% of data points within (+/-)0.8125\n"
    "-  25% of data points within (+/-)0.75\n"
    "Summary of x/y deviations after applying 'linear' transformation:\n"
    "- 100% of data points within (+/-)0.142857\n"
    "-  99% of data points within (+/-)0.125\n"
    "-  95% of data points within (+/-)0.125\n"
    "-  90% of data points within (+/-)0.125\n"
    "-  75% of data points within (+/-)0.125\n"
    "-  50% of data points within (+/-)0.0892857\n"
    "-  25% of data points within (+/-)0.0714286\n\n";
  TEST_STRING_EQUAL(ss.str(), expected);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
