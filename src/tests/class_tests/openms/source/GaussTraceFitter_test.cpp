// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FEATUREFINDER/GaussTraceFitter.h>
///////////////////////////

#include <OpenMS/KERNEL/Peak1D.h>

using namespace OpenMS;
using namespace std;

START_TEST(GaussTraceFitter, "$Id$")

/////////////////////////////////////////////////////////////
typedef GaussTraceFitter GTF;
/////////////////////////////////////////////////////////////
FeatureFinderAlgorithmPickedHelperStructs::MassTraces mts;

FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt1;
mt1.theoretical_int = 0.8;

FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt2;
mt2.theoretical_int = 0.2;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// set up mass traces to fit

Peak1D p1_1;
p1_1.setIntensity(1.08268226589f);
p1_1.setMZ(1000);
mt1.peaks.push_back(make_pair(677.1, &p1_1));
Peak1D p2_1;
p2_1.setIntensity(0.270670566473f);
p2_1.setMZ(1001);
mt2.peaks.push_back(make_pair(677.1, &p2_1));
Peak1D p1_2;
p1_2.setIntensity(1.58318959267f);
p1_2.setMZ(1000);
mt1.peaks.push_back(make_pair(677.4, &p1_2));
Peak1D p2_2;
p2_2.setIntensity(0.395797398167f);
p2_2.setMZ(1001);
mt2.peaks.push_back(make_pair(677.4, &p2_2));
Peak1D p1_3;
p1_3.setIntensity(2.22429840363f);
p1_3.setMZ(1000);
mt1.peaks.push_back(make_pair(677.7, &p1_3));
Peak1D p2_3;
p2_3.setIntensity(0.556074600906f);
p2_3.setMZ(1001);
mt2.peaks.push_back(make_pair(677.7, &p2_3));
Peak1D p1_4;
p1_4.setIntensity(3.00248879081f);
p1_4.setMZ(1000);
mt1.peaks.push_back(make_pair(678, &p1_4));
Peak1D p2_4;
p2_4.setIntensity(0.750622197703f);
p2_4.setMZ(1001);
mt2.peaks.push_back(make_pair(678, &p2_4));
Peak1D p1_5;
p1_5.setIntensity(3.89401804768f);
p1_5.setMZ(1000);
mt1.peaks.push_back(make_pair(678.3, &p1_5));
Peak1D p2_5;
p2_5.setIntensity(0.97350451192f);
p2_5.setMZ(1001);
mt2.peaks.push_back(make_pair(678.3, &p2_5));
Peak1D p1_6;
p1_6.setIntensity(4.8522452777f);
p1_6.setMZ(1000);
mt1.peaks.push_back(make_pair(678.6, &p1_6));
Peak1D p2_6;
p2_6.setIntensity(1.21306131943f);
p2_6.setMZ(1001);
mt2.peaks.push_back(make_pair(678.6, &p2_6));
Peak1D p1_7;
p1_7.setIntensity(5.80919229659f);
p1_7.setMZ(1000);
mt1.peaks.push_back(make_pair(678.9, &p1_7));
Peak1D p2_7;
p2_7.setIntensity(1.45229807415f);
p2_7.setMZ(1001);
mt2.peaks.push_back(make_pair(678.9, &p2_7));
Peak1D p1_8;
p1_8.setIntensity(6.68216169129f);
p1_8.setMZ(1000);
mt1.peaks.push_back(make_pair(679.2, &p1_8));
Peak1D p2_8;
p2_8.setIntensity(1.67054042282f);
p2_8.setMZ(1001);
mt2.peaks.push_back(make_pair(679.2, &p2_8));
Peak1D p1_9;
p1_9.setIntensity(7.38493077109f);
p1_9.setMZ(1000);
mt1.peaks.push_back(make_pair(679.5, &p1_9));
Peak1D p2_9;
p2_9.setIntensity(1.84623269277f);
p2_9.setMZ(1001);
mt2.peaks.push_back(make_pair(679.5, &p2_9));
Peak1D p1_10;
p1_10.setIntensity(7.84158938645f);
p1_10.setMZ(1000);
mt1.peaks.push_back(make_pair(679.8, &p1_10));
Peak1D p2_10;
p2_10.setIntensity(1.96039734661f);
p2_10.setMZ(1001);
mt2.peaks.push_back(make_pair(679.8, &p2_10));
Peak1D p1_11;
p1_11.setIntensity(8.0f);
p1_11.setMZ(1000);
mt1.peaks.push_back(make_pair(680.1, &p1_11));
Peak1D p2_11;
p2_11.setIntensity(2.0f);
p2_11.setMZ(1001);
mt2.peaks.push_back(make_pair(680.1, &p2_11));
Peak1D p1_12;
p1_12.setIntensity(7.84158938645f);
p1_12.setMZ(1000);
mt1.peaks.push_back(make_pair(680.4, &p1_12));
Peak1D p2_12;
p2_12.setIntensity(1.96039734661f);
p2_12.setMZ(1001);
mt2.peaks.push_back(make_pair(680.4, &p2_12));
Peak1D p1_13;
p1_13.setIntensity(7.38493077109f);
p1_13.setMZ(1000);
mt1.peaks.push_back(make_pair(680.7, &p1_13));
Peak1D p2_13;
p2_13.setIntensity(1.84623269277f);
p2_13.setMZ(1001);
mt2.peaks.push_back(make_pair(680.7, &p2_13));
Peak1D p1_14;
p1_14.setIntensity(6.68216169129f);
p1_14.setMZ(1000);
mt1.peaks.push_back(make_pair(681, &p1_14));
Peak1D p2_14;
p2_14.setIntensity(1.67054042282f);
p2_14.setMZ(1001);
mt2.peaks.push_back(make_pair(681, &p2_14));
Peak1D p1_15;
p1_15.setIntensity(5.80919229659f);
p1_15.setMZ(1000);
mt1.peaks.push_back(make_pair(681.3, &p1_15));
Peak1D p2_15;
p2_15.setIntensity(1.45229807415f);
p2_15.setMZ(1001);
mt2.peaks.push_back(make_pair(681.3, &p2_15));
Peak1D p1_16;
p1_16.setIntensity(4.8522452777f);
p1_16.setMZ(1000);
mt1.peaks.push_back(make_pair(681.6, &p1_16));
Peak1D p2_16;
p2_16.setIntensity(1.21306131943f);
p2_16.setMZ(1001);
mt2.peaks.push_back(make_pair(681.6, &p2_16));
Peak1D p1_17;
p1_17.setIntensity(3.89401804768f);
p1_17.setMZ(1000);
mt1.peaks.push_back(make_pair(681.9, &p1_17));
Peak1D p2_17;
p2_17.setIntensity(0.97350451192f);
p2_17.setMZ(1001);
mt2.peaks.push_back(make_pair(681.9, &p2_17));
Peak1D p1_18;
p1_18.setIntensity(3.00248879081f);
p1_18.setMZ(1000);
mt1.peaks.push_back(make_pair(682.2, &p1_18));
Peak1D p2_18;
p2_18.setIntensity(0.750622197703f);
p2_18.setMZ(1001);
mt2.peaks.push_back(make_pair(682.2, &p2_18));
Peak1D p1_19;
p1_19.setIntensity(2.22429840363f);
p1_19.setMZ(1000);
mt1.peaks.push_back(make_pair(682.5, &p1_19));
Peak1D p2_19;
p2_19.setIntensity(0.556074600906f);
p2_19.setMZ(1001);
mt2.peaks.push_back(make_pair(682.5, &p2_19));
Peak1D p1_20;
p1_20.setIntensity(1.58318959267f);
p1_20.setMZ(1000);
mt1.peaks.push_back(make_pair(682.8, &p1_20));
Peak1D p2_20;
p2_20.setIntensity(0.395797398167f);
p2_20.setMZ(1001);
mt2.peaks.push_back(make_pair(682.8, &p2_20));
Peak1D p1_21;
p1_21.setIntensity(1.08268226589f);
p1_21.setMZ(1000);
mt1.peaks.push_back(make_pair(683.1, &p1_21));
Peak1D p2_21;
p2_21.setIntensity(0.270670566473f);
p2_21.setMZ(1001);
mt2.peaks.push_back(make_pair(683.1, &p2_21));

mt1.updateMaximum();
mts.push_back(mt1);

mt2.updateMaximum();
mts.push_back(mt2);

// fix base line to 0 since we have no baseline here
mts.baseline = 0.0;

mts.max_trace = 0;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// setup fitter

Param p;
p.setValue("max_iteration", 500);

GTF gaussian_trace_fitter;
gaussian_trace_fitter.setParameters(p);
gaussian_trace_fitter.fit(mts);

double expected_sigma = 1.5;
double expected_H = 10.0;
double expected_x0 = 680.1;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
TOLERANCE_RELATIVE(1.001)


GTF* ptr = nullptr;
GTF* nullPointer = nullptr;
START_SECTION(GaussTraceFitter())
{
	ptr = new GTF();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~GaussTraceFitter())
{
	delete ptr;
}
END_SECTION

START_SECTION((GaussTraceFitter(const GaussTraceFitter& other)))
{
  GTF gtf2(gaussian_trace_fitter);

  TEST_REAL_SIMILAR(gaussian_trace_fitter.getCenter(), gtf2.getCenter())
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getHeight(), gtf2.getHeight())
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getSigma(), gtf2.getSigma())
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getLowerRTBound(), gtf2.getLowerRTBound())
}
END_SECTION

START_SECTION((GaussTraceFitter& operator=(const GaussTraceFitter& source)))
{
  GTF gtf3 = gaussian_trace_fitter;

  TEST_REAL_SIMILAR(gaussian_trace_fitter.getCenter(), gtf3.getCenter())
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getHeight(), gtf3.getHeight())
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getSigma(), gtf3.getSigma())
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getLowerRTBound(), gtf3.getLowerRTBound())
}
END_SECTION

START_SECTION((void fit(FeatureFinderAlgorithmPickedHelperStructs::MassTraces& traces)))
{
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getCenter(), expected_x0)
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getHeight(), expected_H)
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getSigma(), expected_sigma)
  GaussTraceFitter weighted_fitter;
  Param params = weighted_fitter.getDefaults();
  params.setValue("weighted", "true");
  weighted_fitter.setParameters(params);
  weighted_fitter.fit(mts);
  TEST_REAL_SIMILAR(weighted_fitter.getCenter(), expected_x0)
  TEST_REAL_SIMILAR(weighted_fitter.getHeight(), expected_H - 0.0035)
  mts[0].theoretical_int = 0.4;
  mts[1].theoretical_int = 0.6;
  weighted_fitter.fit(mts);
  TEST_REAL_SIMILAR(weighted_fitter.getCenter(), expected_x0)
    TEST_REAL_SIMILAR(weighted_fitter.getHeight(), 6.0847)
}
END_SECTION

START_SECTION((double getLowerRTBound() const))
{
  // given sigma this should be
  // x0_ - 2.5 * sigma_;
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getLowerRTBound(), expected_x0 - 2.5 * expected_sigma)
}
END_SECTION

START_SECTION((double getUpperRTBound() const))
{
  // given sigma this should be
  // x0_ + 2.5 * sigma_;
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getUpperRTBound(), expected_x0 + 2.5 * expected_sigma)
}
END_SECTION

START_SECTION((double getHeight() const))
{
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getHeight(), 10.0)
}
END_SECTION

START_SECTION((double getCenter() const))
{
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getCenter(), 680.1)
}
END_SECTION

START_SECTION((double getSigma() const))
{
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getSigma(), 1.5)
}
END_SECTION

START_SECTION((bool checkMaximalRTSpan(const double max_rt_span)))
{
  // Maximum RT span in relation to extended area that the model is allowed to have
  // 5.0 * sigma_ > max_rt_span * region_rt_span_

  double region_rt_span = mt1.peaks[mt1.peaks.size() - 1].first - mt1.peaks[0].first;
  double max_rt_span = 5.0 * gaussian_trace_fitter.getSigma() / region_rt_span + 0.00000000000001; // we add some small number to overcome precision problems on 32-bit machines

  TEST_EQUAL(gaussian_trace_fitter.checkMaximalRTSpan(max_rt_span), false);
  max_rt_span -= 0.1; // accept only smaller regions
  TEST_EQUAL(gaussian_trace_fitter.checkMaximalRTSpan(max_rt_span), true);
}
END_SECTION

START_SECTION((bool checkMinimalRTSpan(const std::pair<double, double> &rt_bounds, const double min_rt_span)))
{
  // is
  // (rt_bounds.second-rt_bounds.first) < min_rt_span * 5.0 * sigma_;
  // Minimum RT span in relation to extended area that has to remain after model fitting.

  pair<double, double> rt_bounds = make_pair(0.0,4.0);
  double min_rt_span = 0.5;

  TEST_EQUAL(gaussian_trace_fitter.checkMinimalRTSpan(rt_bounds, min_rt_span), false)
  min_rt_span += 0.5;
  TEST_EQUAL(gaussian_trace_fitter.checkMinimalRTSpan(rt_bounds, min_rt_span), true)
}
END_SECTION

START_SECTION((double getValue(double rt) const))
{
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getValue(expected_x0), expected_H)
}
END_SECTION

START_SECTION((double computeTheoretical(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, Size k)))
{
  FeatureFinderAlgorithmPickedHelperStructs::MassTrace mt;
  mt.theoretical_int = 0.8;

  Peak1D peak;
  peak.setIntensity(8.0);

  mt.peaks.push_back(make_pair(expected_x0, &peak));

  // theoretical should be expected_H * theoretical_int at position expected_x0
  TEST_REAL_SIMILAR(gaussian_trace_fitter.computeTheoretical(mt, 0), mt.theoretical_int * expected_H)
}
END_SECTION

START_SECTION((virtual double getArea()))
{
  // is 2.5... * height_ * sigma_;
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getArea(), 2.506628 * expected_sigma * expected_H)
}
END_SECTION

START_SECTION((virtual String getGnuplotFormula(const FeatureFinderAlgorithmPickedHelperStructs::MassTrace& trace, const char function_name, const double baseline, const double rt_shift)))
{
  String formula = gaussian_trace_fitter.getGnuplotFormula(mts[0], 'f', 0.0, 0.0);
  // should look like -- f(x)= 0 + 7.99996 * exp(-0.5*(x-680.1)**2/(1.50001)**2) --
  TEST_EQUAL(formula.hasPrefix("f(x)= 0 + "), true)
  TEST_EQUAL(formula.hasSubstring("exp(-0.5*(x-"), true)
  TEST_EQUAL(formula.hasSubstring(")**2/("), true)
  TEST_EQUAL(formula.hasSuffix(")**2)"), true)
}
END_SECTION

START_SECTION((double getFWHM() const))
{
  TEST_REAL_SIMILAR(gaussian_trace_fitter.getFWHM(), 2.35482 * expected_sigma)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



