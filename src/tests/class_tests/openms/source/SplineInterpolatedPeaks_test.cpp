// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>

#include <OpenMS/PROCESSING/MISC/SplineInterpolatedPeaks.h>

#include <algorithm>
#include <random>

using namespace OpenMS;

double Gauss1(double x)
{
    return exp(- pow(x-416.8, 2)/(2*0.15*0.15));
}

double Gauss2(double x)
{
    return exp(- pow(x-418.7, 2)/(2*0.15*0.15));
}

START_TEST(SplineInterpolatedPeaks, "$Id$")

std::vector<double> pos;
std::vector<double> intensity;
for (int i=0; i < 11; ++i)
{
    pos.push_back(416.3 + 0.1*i);
    intensity.push_back(Gauss1(416.3+0.1*i));
}
for (int i=0; i < 11; ++i)
{
    pos.push_back(418.2 + 0.1*i);
    intensity.push_back(Gauss2(418.2+0.1*i));
}

MSSpectrum spectrum;
Peak1D peak;
spectrum.setRT(1789.0714);
for (size_t i=0; i < pos.size(); ++i)
{
    peak.setMZ(pos[i]);
    peak.setIntensity(intensity[i]);
    spectrum.push_back(peak);
}

MSChromatogram chromatogram;
ChromatogramPeak peak_c;
for (size_t i=0; i < pos.size(); ++i)
{
    peak_c.setRT(pos[i]);
    peak_c.setIntensity(intensity[i]);
    chromatogram.push_back(peak_c);
}

SplineInterpolatedPeaks* nullPointer = nullptr;
SplineInterpolatedPeaks* ptr;

START_SECTION(SplineInterpolatedPeaks(const std::vector<double>& pos, const std::vector<double>& intensity))
    SplineInterpolatedPeaks spline(pos, intensity);
    TEST_REAL_SIMILAR(spline.getPosMin(), 416.3);
    ptr = new SplineInterpolatedPeaks(pos, intensity);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(SplineInterpolatedPeaks(const MSSpectrum& raw_spectrum))
	SplineInterpolatedPeaks spline(spectrum);
    TEST_REAL_SIMILAR(spline.getPosMin(), 416.3)
    ptr = new SplineInterpolatedPeaks(spectrum);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

START_SECTION(SplineInterpolatedPeaks(const MSChromatogram& raw_chromatogram))
	SplineInterpolatedPeaks spline(chromatogram);
    TEST_REAL_SIMILAR(spline.getPosMin(), 416.3)
    ptr = new SplineInterpolatedPeaks(chromatogram);
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

SplineInterpolatedPeaks spectrum2(pos, intensity);

START_SECTION(double getPosMin() const)
  TEST_EQUAL(spectrum2.getPosMin(), 416.3);
END_SECTION

START_SECTION(double getPosMax() const)
  TEST_EQUAL(spectrum2.getPosMax(), 419.2);
END_SECTION

START_SECTION(size_t size() const)
  TEST_EQUAL(spectrum2.size(), 2)
END_SECTION

START_SECTION(SplineInterpolatedPeaks::Navigator getNavigator(double scaling))
  // just to test if it can be called
  SplineInterpolatedPeaks::Navigator nav = spectrum2.getNavigator();
END_SECTION

START_SECTION(double SplineInterpolatedPeaks::Navigator::eval(double pos))
  // outside range of Gaussians
  TEST_EQUAL(spectrum2.getNavigator().eval(400.0), 0);
  TEST_EQUAL(spectrum2.getNavigator().eval(417.8), 0);
  TEST_EQUAL(spectrum2.getNavigator().eval(500.0), 0);
  // near the edge
  TEST_REAL_SIMILAR(spectrum2.getNavigator().eval(416.33), 0.007848195698809);  // expected 0.00738068453767004 differs by 6%
  // near the maximum
  TEST_REAL_SIMILAR(spectrum2.getNavigator().eval(416.81), 0.997572728799559);  // expected 0.99778024508561 differs by 0.02%
  // evaluation in first package, then search in last package
  SplineInterpolatedPeaks::Navigator nav = spectrum2.getNavigator();
  TEST_REAL_SIMILAR(nav.eval(416.81), 0.997572728799559);
  TEST_REAL_SIMILAR(nav.eval(418.75), 0.944147611428987);
  // evaluation in last package, then search in first package
  SplineInterpolatedPeaks::Navigator nav2 = spectrum2.getNavigator();
  TEST_REAL_SIMILAR(nav2.eval(418.75), 0.944147611428987);
  TEST_REAL_SIMILAR(nav2.eval(416.81), 0.997572728799559);
END_SECTION

START_SECTION(double SplineInterpolatedPeaks::Navigator::getNextPos(double pos))
  // advancing within package
  TEST_EQUAL(spectrum2.getNavigator().getNextPos(417.0), 417.07);
  // advancing to next package
  TEST_EQUAL(spectrum2.getNavigator().getNextPos(417.29), 418.2);
  // advancing beyond range
  TEST_REAL_SIMILAR(spectrum2.getNavigator().getNextPos(500.0), 419.2);
END_SECTION

// Each SplinePackage in a SplineInterpolatedPeaks must contain two or more data points. If this is not the case, the interpolation might lead to unexpected results.
// In the example below, a single data point @ 407.5 is placed between two packages. It does not form a SplinePackage on its own, but is instead part of the second SplinePackage.
std::vector<double> pos3;
std::vector<double> intensity3;
for (size_t i=0; i<4; ++i)
{
    pos3.push_back(400+i*0.5);
    intensity3.push_back(10.0);
}
pos3.push_back(407.5);
intensity3.push_back(10.0);
for (size_t i=0; i<4; ++i)
{
    pos3.push_back(410+i*0.5);
    intensity3.push_back(10.0);
}
SplineInterpolatedPeaks spectrum3(pos3, intensity3);

START_SECTION(double SplineInterpolatedPeaks::Navigator::eval(double pos))
  TEST_EQUAL(spectrum3.size(),2);
  TEST_EQUAL(spectrum3.getNavigator().eval(405),0);    // Zero as expected, since 405 is between packages.
  TEST_EQUAL(spectrum3.getNavigator().eval(408),10);    // One might expect zero, but 407.5 is part of the second package.
END_SECTION

std::vector<double> pos4;
std::vector<double> intensity4;
pos4.push_back(407.5);
intensity4.push_back(10.0);
START_SECTION(SplineInterpolatedPeaks(const std::vector<double>& pos, const std::vector<double>& intensity))
  TEST_EXCEPTION(Exception::IllegalArgument, new SplineInterpolatedPeaks(pos4,intensity4));
END_SECTION

// A navigator starts its search at the package it used last. Its results must not depend on that, i.e. on the
// positions requested before: on many packages, every position (left of, between, on the borders of, inside and right
// of the packages) gives the same as with a fresh navigator, requested in increasing, decreasing and random order.
std::vector<double> pos5;
std::vector<double> intensity5;
for (int package = 0; package < 40; ++package)
{
  const double start = 400.0 + package * 2.0 + 0.1 * (package % 3);
  for (int i = 0; i < 6 + package % 4; ++i)
  {
    pos5.push_back(start + 0.1 * i);
    intensity5.push_back(1.0 + (i % 3) + 0.01 * package);
  }
}
SplineInterpolatedPeaks spectrum5(pos5, intensity5);

START_SECTION([EXTRA] Navigator results do not depend on the positions requested before)
  TEST_EQUAL(spectrum5.size(), 40)
  std::vector<double> queries;
  for (double q = 398.0; q < spectrum5.getPosMax() + 2.0; q += 0.0137)
  {
    queries.push_back(q);
  }
  queries.insert(queries.end(), pos5.begin(), pos5.end()); // borders and data points
  std::vector<double> expected_eval, expected_next;
  for (double q : queries)
  {
    expected_eval.push_back(spectrum5.getNavigator().eval(q));
    expected_next.push_back(spectrum5.getNavigator().getNextPos(q));
  }
  std::vector<size_t> order(queries.size());
  for (size_t i = 0; i < order.size(); ++i) { order[i] = i; }
  std::vector<std::vector<size_t>> orders{order, std::vector<size_t>(order.rbegin(), order.rend())};
  std::mt19937 rng(42);
  std::shuffle(order.begin(), order.end(), rng);
  orders.push_back(order);
  for (const std::vector<size_t>& o : orders)
  {
    SplineInterpolatedPeaks::Navigator nav = spectrum5.getNavigator();
    Size mismatches = 0;
    for (size_t i : o)
    {
      mismatches += (nav.eval(queries[i]) != expected_eval[i]);
      mismatches += (nav.getNextPos(queries[i]) != expected_next[i]);
    }
    TEST_EQUAL(mismatches, 0)
  }
  // stepping through all packages, as the profile filter of FeatureFinderMultiplex does, visits every package
  SplineInterpolatedPeaks::Navigator nav = spectrum5.getNavigator();
  Size steps = 0, not_increasing = 0;
  double last = 0.0;
  for (double q = 398.0; q < spectrum5.getPosMax(); q = nav.getNextPos(q))
  {
    not_increasing += !(q > last);
    last = q;
    ++steps;
  }
  TEST_EQUAL(not_increasing, 0)
  TEST_EQUAL(steps > pos5.size(), true)
END_SECTION

END_TEST
