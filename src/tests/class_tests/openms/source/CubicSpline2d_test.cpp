// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/MATH/MISC/CubicSpline2d.h>

using namespace OpenMS;

START_TEST(CubicSpline2d, "$Id$")

std::vector<double> mz;
mz.push_back(486.784);
mz.push_back(486.787);
mz.push_back(486.790);
mz.push_back(486.793);
mz.push_back(486.795);
mz.push_back(486.797);
mz.push_back(486.800);
mz.push_back(486.802);
mz.push_back(486.805);
mz.push_back(486.808);
mz.push_back(486.811);
std::vector<double> intensity;
intensity.push_back(0.0);
intensity.push_back(154683.17);
intensity.push_back(620386.5);
intensity.push_back(1701390.12);
intensity.push_back(2848879.25);
intensity.push_back(3564045.5);
intensity.push_back(2744585.7);
intensity.push_back(1605583.0);
intensity.push_back(1518984.0);
intensity.push_back(1591352.21);
intensity.push_back(1691345.1);

double x_min = -0.5;
double x_max = 1.5;
Size n = 10;
std::vector<double> x;
std::vector<double> y;
for (Size i=0; i<(n+1); ++i)
{
    x.push_back(x_min + (double)i/10*(x_max-x_min));
    y.push_back(sin(x_min + (double)i/10*(x_max-x_min)));
}

std::map<double,double> mz_intensity;
for (Size i=0; i<mz.size(); ++i)
{
    mz_intensity.insert(std::pair<double,double>(mz[i], intensity[i]));
}

std::map<double,double> x_y;
for (Size i=0; i<x.size(); ++i)
{
    x_y.insert(std::pair<double,double>(x[i], y[i]));
}

CubicSpline2d sp1(mz, intensity);
CubicSpline2d sp2(mz_intensity);
CubicSpline2d sp5(x, y);
CubicSpline2d sp6(x_y);

CubicSpline2d* nullPointer = nullptr;

START_SECTION(CubicSpline2d(const std::vector<double>& x, const std::vector<double>& y))
  CubicSpline2d* sp3 = new CubicSpline2d(mz, intensity);
  TEST_NOT_EQUAL(sp3, nullPointer)
  delete sp3;
END_SECTION

START_SECTION(CubicSpline2d(const std::map<double, double>& m))
  CubicSpline2d* sp4 = new CubicSpline2d(mz_intensity);
  TEST_NOT_EQUAL(sp4, nullPointer)
  delete sp4;
END_SECTION

START_SECTION(double eval(double x))
  // near border of spline range
  TEST_REAL_SIMILAR(sp1.eval(486.785), 35173.1841778984);
  TEST_REAL_SIMILAR(sp2.eval(486.785), 35173.1841778984);
  // inside spline range
  TEST_REAL_SIMILAR(sp1.eval(486.794), 2271426.93316241);
  TEST_REAL_SIMILAR(sp2.eval(486.794), 2271426.93316241);
  // at the input nodes
  TEST_REAL_SIMILAR(sp1.eval(486.784), 0.0);
  TEST_REAL_SIMILAR(sp1.eval(486.790), 620386.5);
  TEST_REAL_SIMILAR(sp2.eval(486.808), 1591352.21);
  TEST_REAL_SIMILAR(sp2.eval(486.811), 1691345.1);
  // test sine at nodes
  for (Size i=0; i<(n+1); ++i)
  {
    TEST_REAL_SIMILAR(sp5.eval(x[i]), y[i]);
    TEST_REAL_SIMILAR(sp6.eval(x[i]), y[i]);
  }
  // test sine between nodes
  // The cubic spline is a third order approximation of the (co)sines.
  TOLERANCE_RELATIVE(1.005);
  for (Size i=0; i<(n+6); ++i)
  {
    double xx = x_min + (double)i/(n+5)*(x_max-x_min);
    TEST_REAL_SIMILAR(sp5.eval(xx), sin(xx));
    TEST_REAL_SIMILAR(sp6.eval(xx), sin(xx));
  }
END_SECTION

START_SECTION(double derivatives(double x, unsigned order))
  // near border of spline range
  TEST_REAL_SIMILAR(sp1.derivatives(486.785,1), 39270152.2996247)
  TEST_REAL_SIMILAR(sp1.derivatives(486.785,2), 12290904368.2736);
  TEST_REAL_SIMILAR(sp2.derivatives(486.785,1), 39270152.2996247);
  TEST_REAL_SIMILAR(sp2.derivatives(486.785,2), 12290904368.2736);
  // inside spline range
  TEST_REAL_SIMILAR(sp1.derivatives(486.794,1), 594825947.154264);
  TEST_REAL_SIMILAR(sp1.derivatives(486.794,2), 7415503644.8958);
  TEST_REAL_SIMILAR(sp2.derivatives(486.794,1), 594825947.154264);
  TEST_REAL_SIMILAR(sp2.derivatives(486.794,2), 7415503644.8958);
  // test cosine at nodes
  // No tests near boundaries, since deviation from cos(x) large and expected.
  TOLERANCE_RELATIVE(1.01);
  for (Size i=2; i<n-1; ++i)
  {
    TEST_REAL_SIMILAR(sp5.derivatives(x[i],1), cos(x[i]));
    TEST_REAL_SIMILAR(sp6.derivatives(x[i],1), cos(x[i]));
  }
  // test cosine between nodes
  for (Size i=2; i<(n+4); ++i)
  {
    double xx = x_min + (double)i/(n+5)*(x_max-x_min);
    TEST_REAL_SIMILAR(sp5.derivatives(xx,1), cos(xx));
    TEST_REAL_SIMILAR(sp6.derivatives(xx,1), cos(xx));
  }
  // test boundary conditions y"=0
  TEST_REAL_SIMILAR(sp5.derivatives(x[0],2), 0);
  TEST_REAL_SIMILAR(sp6.derivatives(x[0],2), 0);
  TEST_REAL_SIMILAR(sp5.derivatives(x[n],2), 0);
  TEST_REAL_SIMILAR(sp6.derivatives(x[n],2), 0);
END_SECTION

END_TEST
