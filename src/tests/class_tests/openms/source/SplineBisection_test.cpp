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
#include <OpenMS/MATH/MISC/SplineBisection.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace
{
  // A minimal spline-like model exposing the eval()/derivative() interface that
  // Math::spline_bisection requires: a downward parabola
  //   f(x) = -a*(x - peak)^2 + height
  // with a single maximum of value `height` at x == `peak`.
  struct ParabolaSpline
  {
    double peak;
    double height;
    double a;
    double eval(double x) const { return -a * (x - peak) * (x - peak) + height; }
    double derivative(double x) const { return -2.0 * a * (x - peak); }
  };
}

START_TEST(SplineBisection, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((template <class T> void spline_bisection(const T& peak_spline, double const left_neighbor_mz, double const right_neighbor_mz, double& max_peak_mz, double& max_peak_int, double const threshold = 1e-6)))
{
  TOLERANCE_ABSOLUTE(1e-4)

  // a peak centered in the bracket -> bisection finds the apex (mz and intensity)
  {
    ParabolaSpline s{500.0, 1000.0, 1.0};
    double mz = 0.0, intensity = 0.0;
    Math::spline_bisection(s, 499.0, 501.0, mz, intensity);
    TEST_REAL_SIMILAR(mz, 500.0)
    TEST_REAL_SIMILAR(intensity, 1000.0)
  }

  // an off-center peak (and different curvature/height) is still located
  {
    ParabolaSpline s{500.3, 750.0, 2.0};
    double mz = 0.0, intensity = 0.0;
    Math::spline_bisection(s, 499.0, 501.0, mz, intensity);
    TEST_REAL_SIMILAR(mz, 500.3)
    TEST_REAL_SIMILAR(intensity, 750.0)
  }

  // a tighter threshold yields a more precise apex
  {
    ParabolaSpline s{500.123456, 100.0, 5.0};
    double mz = 0.0, intensity = 0.0;
    Math::spline_bisection(s, 499.0, 501.0, mz, intensity, 1e-9);
    TEST_REAL_SIMILAR(mz, 500.123456)
    TEST_REAL_SIMILAR(intensity, 100.0)
  }

  // if the true apex lies outside the bracket (here to the right of it), the
  // derivative never changes sign and the result converges to the right border
  {
    ParabolaSpline s{505.0, 50.0, 1.0}; // apex at 505, outside [499, 501]
    double mz = 0.0, intensity = 0.0;
    Math::spline_bisection(s, 499.0, 501.0, mz, intensity);
    TEST_REAL_SIMILAR(mz, 501.0)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
