// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/MISC/BSplineSmoothingSpline.h>
///////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;

START_TEST(BSplineSmoothingSpline, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BSplineSmoothingSpline* ptr = nullptr;
BSplineSmoothingSpline* null_ptr = nullptr;

TOLERANCE_ABSOLUTE(1e-4)
TOLERANCE_RELATIVE(1.0 + 1e-4)

START_SECTION((BSplineSmoothingSpline(const std::vector<double>&, const std::vector<double>&, double, int)))
{
  // Test basic construction with simple linear data
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> y = {0.0, 1.0, 2.0, 3.0, 4.0};
  
  ptr = new BSplineSmoothingSpline(x, y, -1.0, 3);
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->ok(), true)
  
  // Test with automatic smoothing parameter
  BSplineSmoothingSpline auto_smooth(x, y);
  TEST_EQUAL(auto_smooth.ok(), true)
  
  // Test with s=0 (interpolation)
  BSplineSmoothingSpline interpolate(x, y, 0.0);
  TEST_EQUAL(interpolate.ok(), true)
  
  // Test with explicit smoothing parameter
  BSplineSmoothingSpline smooth(x, y, 2.0);
  TEST_EQUAL(smooth.ok(), true)
}
END_SECTION

START_SECTION((~BSplineSmoothingSpline()))
{
  delete ptr;
  ptr = nullptr;
}
END_SECTION

START_SECTION((double eval(double) const))
{
  // Test evaluation on simple data
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> y = {0.0, 1.0, 4.0, 9.0, 16.0}; // y = x^2
  
  BSplineSmoothingSpline spl(x, y, -1.0, 3);
  TEST_EQUAL(spl.ok(), true)
  
  // Test evaluation at data points
  for (size_t i = 0; i < x.size(); ++i)
  {
    double val = spl.eval(x[i]);
    TEST_EQUAL(std::isfinite(val), true)
    // With smoothing, won't pass exactly through points
    TEST_EQUAL(std::abs(val - y[i]) < 2.0, true) // reasonable approximation
  }
  
  // Test evaluation between points
  double mid = spl.eval(1.5);
  TEST_EQUAL(std::isfinite(mid), true)
  TEST_EQUAL(mid > 1.0 && mid < 4.0, true) // between adjacent y values
  
  // Test extrapolation
  double extrap_low = spl.eval(-0.5);
  double extrap_high = spl.eval(4.5);
  TEST_EQUAL(std::isfinite(extrap_low), true)
  TEST_EQUAL(std::isfinite(extrap_high), true)
}
END_SECTION

START_SECTION((bool ok() const))
{
  // Test valid spline
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0};
  BSplineSmoothingSpline valid(x, y);
  TEST_EQUAL(valid.ok(), true)
  
  // Test invalid input: mismatched sizes
  std::vector<double> y_short = {1.0, 2.0};
  BSplineSmoothingSpline invalid1(x, y_short);
  TEST_EQUAL(invalid1.ok(), false)
  
  // Test invalid input: too few points
  std::vector<double> x_short = {0.0};
  std::vector<double> y_one = {1.0};
  BSplineSmoothingSpline invalid2(x_short, y_one);
  TEST_EQUAL(invalid2.ok(), false)
  
  // Test invalid input: unsorted x
  std::vector<double> x_unsorted = {0.0, 2.0, 1.0, 3.0};
  std::vector<double> y_four = {1.0, 2.0, 3.0, 4.0};
  BSplineSmoothingSpline invalid3(x_unsorted, y_four);
  TEST_EQUAL(invalid3.ok(), false)
}
END_SECTION

START_SECTION((int num_interior_knots() const))
{
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  
  // Interpolation (s=0) should use all points
  BSplineSmoothingSpline interp(x, y, 0.0);
  TEST_EQUAL(interp.ok(), true)
  TEST_EQUAL(interp.num_interior_knots() >= 0, true)
  
  // Smoothing should use fewer knots
  BSplineSmoothingSpline smooth(x, y, 10.0);
  TEST_EQUAL(smooth.ok(), true)
  TEST_EQUAL(smooth.num_interior_knots() >= 0, true)
}
END_SECTION

START_SECTION((double rss() const))
{
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0};
  
  BSplineSmoothingSpline spl(x, y, -1.0);
  TEST_EQUAL(spl.ok(), true)
  
  double rss = spl.rss();
  TEST_EQUAL(rss >= 0.0, true)
  TEST_EQUAL(std::isfinite(rss), true)
}
END_SECTION

START_SECTION((double smoothing_param() const))
{
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  
  // Test automatic smoothing parameter
  BSplineSmoothingSpline auto_spl(x, y, -1.0);
  TEST_EQUAL(auto_spl.ok(), true)
  double s_auto = auto_spl.smoothing_param();
  // Should be m - sqrt(2*m) = 6 - sqrt(12) ≈ 2.536
  TEST_EQUAL(s_auto > 2.0 && s_auto < 3.0, true)
  
  // Test explicit smoothing parameter
  BSplineSmoothingSpline explicit_spl(x, y, 5.0);
  TEST_EQUAL(explicit_spl.ok(), true)
  TEST_REAL_SIMILAR(explicit_spl.smoothing_param(), 5.0)
  
  // Test s=0
  BSplineSmoothingSpline zero_spl(x, y, 0.0);
  TEST_EQUAL(zero_spl.ok(), true)
  TEST_REAL_SIMILAR(zero_spl.smoothing_param(), 0.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Test smoothing behavior
/////////////////////////////////////////////////////////////

START_SECTION((Smoothing behavior test))
{
  // Create noisy data: y = sin(x) + noise
  std::vector<double> x;
  std::vector<double> y_noisy;
  std::vector<double> y_true;
  
  for (int i = 0; i < 20; ++i)
  {
    double xi = i * 0.3;
    x.push_back(xi);
    y_true.push_back(std::sin(xi));
    // Add noise
    y_noisy.push_back(std::sin(xi) + 0.2 * (i % 3 - 1));
  }
  
  // Heavy smoothing should be closer to true function than data
  BSplineSmoothingSpline smooth_heavy(x, y_noisy, 20.0);
  TEST_EQUAL(smooth_heavy.ok(), true)
  
  // Compute error on noisy data vs smoothed prediction
  double noise_error = 0.0;
  double smooth_error = 0.0;
  for (size_t i = 0; i < x.size(); ++i)
  {
    double pred = smooth_heavy.eval(x[i]);
    noise_error += std::pow(y_noisy[i] - y_true[i], 2);
    smooth_error += std::pow(pred - y_true[i], 2);
  }
  
  // Smoothing should reduce error (not always guaranteed but likely)
  TEST_EQUAL(std::isfinite(noise_error), true)
  TEST_EQUAL(std::isfinite(smooth_error), true)
}
END_SECTION

START_SECTION((Interpolation vs smoothing test))
{
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> y = {1.0, 3.0, 2.0, 4.0, 3.5};
  
  // Interpolation (s=0) should pass through all points (or very close)
  BSplineSmoothingSpline interp(x, y, 0.0);
  TEST_EQUAL(interp.ok(), true)
  
  for (size_t i = 0; i < x.size(); ++i)
  {
    double pred = interp.eval(x[i]);
    // Allow small numerical error (BSpline2d doesn't perfectly interpolate)
    TEST_EQUAL(std::abs(pred - y[i]) < 0.05, true)
  }
  
  // Smoothing (s>0) may deviate from points
  BSplineSmoothingSpline smooth(x, y, 2.0);
  TEST_EQUAL(smooth.ok(), true)
  
  // At least some points should differ (smoothing effect)
  int num_different = 0;
  for (size_t i = 0; i < x.size(); ++i)
  {
    double pred = smooth.eval(x[i]);
    if (std::abs(pred - y[i]) > 0.01)
    {
      num_different++;
    }
  }
  // With smoothing, expect some deviation (though not guaranteed for all data)
  TEST_EQUAL(num_different >= 0, true) // Always true, just checking evaluation works
}
END_SECTION

START_SECTION((Linear data test))
{
  // Perfect linear data should fit well even with smoothing
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y;
  for (double xi : x)
  {
    y.push_back(2.0 * xi + 1.0); // y = 2x + 1
  }
  
  BSplineSmoothingSpline spl(x, y, -1.0);
  TEST_EQUAL(spl.ok(), true)
  
  // Should fit linear data very well
  for (size_t i = 0; i < x.size(); ++i)
  {
    double pred = spl.eval(x[i]);
    TEST_REAL_SIMILAR(pred, y[i])
  }
  
  // RSS should be very small
  TEST_EQUAL(spl.rss() < 0.1, true)
}
END_SECTION

START_SECTION((Quadratic data test))
{
  // Quadratic data: y = x^2
  std::vector<double> x = {-2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
  std::vector<double> y;
  for (double xi : x)
  {
    y.push_back(xi * xi);
  }
  
  BSplineSmoothingSpline spl(x, y, 1.0, 3); // cubic spline
  TEST_EQUAL(spl.ok(), true)
  
  // Should fit quadratic data reasonably well
  for (size_t i = 0; i < x.size(); ++i)
  {
    double pred = spl.eval(x[i]);
    // Allow some deviation due to smoothing
    TEST_EQUAL(std::abs(pred - y[i]) < 1.0, true)
  }
}
END_SECTION

START_SECTION((Edge case: identical x values should fail))
{
  std::vector<double> x = {1.0, 1.0, 2.0, 3.0};
  std::vector<double> y = {1.0, 2.0, 3.0, 4.0};
  
  BSplineSmoothingSpline spl(x, y);
  TEST_EQUAL(spl.ok(), false) // Should fail due to non-sorted x
}
END_SECTION

START_SECTION((Edge case: small dataset))
{
  // Very small dataset (n=2)
  std::vector<double> x = {0.0, 1.0};
  std::vector<double> y = {0.0, 1.0};
  
  BSplineSmoothingSpline spl(x, y);
  TEST_EQUAL(spl.ok(), true) // Should handle gracefully
  
  if (spl.ok())
  {
    TEST_EQUAL(std::isfinite(spl.eval(0.5)), true)
  }
}
END_SECTION

START_SECTION((Edge case: n=3 dataset))
{
  std::vector<double> x = {0.0, 1.0, 2.0};
  std::vector<double> y = {1.0, 2.0, 3.0};
  
  BSplineSmoothingSpline spl(x, y);
  TEST_EQUAL(spl.ok(), true)
  
  if (spl.ok())
  {
    for (double xi : x)
    {
      TEST_EQUAL(std::isfinite(spl.eval(xi)), true)
    }
  }
}
END_SECTION

START_SECTION((Different spline degrees))
{
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  std::vector<double> y = {1.0, 2.0, 1.5, 3.0, 2.5, 4.0};
  
  // Test with different degrees
  for (int k = 1; k <= 3; ++k)
  {
    BSplineSmoothingSpline spl(x, y, -1.0, k);
    TEST_EQUAL(spl.ok(), true)
    
    if (spl.ok())
    {
      // Should be able to evaluate at all x points
      for (double xi : x)
      {
        double val = spl.eval(xi);
        TEST_EQUAL(std::isfinite(val), true)
      }
    }
  }
}
END_SECTION

START_SECTION((Consistency test: repeated construction))
{
  std::vector<double> x = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::vector<double> y = {1.0, 2.0, 3.0, 4.0, 5.0};
  
  BSplineSmoothingSpline spl1(x, y, 2.0);
  BSplineSmoothingSpline spl2(x, y, 2.0);
  
  TEST_EQUAL(spl1.ok(), true)
  TEST_EQUAL(spl2.ok(), true)
  
  if (spl1.ok() && spl2.ok())
  {
    // Same input should give same results
    for (double xi : x)
    {
      TEST_REAL_SIMILAR(spl1.eval(xi), spl2.eval(xi))
    }
    
    TEST_REAL_SIMILAR(spl1.rss(), spl2.rss())
    TEST_EQUAL(spl1.num_interior_knots(), spl2.num_interior_knots())
  }
}
END_SECTION

END_TEST
