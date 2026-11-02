// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/STATISTICS/GumbelMaxLikelihoodFitter.h>
///////////////////////////

#include <OpenMS/CONCEPT/Exception.h>

#include <cmath>
#include <limits>
#include <vector>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

using FitResult = GumbelMaxLikelihoodFitter::GumbelDistributionFitResult;

START_TEST(GumbelMaxLikelihoodFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GumbelMaxLikelihoodFitter* ptr = nullptr;
GumbelMaxLikelihoodFitter* null_ptr = nullptr;

START_SECTION((GumbelMaxLikelihoodFitter()))
{
  ptr = new GumbelMaxLikelihoodFitter();
  TEST_NOT_EQUAL(ptr, null_ptr)
  // with no data points, fitWeighted returns the (default) initial parameters
  FitResult r = ptr->fitWeighted({}, {});
  TEST_REAL_SIMILAR(r.a, 0.25)
  TEST_REAL_SIMILAR(r.b, 0.1)
}
END_SECTION

START_SECTION((~GumbelMaxLikelihoodFitter()))
{
  delete ptr;
}
END_SECTION

START_SECTION((GumbelMaxLikelihoodFitter(GumbelDistributionFitResult init)))
{
  GumbelMaxLikelihoodFitter f(FitResult(3.0, 0.5));
  // empty data -> the supplied initial parameters are returned unchanged
  FitResult r = f.fitWeighted({}, {});
  TEST_REAL_SIMILAR(r.a, 3.0)
  TEST_REAL_SIMILAR(r.b, 0.5)
}
END_SECTION

START_SECTION((void setInitialParameters(const GumbelDistributionFitResult & result)))
{
  GumbelMaxLikelihoodFitter f;
  f.setInitialParameters(FitResult(4.0, 0.7));
  FitResult r = f.fitWeighted({}, {});
  TEST_REAL_SIMILAR(r.a, 4.0)
  TEST_REAL_SIMILAR(r.b, 0.7)
}
END_SECTION

START_SECTION((GumbelDistributionFitResult fitWeighted(const std::vector<double> & x, const std::vector<double> & w)))
{
  // Build a Gumbel(a=2.0, b=0.8) probability density on a grid and use the
  // densities as weights -> the MLE fit must recover the generating parameters.
  const double A = 2.0, B = 0.8;
  std::vector<double> x, w;
  for (double xi = -2.0; xi <= 8.0; xi += 0.1)
  {
    const double z = (xi - A) / B;
    const double dens = (1.0 / B) * std::exp(-(z + std::exp(-z)));
    x.push_back(xi);
    w.push_back(dens);
  }

  GumbelMaxLikelihoodFitter f(FitResult(1.0, 1.0));
  FitResult r = f.fitWeighted(x, w);

  TOLERANCE_ABSOLUTE(0.05)
  TEST_REAL_SIMILAR(r.a, A)
  TEST_REAL_SIMILAR(r.b, B)
  TEST_EQUAL(r.b > 0.0, true)

  // Narrow distributions (b < ~0.2) have a positive log-likelihood at the optimum.
  // A least-squares minimization of the log-likelihood would converge to LL = 0 instead of the MLE.
  {
    const double A2 = -1.5, B2 = 0.05;
    std::vector<double> x2, w2;
    for (double xi = -2.0; xi <= 0.0; xi += 0.002)
    {
      const double z = (xi - A2) / B2;
      x2.push_back(xi);
      w2.push_back((1.0 / B2) * std::exp(-(z + std::exp(-z))));
    }
    GumbelMaxLikelihoodFitter f2(FitResult(0.0, 1.0)); // deliberately poor start
    FitResult r2 = f2.fitWeighted(x2, w2);
    TOLERANCE_ABSOLUTE(0.002)
    TEST_REAL_SIMILAR(r2.a, A2)
    TEST_REAL_SIMILAR(r2.b, B2)
  }

  // the MLE must satisfy the closed-form location equation: sum w * exp(-(x - a)/b) = sum w
  {
    std::vector<double> x3 = {0.1, 0.5, 0.7, 1.3, 2.9, 0.4, 0.8, 1.1};
    std::vector<double> w3 = {1.0, 2.0, 1.0, 0.5, 0.25, 1.0, 1.0, 3.0};
    GumbelMaxLikelihoodFitter f3;
    FitResult r3 = f3.fitWeighted(x3, w3);
    double lhs = 0.0, sum_w = 0.0;
    for (size_t i = 0; i < x3.size(); ++i)
    {
      lhs += w3[i] * std::exp(-(x3[i] - r3.a) / r3.b);
      sum_w += w3[i];
    }
    TOLERANCE_ABSOLUTE(1e-8)
    TEST_REAL_SIMILAR(lhs, sum_w)
    // the fitted parameters are the new initial parameters
    FitResult r3b = f3.fitWeighted({}, {});
    TEST_REAL_SIMILAR(r3b.a, r3.a)
    TEST_REAL_SIMILAR(r3b.b, r3.b)
  }

  // non-finite values and non-positive weights are ignored
  {
    std::vector<double> x4 = {0.1, 0.5, 0.7, 1.3, 2.9, std::numeric_limits<double>::infinity(), 100.0};
    std::vector<double> w4 = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0};
    GumbelMaxLikelihoodFitter f4;
    FitResult r4 = f4.fitWeighted(x4, w4);
    GumbelMaxLikelihoodFitter f5;
    FitResult r5 = f5.fitWeighted({0.1, 0.5, 0.7, 1.3, 2.9}, {1.0, 1.0, 1.0, 1.0, 1.0});
    TOLERANCE_ABSOLUTE(1e-10)
    TEST_REAL_SIMILAR(r4.a, r5.a)
    TEST_REAL_SIMILAR(r4.b, r5.b)
  }

  // data without spread cannot be fitted; mismatching sizes are rejected
  {
    GumbelMaxLikelihoodFitter f6;
    TEST_EXCEPTION(Exception::UnableToFit, f6.fitWeighted({1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}))
    TEST_EXCEPTION(Exception::InvalidParameter, f6.fitWeighted({1.0, 2.0}, {1.0}))
  }
}
END_SECTION

START_SECTION(([EXTRA] GumbelDistributionFitResult::eval and mode))
{
  FitResult r(2.0, 0.5);
  TEST_REAL_SIMILAR(r.mode(), 2.0)
  TEST_REAL_SIMILAR(r.eval(2.0), std::exp(-1.0) / 0.5)
  TEST_REAL_SIMILAR(std::log(r.eval(3.1)), r.log_eval_no_normalize(3.1))
}
END_SECTION

START_SECTION(([EXTRA] GumbelDistributionFitResult(double a, double b) + log_eval_no_normalize))
{
  // the struct stores a and b
  FitResult r(2.0, 1.0);
  TEST_REAL_SIMILAR(r.a, 2.0)
  TEST_REAL_SIMILAR(r.b, 1.0)

  // log Gumbel density: -log(b) - (x-a)/b - exp(-(x-a)/b)
  auto logGumbel = [](double x, double a, double b) {
    const double diff = (x - a) / b;
    return -std::log(b) - diff - std::exp(-diff);
  };
  TEST_REAL_SIMILAR(r.log_eval_no_normalize(2.0), logGumbel(2.0, 2.0, 1.0)) // = -1
  TEST_REAL_SIMILAR(r.log_eval_no_normalize(3.0), logGumbel(3.0, 2.0, 1.0))
  TEST_REAL_SIMILAR(r.log_eval_no_normalize(0.0), logGumbel(0.0, 2.0, 1.0))

  // the maximum of the log-density is at x = a (mode of the Gumbel)
  TEST_EQUAL(r.log_eval_no_normalize(2.0) > r.log_eval_no_normalize(3.0), true)
  TEST_EQUAL(r.log_eval_no_normalize(2.0) > r.log_eval_no_normalize(1.0), true)

  // scale parameter shifts the density by -log(b)
  FitResult r2(2.0, 2.0);
  TEST_REAL_SIMILAR(r2.log_eval_no_normalize(2.0), logGumbel(2.0, 2.0, 2.0))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
