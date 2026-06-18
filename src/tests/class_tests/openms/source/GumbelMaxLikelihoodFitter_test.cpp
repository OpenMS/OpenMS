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

#include <cmath>
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
