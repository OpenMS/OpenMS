// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ML/RANSAC/RANSAC.h>
#include <OpenMS/ML/RANSAC/RANSACModelQuadratic.h>
///////////////////////////

using namespace std;
using namespace OpenMS;
using namespace Math;

///////////////////////////

// random number generator using srand (used in std::random_shuffle())
int myRNG(int n) {
  return std::rand() / (1.0 + RAND_MAX) * n;
}

START_TEST(RANSACModelQuadratic, "$Id$")

// fixed seed across all platforms
srand(123);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RansacModelQuadratic mod;
/*
R code to produce the test data:
A = 15
B = 1.04
C = 0.00001

ppmSD = 2
p_outlier = 0.0

quad = function(xi, withError = 1, p_outlier = p_outlier, A=A, B=B,C=C)
{
errPPM = rnorm(1, 0, ppmSD)
errDa = errPPM / 1e6 *xi
yi = A + B*xi + C*xi*xi + errDa * withError
## make it an outlier (up to x10 times off)?
if (p_outlier > runif(1)) yi = yi * runif(1, 0.1, 10) ## upscaling bias, but hey...
yi
}

x = seq(300, 1200, by = 50)
y = sapply(x, quad)
yperfect = sapply(x, quad, withError=0, p_outlier = 0.0, A= 15.000922133127460, B=1.0399979867717661, C=1.0000728397741021e-005)

plot(x, y-yperfect)

plot(x,y)
paste(x, collapse=", ")
paste(y, collapse=", ")

*/

double tx[] = {300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1150, 1200};
double ty[] = {327.899691196695, 380.224568059509, 432.60001240954, 485.025329595228, 537.500950780092, 590.025780790517, 642.599328692568, 695.226226189721, 747.901404527071, 800.628123375242, 853.398402241557, 906.221690051125, 959.098650822155, 1012.02444287071, 1064.99750095039, 1118.02244136503, 1171.10218429403, 1224.22697636988, 1277.3989501643};
// with outliers
double tyo[] = {327.899129765285, 380.22505352209, 3107.92239745832, 485.025787154647, 4105.48012991713, 590.022752890109, 642.600208487572, 695.225450449224, 747.90063358956, 800.627835435032, 853.39930288619, 906.221242296453, 4676.76961970711, 1012.02448829764, 1765.7626054147, 1118.02258525226, 1171.09984462503, 1224.22603819705, 1277.40083805864};

std::vector<std::pair<double, double> > test_pairs;
test_pairs.reserve(19);
for (int i = 0; i < 19; ++i)
{
  test_pairs.push_back(make_pair(tx[i], ty[i]));
}

std::vector<std::pair<double, double> > test_pairs_o;
test_pairs_o.reserve(19);
for (int i = 0; i < 19; ++i)
{
  test_pairs_o.push_back(make_pair(tx[i], tyo[i]));
}

START_SECTION((static ModelParameters rm_fit_impl(const DVecIt& begin, const DVecIt& end)))
{


  RansacModel<>::ModelParameters coeff = mod.rm_fit_impl(test_pairs.begin(), test_pairs.end());
  TEST_REAL_SIMILAR( coeff[0], 15.0009) // should be 15.0
  TEST_REAL_SIMILAR( coeff[1], 1.04)
  TEST_REAL_SIMILAR( coeff[2], 0.00001)

  double rss = mod.rm_rss_impl(test_pairs.begin(), test_pairs.end(), coeff);
  TEST_REAL_SIMILAR( rss, 5.2254915523925468e-005)

  RansacModel<>::DVec inliers = mod.rm_inliers(test_pairs.begin(), test_pairs.end(), coeff, 0.5);
  TEST_EQUAL(inliers.size(), test_pairs.size()) // all should be inliers

  inliers = mod.rm_inliers(test_pairs_o.begin(), test_pairs_o.end(), coeff, 0.5);
  TEST_EQUAL( inliers.size(), 15); // 19-15 = 4 outliers
  // just test the gaps
  TEST_REAL_SIMILAR( inliers[2].first, 450.0)
  TEST_REAL_SIMILAR( inliers[3].first, 550.0)
  TEST_REAL_SIMILAR( inliers[10].first, 950.0)
  TEST_REAL_SIMILAR( inliers[11].first, 1050)
}
END_SECTION

START_SECTION((static double rm_rsq_impl(const DVecIt& begin, const DVecIt& end)))
  NOT_TESTABLE  // tested above in rm_fit_impl
END_SECTION

START_SECTION((static double rm_rss_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients)))
  NOT_TESTABLE // tested above in rm_fit_impl
END_SECTION

START_SECTION((static DVec rm_inliers_impl(const DVecIt& begin, const DVecIt& end, const ModelParameters& coefficients, double max_threshold)))
  NOT_TESTABLE // tested above in rm_fit_impl
END_SECTION

START_SECTION([EXTRA](static Math::RANSAC<Math::RANSACModelQuadratic>::ransac(const std::vector<std::pair<double, double> >& pairs, 
                                                                           size_t n, 
                                                                           size_t k, 
                                                                           double t, 
                                                                           size_t d, 
                                                                           bool relative_d = false,
                                                                           int (*rng)(int) = NULL)))
{
  // full RANSAC with outliers
  /*@param n The minimum number of data points required to fit the model
    @param k The maximum number of iterations allowed in the algorithm 
    @param t Threshold value for determining when a data point fits a
             model. Corresponds to the maximal squared deviation in units of the
             _second_ dimension (dim2).
    @param d The number of close data values (according to 't') required to assert that a model fits well to data
    @param rng Custom RNG function (useful for testing with fixed seeds)
  */
  RANSAC<RansacModelQuadratic> r{0};
  std::vector<std::pair<double, double> > test_pairs_out = r.ransac(test_pairs_o, 5, 50, 2.0, 3, false);

  TEST_EQUAL( test_pairs_out.size(), 15)
  ABORT_IF(test_pairs_out.size() != 15)
  // just test the gaps
  std::sort(test_pairs_out.begin(), test_pairs_out.end());
  TEST_REAL_SIMILAR( test_pairs_out[2].first, 450.0)
  TEST_REAL_SIMILAR( test_pairs_out[3].first, 550.0)
  TEST_REAL_SIMILAR( test_pairs_out[10].first, 950.0)
  TEST_REAL_SIMILAR( test_pairs_out[11].first, 1050)

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

