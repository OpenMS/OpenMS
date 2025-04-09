// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/PROCESSING/SMOOTHING/FastLowessSmoothing.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

double targetFunction(double x)
{
  return 10 + 20*x + 40*x*x;
}

START_TEST(FastLowessSmoothing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION([FastLowessSmoothing_Original tests]void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&))
{

/*
* These are the original tests described in the FORTRAN code. We should be able to reproduce those.
*
*        X values:
*          1  2  3  4  5  (10)6  8  10  12  14  50
* 
*        Y values:
*           18  2  15  6  10  4  16  11  7  3  14  17  20  12  9  13  1  8  5  19
* 
* 
*        YS values with F = .25, NSTEPS = 0, DELTA = 0.0
*         13.659  11.145  8.701  9.722  10.000  (10)11.300  13.000  6.440  5.596
*           5.456  18.998
* 
*        YS values with F = .25, NSTEPS = 0 ,  DELTA = 3.0
*          13.659  12.347  11.034  9.722  10.511  (10)11.300  13.000  6.440  5.596
*            5.456  18.998
* 
*        YS values with F = .25, NSTEPS = 2, DELTA = 0.0
*          14.811  12.115  8.984  9.676  10.000  (10)11.346  13.000  6.734  5.744
*            5.415  18.998
*/

  double xval[] = {1, 2, 3, 4, 5, 6,  6,  6,  6,  6,  6,  6,  6,  6,  6, 8, 10, 12, 14, 50 };
  double yval[] = { 18, 2, 15, 6, 10, 4, 16, 11, 7, 3, 14, 17, 20, 12, 9, 13, 1, 8, 5, 19};

  double ys_1[] = { 13.659, 11.145, 8.701, 9.722, 10.000, 11.300, 11.300,
    11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 13.000,
    6.440, 5.596, 5.456, 18.998};

  double ys_2[] = {13.659, 12.347, 11.034, 9.722, 10.511, 11.300, 11.300,
    11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 11.300, 13.000,
    6.440, 5.596, 5.456, 18.998};

  double ys_3[] = { 14.811, 12.115, 8.984, 9.676, 10.000, 11.346, 11.346,
    11.346, 11.346, 11.346, 11.346, 11.346, 11.346, 11.346, 11.346, 13.000,
    6.734, 5.744, 5.415, 18.998};

  // the original test has limited numerical accuracy
  TOLERANCE_RELATIVE(1e-4);
  TOLERANCE_ABSOLUTE(1e-3);
  {
    std::vector< double > v_xval;
    std::vector< double > v_yval;
    for (size_t i = 0; i < 20; i++)
    {
      v_xval.push_back(xval[i]);
      v_yval.push_back(yval[i]);
    }

    // YS values with F = .25, NSTEPS = 0, DELTA = 0.0
    {
      std::vector< double > out(20), tmp1(20), tmp2(20);
      FastLowessSmoothing::lowess(v_xval, v_yval, 0.25, 0, 0.0, out);
      for (size_t i = 0; i < 20; i++)
      {
        TEST_REAL_SIMILAR(out[i], ys_1[i]);
      }
    }

    // YS values with F = .25, NSTEPS = 0 ,  DELTA = 3.0
    {
      std::vector< double > out(20), tmp1(20), tmp2(20);
      FastLowessSmoothing::lowess(v_xval, v_yval, 0.25, 0, 3.0, out);
      for (size_t i = 0; i < 20; i++)
      {
        TEST_REAL_SIMILAR(out[i], ys_2[i]);
      }
    }

    // YS values with F = .25, NSTEPS = 2, DELTA = 0.0
    {
      std::vector< double > out(20), tmp1(20), tmp2(20);
      FastLowessSmoothing::lowess(v_xval, v_yval, 0.25, 2, 0.0, out);
      for (size_t i = 0; i < 20; i++)
      {
        TEST_REAL_SIMILAR(out[i], ys_3[i]);
      }
    }
  }

}
END_SECTION


START_SECTION([FastLowessSmoothing_cars]void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&))
{

  /*
  In R

     require(graphics)
     
     plot(cars, main = "lowess(cars)")
     lines(lowess(cars), col = 2)
     lines(lowess(cars, f = .2), col = 3)
     legend(5, 120, c(paste("f = ", c("2/3", ".2"))), lty = 1, col = 2:3)

     The data below is what we expect from the R function when running the cars example.
     
  */
  TOLERANCE_RELATIVE(4e-7);

  int speed[] = {4, 4, 7, 7, 8, 9, 10, 10, 10, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14, 14, 15, 15, 15, 16, 16, 17, 17, 17, 18, 18, 18, 18, 19, 19, 19, 20, 20, 20, 20, 20, 22, 23, 24, 24, 24, 24, 25};
  int dist[] = { 2, 10, 4, 22, 16, 10, 18, 26, 34, 17, 28, 14, 20, 24, 28, 26, 34, 34, 46, 26, 36, 60, 80, 20, 26, 54, 32, 40, 32, 40, 50, 42, 56, 76, 84, 36, 46, 68, 32, 48, 52, 56, 64, 66, 54, 70, 92, 93, 120, 85};

  double expected_1[] = {4.965459, 4.965459, 13.124495, 13.124495, 15.858633, 18.579691, 21.280313, 21.280313, 21.280313, 24.129277, 24.129277, 27.119549, 27.119549, 27.119549, 27.119549, 30.027276, 30.027276, 30.027276, 30.027276, 32.962506, 32.962506, 32.962506, 32.962506, 36.757728, 36.757728, 36.757728, 40.435075, 40.435075, 43.463492, 43.463492, 43.463492, 46.885479, 46.885479, 46.885479, 46.885479, 50.793152, 50.793152, 50.793152, 56.491224, 56.491224, 56.491224, 56.491224, 56.491224, 67.585824, 73.079695, 78.643164, 78.643164, 78.643164, 78.643164, 84.328698};
  double expected_2[] = {6.030408, 6.030408, 12.678893, 12.678893, 15.383796, 18.668847, 22.227571, 22.227571, 22.227571, 23.306483, 23.306483, 21.525372, 21.525372, 21.525372, 21.525372, 34.882735, 34.882735, 34.882735, 34.882735, 47.059947, 47.059947, 47.059947, 47.059947, 37.937118, 37.937118, 37.937118, 36.805260, 36.805260, 46.267862, 46.267862, 46.267862, 65.399825, 65.399825, 65.399825, 65.399825, 48.982482, 48.982482, 48.982482, 51.001919, 51.001919, 51.001919, 51.001919, 51.001919, 66.000000, 71.873554, 82.353574, 82.353574, 82.353574, 82.353574, 92.725141};

  std::vector<double> x, y, y_noisy, out;
  for (Size i = 0; i < 50; i++)
  {
    x.push_back(speed[i]);
    y.push_back(dist[i]);
  }

  FastLowessSmoothing::lowess(x, y, out);
  for (Size i = 0; i < out.size(); ++i)
  {
    TEST_REAL_SIMILAR(out[i], expected_1[i]);
  }

  out.clear();
  double delta = 0.01 * (x[ x.size()-1 ] - x[0]); // x is sorted
  FastLowessSmoothing::lowess(x, y, 0.2,  3, delta, out);
  for (Size i = 0; i < out.size(); ++i)
  {
    TEST_REAL_SIMILAR(out[i], expected_2[i]);
  }

  // numerical identity with the C++ function
  TOLERANCE_RELATIVE(4e-7);
  TOLERANCE_ABSOLUTE(1e-4);
  double expected_3[] = {4.96545927718688, 4.96545927718688, 13.1244950396665, 13.1244950396665, 15.8586333820983, 18.5796905142177, 21.2803125785285, 21.2803125785285, 21.2803125785285, 24.1292771489265, 24.1292771489265, 27.1195485506035, 27.1195485506035, 27.1195485506035, 27.1195485506035, 30.027276331154, 30.027276331154, 30.027276331154, 30.027276331154, 32.9625061361576, 32.9625061361576, 32.9625061361576, 32.9625061361576, 36.7577283416497, 36.7577283416497, 36.7577283416497, 40.4350745619887, 40.4350745619887, 43.4634917818176, 43.4634917818176, 43.4634917818176, 46.885478946024, 46.885478946024, 46.885478946024, 46.885478946024, 50.7931517254206, 50.7931517254206, 50.7931517254206, 56.4912240928772, 56.4912240928772, 56.4912240928772, 56.4912240928772, 56.4912240928772, 67.5858242314312, 73.0796952693701, 78.6431635544, 78.6431635544, 78.6431635544, 78.6431635544, 84.3286980968344};
  double expected_4[] = {6.03040788454055, 6.03040788454055, 12.6788932684282, 12.6788932684282, 15.3837960614806, 18.6688467170581, 22.2275706232724, 22.2275706232724, 22.2275706232724, 23.3064828196959, 23.3064828196959, 21.52537248518, 21.52537248518, 21.52537248518, 21.52537248518, 34.8827348652577, 34.8827348652577, 34.8827348652577, 34.8827348652577, 47.0599472320042, 47.0599472320042, 47.0599472320042, 47.0599472320042, 37.9371179560115, 37.9371179560115, 37.9371179560115, 36.8052597644327, 36.8052597644327, 46.2678618410954, 46.2678618410954, 46.2678618410954, 65.3998245907766, 65.3998245907766, 65.3998245907766, 65.3998245907766, 48.9824817807382, 48.9824817807382, 48.9824817807382, 51.0019185064708, 51.0019185064708, 51.0019185064708, 51.0019185064708, 51.0019185064708, 65.9999999999999, 71.8735541744287, 82.3535742388261, 82.3535742388261, 82.3535742388261, 82.3535742388261, 92.7251407107177};

  out.clear();
  FastLowessSmoothing::lowess(x, y, out);
  for (Size i = 0; i < out.size(); ++i)
  {
    TEST_REAL_SIMILAR(out[i], expected_3[i]);
  }
  out.clear();
  FastLowessSmoothing::lowess(x, y, 0.2,  3, delta, out);
  for (Size i = 0; i < out.size(); ++i)
  {
    TEST_REAL_SIMILAR(out[i], expected_4[i]);
  }
}
END_SECTION

// trying to fit a quadratic function -> wont work so well, obviously
TOLERANCE_RELATIVE(1.06);
START_SECTION([FastLowessSmoothing]void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&))
{
  std::vector<double> x, y, y_noisy, out, expect;

  //exact data -> sample many points
  for (Size i = 1; i <= 10000; ++i)
  {
    x.push_back(i/500.0);
    y.push_back(targetFunction(i/500.0));
    expect.push_back(targetFunction(i/500.0));
  }

  //noisy data
  // make some noise
  boost::random::mt19937 rnd_gen_;
  for (Size i = 0; i < y.size(); ++i)
  {
    boost::normal_distribution<float> udist (y.at(i), 0.05);
    y_noisy.push_back( udist(rnd_gen_) );
  }

  FastLowessSmoothing::lowess(x, y, 0.02, 3, 0.2, out);
  for (Size i = 0; i < out.size(); ++i)
  {
    TEST_REAL_SIMILAR(out[i], expect[i]);
  }

  out.clear();

  FastLowessSmoothing::lowess(x, y_noisy, 0.02, 3, 0.2, out);
  for (Size i = 0; i < out.size(); ++i)
  {
    TEST_REAL_SIMILAR(out[i], expect[i]);
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


