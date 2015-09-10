// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FILTERING/SMOOTHING/FastLowessSmoothing.h>
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

TOLERANCE_RELATIVE(1.025);
START_SECTION([FastLowessSmoothing_cars]void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&))
{

  /*
  In R

     require(graphics)
     
     plot(cars, main = "lowess(cars)")
     lines(lowess(cars), col = 2)
     lines(lowess(cars, f = .2), col = 3)
     legend(5, 120, c(paste("f = ", c("2/3", ".2"))), lty = 1, col = 2:3)
     
  */

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
}
END_SECTION

TOLERANCE_RELATIVE(1.06);
START_SECTION([FastLowessSmoothing]void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&))
{
#if 1
  std::vector<double> x, y, y_noisy, out, expect;

  //exact data
  // for (double i = 1.0; i <= 20.0; i += 1.0)
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

  Size idx = 1;
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

#endif
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



