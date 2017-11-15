// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/OpenSwathAlgoConfig.h"

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/Scoring.h"

#ifdef USE_BOOST_UNIT_TEST

// include boost unit test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE MyTest
#include <boost/test/unit_test.hpp>
// macros for boost
#define EPS_05 boost::test_tools::fraction_tolerance(1.e-5)
#define TEST_REAL_SIMILAR(val1, val2) \
  BOOST_CHECK ( boost::test_tools::check_is_close(val1, val2, EPS_05 ));
#define TEST_EQUAL(val1, val2) BOOST_CHECK_EQUAL(val1, val2);
#define END_SECTION
#define START_TEST(var1, var2)
#define END_TEST

#else

#include <OpenMS/CONCEPT/ClassTest.h>
#define BOOST_AUTO_TEST_CASE START_SECTION
using namespace OpenMS;

#endif

using namespace std;
using namespace OpenSwath;

///////////////////////////

START_TEST(Scoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_CASE(double_NormalizedManhattanDist_test)
{
  // Numpy 
  // arr1 = [ 0,1,3,5,2,0 ];
  // arr2 = [ 1,3,5,2,0,0 ];
  // arr1 = (arr1 / (sum(arr1) *1.0) ) 
  // arr2 = (arr2 / (sum(arr2) *1.0) ) 
  // deltas = [ abs(a-b) for (a,b) in zip(arr1, arr2) ]
  // sum(deltas) / 6

  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  TEST_REAL_SIMILAR (Scoring::NormalizedManhattanDist(&data1[0], &data2[0], 6), 0.15151515)
}
END_SECTION

BOOST_AUTO_TEST_CASE(double_RootMeanSquareDeviation_test)
{
  // Numpy 
  // arr1 = [ 0,1,3,5,2,0 ];
  // arr2 = [ 1,3,5,2,0,0 ];
  // res = [ (a-b)*(a-b) for (a,b) in zip(arr1, arr2) ]
  // sqrt(sum(res)/6.0)


  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  TEST_REAL_SIMILAR (Scoring::RootMeanSquareDeviation(&data1[0], &data2[0], 6), 1.91485421551)
}
END_SECTION

BOOST_AUTO_TEST_CASE(double_SpectralAngle_test)
{
  // import math 
  // arr1 = [ 0,1,3,5,2,0 ];
  // arr2 = [ 1,3,5,2,0,0 ];
  // dotprod = sum([ (a*b) for (a,b) in zip(arr1, arr2) ])
  // lenx = sqrt(sum([ (a*a) for (a,b) in zip(arr1, arr2) ]))
  // leny = sqrt(sum([ (b*b) for (a,b) in zip(arr1, arr2) ]))
  // math.acos(dotprod/(lenx*leny))


  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  TEST_REAL_SIMILAR (Scoring::SpectralAngle(&data1[0], &data2[0], 6), 0.7699453419277419)

      /*
      normalize_sum(x, n);
      normalize_sum(y, n);
      */
}
END_SECTION

BOOST_AUTO_TEST_CASE(void_normalize_sum_test)
// void normalize_sum(double x[], unsigned int n)
{
  // arr1 = [ 0,1,3,5,2,0 ];
  // n_arr1 = (arr1 / (sum(arr1) *1.0) )
  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  Scoring::normalize_sum(&data1[0], 6);
  TEST_REAL_SIMILAR (data1[0], 0.0)
  TEST_REAL_SIMILAR (data1[1], 0.09090909)
  TEST_REAL_SIMILAR (data1[2], 0.27272727)
  TEST_REAL_SIMILAR (data1[3], 0.45454545)
  TEST_REAL_SIMILAR (data1[4], 0.18181818)
  TEST_REAL_SIMILAR (data1[5], 0.0)
}
END_SECTION

BOOST_AUTO_TEST_CASE(standardize_data_test)
//START_SECTION((void MRMFeatureScoring::standardize_data(std::vector<double>& data)))
{
  // Numpy 
  // arr1 = [ 0,1,3,5,2,0 ];
  // arr2 = [ 1,3,5,2,0,0 ];
  // (arr1 - mean(arr1) ) / std(arr1)
  // (arr2 - mean(arr2) ) / std(arr2)
  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  Scoring::standardize_data(data1);
  Scoring::standardize_data(data2);

  TEST_REAL_SIMILAR (data1[0], -1.03479296);
  TEST_REAL_SIMILAR (data1[1], -0.47036043);
  TEST_REAL_SIMILAR (data1[2], 0.65850461);
  TEST_REAL_SIMILAR (data1[3], 1.78736965);
  TEST_REAL_SIMILAR (data1[4], 0.09407209);
  TEST_REAL_SIMILAR (data1[5], -1.03479296);

  TEST_REAL_SIMILAR (data2[0], -0.47036043);
  TEST_REAL_SIMILAR (data2[1], 0.65850461);
  TEST_REAL_SIMILAR (data2[2], 1.78736965);
  TEST_REAL_SIMILAR (data2[3], 0.09407209);
  TEST_REAL_SIMILAR (data2[4], -1.03479296);
  TEST_REAL_SIMILAR (data2[5], -1.03479296);

}
END_SECTION

BOOST_AUTO_TEST_CASE(test_calculateCrossCorrelation)
//START_SECTION((MRMFeatureScoring::XCorrArrayType MRMFeatureScoring::calculateCrossCorrelation(std::vector<double>& data1, std::vector<double>& data2, int maxdelay, int lag)))
{

  // Numpy 
  // arr1 = [ 0,1,3,5,2,0 ];
  // arr2 = [ 1,3,5,2,0,0 ];
  // data1 = (arr1 - mean(arr1) ) / std(arr1)
  // data2 = (arr2 - mean(arr2) ) / std(arr2)
  // correlate(data1, data2, "same") / 6.0

  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  Scoring::standardize_data(data1);
  Scoring::standardize_data(data2);

  OpenSwath::Scoring::XCorrArrayType result = Scoring::calculateCrossCorrelation(data1, data2, 2, 1);
  for(OpenSwath::Scoring::XCorrArrayType::iterator it = result.begin(); it != result.end(); it++)
  {
    it->second = it->second / 6.0;
  }

  TEST_REAL_SIMILAR (result.data[4].second, -0.7374631);    // find( 2)->
  TEST_REAL_SIMILAR (result.data[3].second, -0.567846);     // find( 1)->
  TEST_REAL_SIMILAR (result.data[2].second,  0.4159292);    // find( 0)->
  TEST_REAL_SIMILAR (result.data[1].second,  0.8215339);    // find(-1)->
  TEST_REAL_SIMILAR (result.data[0].second,  0.15634218);   // find(-2)->

  TEST_EQUAL (result.data[4].first, 2)
  TEST_EQUAL (result.data[3].first, 1)
  TEST_EQUAL (result.data[2].first, 0)
  TEST_EQUAL (result.data[1].first, -1)
  TEST_EQUAL (result.data[0].first, -2)
}
END_SECTION

BOOST_AUTO_TEST_CASE(test_MRMFeatureScoring_normalizedCrossCorrelation)
//START_SECTION((MRMFeatureScoring::XCorrArrayType MRMFeatureScoring::normalizedCrossCorrelation(std::vector<double>& data1, std::vector<double>& data2, int maxdelay, int lag)))
{

  // Numpy 
  // arr1 = [ 0,1,3,5,2,0 ];
  // arr2 = [ 1,3,5,2,0,0 ];
  // data1 = (arr1 - mean(arr1) ) / std(arr1)
  // data2 = (arr2 - mean(arr2) ) / std(arr2)
  // correlate(data1, data2, "same")

  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  OpenSwath::Scoring::XCorrArrayType result = Scoring::normalizedCrossCorrelation(data1, data2, 2, 1);

  TEST_REAL_SIMILAR (result.data[4].second, -0.7374631);  // .find( 2)
  TEST_REAL_SIMILAR (result.data[3].second, -0.567846);   // .find( 1)
  TEST_REAL_SIMILAR (result.data[2].second,  0.4159292);  // .find( 0)
  TEST_REAL_SIMILAR (result.data[1].second,  0.8215339);  // .find(-1)
  TEST_REAL_SIMILAR (result.data[0].second,  0.15634218); // .find(-2)
    
  TEST_EQUAL (result.data[4].first, 2)
  TEST_EQUAL (result.data[3].first, 1)
  TEST_EQUAL (result.data[2].first, 0)
  TEST_EQUAL (result.data[1].first, -1)
  TEST_EQUAL (result.data[0].first, -2)
}
END_SECTION

BOOST_AUTO_TEST_CASE(test_MRMFeatureScoring_calcxcorr_legacy_mquest_)
//START_SECTION((MRMFeatureScoring::XCorrArrayType MRMFeatureScoring::calcxcorr(std::vector<double>& data1, std::vector<double>& data2, bool normalize)))
{

  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  OpenSwath::Scoring::XCorrArrayType result = Scoring::calcxcorr_legacy_mquest_(data1, data2, true);
  TEST_EQUAL (result.data.size(), 13)

  TEST_REAL_SIMILAR (result.data[4+4].second, -0.7374631);    // .find( 2)
  TEST_REAL_SIMILAR (result.data[3+4].second, -0.567846);     // .find( 1)
  TEST_REAL_SIMILAR (result.data[2+4].second,  0.4159292);    // .find( 0)
  TEST_REAL_SIMILAR (result.data[1+4].second,  0.8215339);    // .find(-1)
  TEST_REAL_SIMILAR (result.data[0+4].second,  0.15634218);   // .find(-2)
    
  TEST_EQUAL (result.data[4+4].first, 2)
  TEST_EQUAL (result.data[3+4].first, 1)
  TEST_EQUAL (result.data[2+4].first, 0)
  TEST_EQUAL (result.data[1+4].first, -1)
  TEST_EQUAL (result.data[0+4].first, -2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
