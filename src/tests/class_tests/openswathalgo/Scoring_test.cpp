// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include "OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h"

#include "OpenMS/OPENSWATHALGO/ALGO/Scoring.h"

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

#include <algorithm>
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
/*
  # example python code of two reference implementations
  # see https://stackoverflow.com/questions/2827393/angles-between-two-n-dimensional-vectors-in-python/13849249#13849249
  
  import numpy as np

  def unit_vector(vector):
      """ Returns the unit vector of the vector.  """
      return np.array(vector) / max(1e-15, np.linalg.norm(vector))

  def angle_between(v1, v2):
      """ Returns the angle in radians between vectors 'v1' and 'v2'::

              >>> angle_between((1, 0, 0), (0, 1, 0))
              1.5707963267948966
              >>> angle_between((1, 0, 0), (1, 0, 0))
              0.0
              >>> angle_between((1, 0, 0), (-1, 0, 0))
              3.141592653589793
              >>> angle_between((0, 0, 0), (0, 0, 0))  # error or pi/2?
              1.5707963267948966
      """
      v1_u = unit_vector(v1)
      v2_u = unit_vector(v2)
      return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

  def spectral_angle(v1, v2):
      """ Returns the angle in radians between vectors 'v1' and 'v2'::

              >>> spectral_angle((1, 0, 0), (0, 1, 0))
              1.5707963267948966
              >>> spectral_angle((1, 0, 0), (1, 0, 0))
              0.0
              >>> spectral_angle((1, 0, 0), (-1, 0, 0))
              3.141592653589793
              >>> spectral_angle((0, 0, 0), (0, 0, 0))  # error or pi/2?
              1.5707963267948966
      """
      numer = np.dot(v1, v2)
      v1_u = np.sqrt(np.dot(v1, v1))
      v2_u = np.sqrt(np.dot(v2, v2))
      denom = v1_u * v2_u
      theta = 0.0 if denom == 0 else numer / denom
      return np.arccos(np.clip(theta, -1.0, 1.0))

  vecs = [
      ((1, 0, 0), (0, 1, 0)),
      ((1, 0, 0), (1, 0, 0)),
      ((1, 0, 0), (-1, 0, 0)),
      ((0, 0, 0), (0, 0, 0)),
  ]
  for i in range(10):
      vecs.append((np.random.uniform(size=3), np.random.uniform(size=3)))
  for v1, v2 in vecs:
      a = angle_between(v1, v2)
      b = spectral_angle(v1, v2)
      if a != b:
          print(f'Failed:\n\tv1 = {v1}\n\tv2 = {v2}\n\ta = {a}\n\tb = {b}\n\ta - b = {a - b}')
*/
  
  static constexpr double pi{3.141592653589793};
  static constexpr double piOver2{0.5 * pi};

  auto spectralAngle = [](
    std::vector<double> d1,
    std::vector<double> d2
  ) -> double {
    return Scoring::SpectralAngle(&d1[0], &d2[0], d1.size());
  };

  // previous unit test
  TEST_REAL_SIMILAR (
    spectralAngle({0,1,3,5,2,0}, {1,3,5,2,0,0}),
    0.7699453419277419
  )

  // zero
  TEST_REAL_SIMILAR (
    spectralAngle({0, 0, 0}, {0, 0, 0}),
    piOver2
  )

  // same
  TEST_REAL_SIMILAR (
    spectralAngle({1, 0, 0}, {1, 0, 0}),
    0.0
  )

  // reversed
  TEST_REAL_SIMILAR (
    spectralAngle({1, 0, 0}, {-1, 0, 0}),
    pi
  )

  // orthogonal
  TEST_REAL_SIMILAR (
    spectralAngle({1, 0, 0}, {0, 1, 0}),
    piOver2
  )

  // random from python
  TEST_REAL_SIMILAR (
    spectralAngle({0.03174064, 0.11582065, 0.63258941}, {0.71882213, 0.00087569, 0.36516896}),
    1.0597217204768459
  )
  TEST_REAL_SIMILAR (
    spectralAngle({0.6608937, 0.0726909, 0.40912141}, {0.52081914, 0.71088, 0.0175557}),
    0.9449782659258582
  )
  TEST_REAL_SIMILAR (
    spectralAngle({0.58858475, 0.08963515, 0.08578046}, {0.76180969, 0.72763536, 0.50090751}),
    0.6547156284689354
  )
  TEST_REAL_SIMILAR (
    spectralAngle({0.08653022, 0.11595108, 0.74268632}, {0.55176333, 0.16783033, 0.70364679}),
    0.5418305329889055
  )
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

BOOST_AUTO_TEST_CASE(test_computeAndAppendRank)
{
/*
* Requires Octave with installed MIToolbox

y = [5.97543668746948 4.2749171257019 3.3301842212677 4.08597040176392 5.50307035446167 5.24326848983765 8.40812492370605 2.83419919013977 6.94378805160522 7.69957494735718 4.08597040176392]';

[~, ~, y_ranking] = unique(y);
% Note: Matlab handles ties differently than Scoring::computeAndAppendRank, but this makes no difference for MI estimation.

*/



  /*static const double arr2[] =
  {
    15.8951349258423, 41.5446395874023, 76.0746307373047, 109.069435119629, 111.90364074707, 169.79216003418,
    121.043930053711, 63.0136985778809, 44.6150207519531, 21.4926776885986, 7.93575811386108
  };*/
  std::vector<double> data1 {5.97543668746948, 4.2749171257019, 3.3301842212677, 4.08597040176392, 5.50307035446167, 5.24326848983765,
                  8.40812492370605, 2.83419919013977, 6.94378805160522, 7.69957494735718, 4.08597040176392};

  std::vector<double> data2 {15.8951349258423, 41.5446395874023, 76.0746307373047, 109.069435119629, 111.90364074707, 169.79216003418,
                             121.043930053711, 63.0136985778809, 44.6150207519531, 21.4926776885986, 7.93575811386108};

  std::vector<unsigned int> result;
  Scoring::computeAndAppendRank(data1, result);

  TEST_EQUAL (result[0],7);
  TEST_EQUAL (result[1],4);
  TEST_EQUAL (result[2],1);
  TEST_EQUAL (result[3],2);
  TEST_EQUAL (result[4],6);
  TEST_EQUAL (result[5],5);
  TEST_EQUAL (result[6],10);
  TEST_EQUAL (result[7],0);
  TEST_EQUAL (result[8],8);
  TEST_EQUAL (result[9],9);
  TEST_EQUAL (result[10],2);

}
END_SECTION

BOOST_AUTO_TEST_CASE(test_rankedMutualInformation)
{
/*
* Requires Octave with installed MIToolbox

y = [5.97543668746948 4.2749171257019 3.3301842212677 4.08597040176392 5.50307035446167 5.24326848983765 8.40812492370605 2.83419919013977 6.94378805160522 7.69957494735718 4.08597040176392]';
x = [15.8951349258423 41.5446395874023 76.0746307373047 109.069435119629 111.90364074707 169.79216003418 121.043930053711 63.0136985778809 44.6150207519531 21.4926776885986 7.93575811386108]';

[~, ~, y_ranking] = unique(y);
[~, ~, x_ranking] = unique(x);

m1 = mi(x_ranking,y_ranking)
*/

  std::vector<double> data1 = {5.97543668746948, 4.2749171257019, 3.3301842212677, 4.08597040176392, 5.50307035446167, 5.24326848983765,
                                 8.40812492370605, 2.83419919013977, 6.94378805160522, 7.69957494735718, 4.08597040176392};
  std::vector<double> data2 = {15.8951349258423, 41.5446395874023, 76.0746307373047, 109.069435119629, 111.90364074707, 169.79216003418,
                               121.043930053711, 63.0136985778809, 44.6150207519531, 21.4926776885986, 7.93575811386108};
  std::vector<unsigned int> rank_vec1, rank_vec2;
  unsigned int max_rank1 = Scoring::computeAndAppendRank(data1, rank_vec1);
  unsigned int max_rank2 = Scoring::computeAndAppendRank(data2, rank_vec2);

  unsigned int max_rank_check1 = *std::max_element(rank_vec1.begin(), rank_vec1.end());
  unsigned int max_rank_check2 = *std::max_element(rank_vec2.begin(), rank_vec2.end());

  TEST_EQUAL (max_rank1, max_rank_check1);
  TEST_EQUAL (max_rank2, max_rank_check2);

  double result = Scoring::rankedMutualInformation(rank_vec1, rank_vec2, max_rank1, max_rank2);

  TEST_REAL_SIMILAR (result, 3.2776);

  rank_vec1 = {0};
  rank_vec2 = {0};
  result = Scoring::rankedMutualInformation(rank_vec1, rank_vec2, 0, 0);
  TEST_REAL_SIMILAR (result, 0);

  rank_vec1 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  rank_vec2 = {0, 1, 5, 4, 4, 2, 3, 1, 0, 2};
  result = Scoring::rankedMutualInformation(rank_vec1, rank_vec2, 0, 5);
  TEST_REAL_SIMILAR (result, 0);

  double result_symmetric = Scoring::rankedMutualInformation(rank_vec2, rank_vec1, 5, 0);
  TEST_REAL_SIMILAR (result, result_symmetric);

  rank_vec1 = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7};
  rank_vec2 = {0, 1, 5, 4, 4, 2, 3, 1, 0, 2, 6, 7, 7, 6, 5, 3};
  result = Scoring::rankedMutualInformation(rank_vec1, rank_vec2, 7, 7);
  TEST_REAL_SIMILAR (result, 2);

  rank_vec1 = {0, 1, 2, 3, 4, 4, 5, 6, 5, 1};
  rank_vec2 = {6, 7, 8, 4, 5, 1, 2, 0, 3, 0};
  result = Scoring::rankedMutualInformation(rank_vec1, rank_vec2, 6, 8);
  TEST_REAL_SIMILAR (result, 2.52193);

  result_symmetric = Scoring::rankedMutualInformation(rank_vec2, rank_vec1, 8, 6);
  TEST_REAL_SIMILAR (result, result_symmetric);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
