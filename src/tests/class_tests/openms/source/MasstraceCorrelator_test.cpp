// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MasstraceCorrelator.h>
///////////////////////////

#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>

using namespace OpenMS;
using namespace std;

class MasstraceCorrelator_facade : MasstraceCorrelator
{
  public:

    void matchMassTraces_(const MasstracePointsType& hull_points1, const MasstracePointsType& hull_points2,
        std::vector<double>& vec1, std::vector<double>& vec2, double mindiff,
        double padEnds = true)
    {
      MasstraceCorrelator::matchMassTraces_(hull_points1, hull_points2, vec1, vec2, mindiff, padEnds);
    }

};
START_TEST(CorrelateMasstraces, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((virtual void matchMassTraces_()))
{

  MasstraceCorrelator_facade mtcorr;
  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  static const double arr3[] = {0,1,2,3,4,5};
  static const double arr4[] = {0,1,2,4,5,6};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  std::vector<double> rt1 (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]) );
  std::vector<double> rt2 (arr4, arr4 + sizeof(arr4) / sizeof(arr4[0]) );

  std::vector<double> vec1;
  std::vector<double> vec2;

  std::vector<std::pair<double, double> > data1_2d;
  std::vector<std::pair<double, double> > data2_2d;

  data1_2d.clear();
  data2_2d.clear();
  for(Size i = 0; i < data1.size(); i++)
  {
    data1_2d.push_back(std::make_pair(rt1[i], data1[i]) );
    data2_2d.push_back(std::make_pair(rt2[i], data2[i]) );
  }

  vec1.clear(); vec2.clear();
  mtcorr.matchMassTraces_(data1_2d, data2_2d, vec1, vec2, 0.1);

  TEST_EQUAL(vec1.size(), 7);
  TEST_EQUAL(vec2.size(), 7);

  TEST_EQUAL(vec1[0], 0);
  TEST_EQUAL(vec1[1], 1);
  TEST_EQUAL(vec1[2], 3);
  TEST_EQUAL(vec1[3], 5);
  TEST_EQUAL(vec1[4], 2);
  TEST_EQUAL(vec1[5], 0);
  TEST_EQUAL(vec1[6], 0);

  TEST_EQUAL(vec2[0], 1);
  TEST_EQUAL(vec2[1], 3);
  TEST_EQUAL(vec2[2], 5);
  TEST_EQUAL(vec2[3], 0);
  TEST_EQUAL(vec2[4], 2);
  TEST_EQUAL(vec2[5], 0);
  TEST_EQUAL(vec2[6], 0);

  vec1.clear(); vec2.clear();
  mtcorr.matchMassTraces_(data2_2d, data1_2d, vec2, vec1, 0.1);

  TEST_EQUAL(vec1.size(), 7);
  TEST_EQUAL(vec2.size(), 7);

  TEST_EQUAL(vec1[0], 0);
  TEST_EQUAL(vec1[1], 1);
  TEST_EQUAL(vec1[2], 3);
  TEST_EQUAL(vec1[3], 5);
  TEST_EQUAL(vec1[4], 2);
  TEST_EQUAL(vec1[5], 0);
  TEST_EQUAL(vec1[6], 0);

  TEST_EQUAL(vec2[0], 1);
  TEST_EQUAL(vec2[1], 3);
  TEST_EQUAL(vec2[2], 5);
  TEST_EQUAL(vec2[3], 0);
  TEST_EQUAL(vec2[4], 2);
  TEST_EQUAL(vec2[5], 0);
  TEST_EQUAL(vec2[6], 0);

  vec1.clear(); vec2.clear();
  mtcorr.matchMassTraces_(data1_2d, data2_2d, vec1, vec2, 1.5);

  TEST_EQUAL(vec1.size(), 6);
  TEST_EQUAL(vec2.size(), 6);

  TEST_EQUAL(vec1[0], 0);
  TEST_EQUAL(vec1[1], 1);
  TEST_EQUAL(vec1[2], 3);
  TEST_EQUAL(vec1[3], 5);
  TEST_EQUAL(vec1[4], 2);
  TEST_EQUAL(vec1[5], 0);

  TEST_EQUAL(vec2[0], 1);
  TEST_EQUAL(vec2[1], 3);
  TEST_EQUAL(vec2[2], 5);
  TEST_EQUAL(vec2[3], 2);
  TEST_EQUAL(vec2[4], 0);
  TEST_EQUAL(vec2[5], 0);

}
END_SECTION

START_SECTION((virtual void match_elution_arrays_no_padding()))
{

  MasstraceCorrelator_facade mtcorr;
  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  static const double arr3[] = {0,1,2,3,4,5};
  static const double arr4[] = {-1,1,2,4,5,6};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  std::vector<double> rt1 (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]) );
  std::vector<double> rt2 (arr4, arr4 + sizeof(arr4) / sizeof(arr4[0]) );

  std::vector<double> vec1;
  std::vector<double> vec2;

  std::vector<std::pair<double, double> > data1_2d;
  std::vector<std::pair<double, double> > data2_2d;

  data1_2d.clear();
  data2_2d.clear();
  for(Size i = 0; i < data1.size(); i++)
  {
    data1_2d.push_back(std::make_pair(rt1[i], data1[i]) );
    data2_2d.push_back(std::make_pair(rt2[i], data2[i]) );
  }

  vec1.clear(); vec2.clear();
  TEST_EQUAL(vec1.size(), 0);
  TEST_EQUAL(vec2.size(), 0);
  // if we do not pad the ends, this means that we do not add zeros to the first vector that is shorter in RT
  bool pad_ends = false;
  mtcorr.matchMassTraces_(data1_2d, data2_2d, vec1, vec2, 0.1, pad_ends);

  TEST_EQUAL(vec1.size(), 5);
  TEST_EQUAL(vec2.size(), 5);

  TEST_EQUAL(vec1[0], 1);
  TEST_EQUAL(vec1[1], 3);
  TEST_EQUAL(vec1[2], 5);
  TEST_EQUAL(vec1[3], 2);
  TEST_EQUAL(vec1[4], 0);

  TEST_EQUAL(vec2[0], 3);
  TEST_EQUAL(vec2[1], 5);
  TEST_EQUAL(vec2[2], 0);
  TEST_EQUAL(vec2[3], 2);
  TEST_EQUAL(vec2[4], 0);

  vec1.clear(); vec2.clear();
  TEST_EQUAL(vec1.size(), 0);
  TEST_EQUAL(vec2.size(), 0);
  // if we do pad the ends, this means that we do add zeros to the first vector that is shorter in RT
  mtcorr.matchMassTraces_(data1_2d, data2_2d, vec1, vec2, 0.1, true);

  TEST_EQUAL(vec1.size(), 8);
  TEST_EQUAL(vec2.size(), 8);

  TEST_EQUAL(vec1[0], 0); // -1 
  TEST_EQUAL(vec1[1], 0); // 0  
  TEST_EQUAL(vec1[2], 1); // 1  
  TEST_EQUAL(vec1[3], 3); // 2  
  TEST_EQUAL(vec1[4], 5); // 3  
  TEST_EQUAL(vec1[5], 2); // 4  
  TEST_EQUAL(vec1[6], 0); // 5  
  TEST_EQUAL(vec1[7], 0); // 6  

  TEST_EQUAL(vec2[0], 1); // -1
  TEST_EQUAL(vec2[1], 0); // 0
  TEST_EQUAL(vec2[2], 3); // 1
  TEST_EQUAL(vec2[3], 5); // 2
  TEST_EQUAL(vec2[4], 0); // 3
  TEST_EQUAL(vec2[5], 2); // 4
  TEST_EQUAL(vec2[6], 0); // 5
  TEST_EQUAL(vec2[7], 0); // 6

}
END_SECTION

START_SECTION((virtual void scoreHullpoints()))
{

  OpenMS::MasstraceCorrelator mtcorr;
  static const double arr1[] = {0,1,3,5,2,0};
  static const double arr2[] = {1,3,5,2,0,0};
  static const double arr3[] = {0,1,2,3,4,5};
  static const double arr4[] = {0,1,2,4,5,6};
  std::vector<double> data1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  std::vector<double> data2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  std::vector<double> rt1 (arr3, arr3 + sizeof(arr3) / sizeof(arr3[0]) );
  std::vector<double> rt2 (arr4, arr4 + sizeof(arr4) / sizeof(arr4[0]) );

  std::vector<std::pair<double, double> > data1_2d;
  std::vector<std::pair<double, double> > data2_2d;

  OpenSwath::Scoring::standardize_data(data1);
  OpenSwath::Scoring::standardize_data(data2);

  for(Size i = 0; i < data1.size(); i++)
  {
    data1_2d.push_back( std::make_pair(rt1[i], data1[i]) );
    data2_2d.push_back( std::make_pair(rt1[i], data2[i]) );
  }

  OpenSwath::Scoring::XCorrArrayType result = OpenSwath::Scoring::calculateCrossCorrelation(data1, data2, 2, 1);
  for(OpenSwath::Scoring::XCorrArrayType::iterator it = result.begin(); it != result.end(); ++it)
  {
    it->second = it->second / 6.0;
  }

  TEST_EQUAL (result.data[0].first, -2);
  TEST_EQUAL (result.data[1].first, -1);
  TEST_EQUAL (result.data[2].first, 0);
  TEST_EQUAL (result.data[3].first, 1);
  TEST_EQUAL (result.data[4].first, 2);

  TEST_REAL_SIMILAR (result.data[4].second, -0.7374631);   // .find( 2)->
  TEST_REAL_SIMILAR (result.data[3].second, -0.567846);    // .find( 1)->
  TEST_REAL_SIMILAR (result.data[2].second,  0.4159292);   // .find( 0)->
  TEST_REAL_SIMILAR (result.data[1].second,  0.8215339);   // .find(-1)->
  TEST_REAL_SIMILAR (result.data[0].second,  0.15634218);  // .find(-2)->

  double min_pearson_score = -1.1; int maxlag = data1_2d.size();
  int lag; double lag_intensity; double pearson_score;
  mtcorr.scoreHullpoints(data1_2d, data2_2d, lag, lag_intensity, pearson_score, min_pearson_score, maxlag);
  TEST_EQUAL(lag, -1);
  TEST_REAL_SIMILAR(lag_intensity, 0.821534);
  TEST_REAL_SIMILAR(pearson_score, 0.41593);

  // now we use different RT data for the 2nd data array
  data1_2d.clear();
  data2_2d.clear();
  for(Size i = 0; i < data1.size(); i++)
  {
    data1_2d.push_back( std::make_pair(rt1[i], data1[i]) );
    data2_2d.push_back( std::make_pair(rt2[i], data2[i]) );
  }

  // if we allow for an RT difference of more than 1, we should get the same result as above
  mtcorr.scoreHullpoints(data1_2d, data2_2d, lag, lag_intensity, pearson_score, min_pearson_score, maxlag, 1.5);
  TEST_EQUAL(lag, -1);
  TEST_REAL_SIMILAR(lag_intensity, 0.821534);
  TEST_REAL_SIMILAR(pearson_score, 0.41593);

  // if the allowed difference in RT is less than 1, the algorithm will substitute zeros
  mtcorr.scoreHullpoints(data1_2d, data2_2d, lag, lag_intensity, pearson_score, min_pearson_score, maxlag, 0.1);
  TEST_EQUAL(lag, -1);
  TEST_REAL_SIMILAR(lag_intensity, 0.625368);
  TEST_REAL_SIMILAR(pearson_score, 0.405604);

}
END_SECTION

START_SECTION((virtual void createPseudoSpectra()))
{

  ConsensusMap masstraces;
  ConsensusXMLFile().load(OPENMS_GET_TEST_DATA_PATH("Masstraces_Testdata.consensusXML"), masstraces);

  MSExperiment pseudo_spectra;
  masstraces.sortByIntensity(true);
  OpenMS::MasstraceCorrelator mtcorr;
  mtcorr.createPseudoSpectra(masstraces, pseudo_spectra, 0, 0.7, 1, 3);

  TEST_EQUAL(pseudo_spectra.size(), 3);
  TEST_EQUAL(pseudo_spectra[1].size(), 1);
  TEST_REAL_SIMILAR(pseudo_spectra[1].getRT(), 4203.0);
  TEST_REAL_SIMILAR(pseudo_spectra[1][0].getMZ(), 668.5);

  pseudo_spectra.clear(true);
  mtcorr.createPseudoSpectra(masstraces, pseudo_spectra, 1, 0.7, 1, 3);
  TEST_EQUAL(pseudo_spectra.size(), 2);

  TEST_EQUAL(pseudo_spectra[0].size(), 2);
  TEST_EQUAL(pseudo_spectra[1].size(), 2);

  TEST_REAL_SIMILAR(pseudo_spectra[0].getRT(), 5201.0);
  TEST_REAL_SIMILAR(pseudo_spectra[0][0].getMZ(), 568.5);

  TEST_REAL_SIMILAR(pseudo_spectra[1].getRT(), 5203.0);
  TEST_REAL_SIMILAR(pseudo_spectra[1][0].getMZ(), 768.5);

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



