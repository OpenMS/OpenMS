// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include "OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h"

#include "OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/MockObjects.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h"

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

void fill_mock_objects(MockMRMFeature * imrmfeature, std::vector<std::string>& native_ids)
{
  native_ids.push_back("group1");
  native_ids.push_back("group2");

  
  std::vector<double> intensity1 {5.97543668746948, 4.2749171257019, 3.3301842212677, 4.08597040176392, 5.50307035446167, 5.24326848983765,
                                  8.40812492370605, 2.83419919013977, 6.94378805160522, 7.69957494735718, 4.08597040176392};

  std::vector<double> intensity2 {15.8951349258423, 41.5446395874023, 76.0746307373047, 109.069435119629, 111.90364074707, 169.79216003418,
                                  121.043930053711, 63.0136985778809, 44.6150207519531, 21.4926776885986, 7.93575811386108};

  std::vector<double> ms1intensity {0.0, 110.0, 200.0, 270.0, 320.0, 350.0, 360.0, 350.0, 320.0, 270.0, 200.0};

  boost::shared_ptr<MockFeature> f1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> f2_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> ms1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  f1_ptr->m_intensity_vec = intensity1;
  f2_ptr->m_intensity_vec = intensity2;
  ms1_ptr->m_intensity_vec = ms1intensity;
  std::map<std::string, boost::shared_ptr<MockFeature> > features;
  features["group1"] = f1_ptr;
  features["group2"] = f2_ptr;
  imrmfeature->m_features = features; // add features

  std::map<std::string, boost::shared_ptr<MockFeature> > ms1_features;
  ms1_features["ms1trace"] = ms1_ptr;
  imrmfeature->m_precursor_features = ms1_features; // add ms1 feature
}

void fill_mock_objects2(MockMRMFeature * imrmfeature, std::vector<std::string>& precursor_ids, std::vector<std::string>& native_ids)
{
  precursor_ids.push_back("ms1trace1");
  precursor_ids.push_back("ms1trace2");
  precursor_ids.push_back("ms1trace3");

  native_ids.push_back("group1");
  native_ids.push_back("group2");


  std::vector<double> intensity1 {5.97543668746948, 4.2749171257019, 3.3301842212677, 4.08597040176392, 5.50307035446167, 5.24326848983765,
                                  8.40812492370605, 2.83419919013977, 6.94378805160522, 7.69957494735718, 4.08597040176392};

  std::vector<double> intensity2 {15.8951349258423, 41.5446395874023, 76.0746307373047, 109.069435119629, 111.90364074707, 169.79216003418,
                                  121.043930053711, 63.0136985778809, 44.6150207519531, 21.4926776885986, 7.93575811386108};

  std::vector<double> ms1intensity1 {0.0, 110.0, 200.0, 270.0, 320.0, 350.0, 360.0, 350.0, 320.0, 270.0, 200.0};

  std::vector<double> ms1intensity2 {10.0, 115.0, 180.0, 280.0, 330.0, 340.0, 390.0, 320.0, 300.0, 250.0, 100.0};

  std::vector<double> ms1intensity3 {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  boost::shared_ptr<MockFeature> f1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> f2_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> ms1_f1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> ms1_f2_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> ms1_f3_ptr = boost::shared_ptr<MockFeature>(new MockFeature());

  f1_ptr->m_intensity_vec = intensity1;
  f2_ptr->m_intensity_vec = intensity2;
  ms1_f1_ptr->m_intensity_vec = ms1intensity1;
  ms1_f2_ptr->m_intensity_vec = ms1intensity2;
  ms1_f3_ptr->m_intensity_vec = ms1intensity3;

  std::map<std::string, boost::shared_ptr<MockFeature> > features;
  features["group1"] = f1_ptr;
  features["group2"] = f2_ptr;
  imrmfeature->m_features = features; // add features

  std::map<std::string, boost::shared_ptr<MockFeature> > ms1_features;
  ms1_features["ms1trace1"] = ms1_f1_ptr;
  ms1_features["ms1trace2"] = ms1_f2_ptr;
  ms1_features["ms1trace3"] = ms1_f3_ptr;
  imrmfeature->m_precursor_features = ms1_features; // add ms1 feature
}

///////////////////////////

START_TEST(MRMScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#ifndef USE_BOOST_UNIT_TEST
    {
      MRMScoring* ptr = nullptr;
      MRMScoring* nullPointer = nullptr;

      START_SECTION(MRMScoring())
          {
            ptr = new MRMScoring();
            TEST_NOT_EQUAL(ptr, nullPointer)
          }
      END_SECTION

      START_SECTION(~MRMScoring())
          {
            delete ptr;
          }
      END_SECTION
    }
#endif


/*
 * Validation of the cross-correlation in Python
 *

from numpy import *

data1 = [5.97543668746948, 4.2749171257019, 3.3301842212677, 4.08597040176392, 5.50307035446167, 5.24326848983765,
       8.40812492370605, 2.83419919013977, 6.94378805160522, 7.69957494735718, 4.08597040176392]
data2 = [15.8951349258423, 41.5446395874023, 76.0746307373047, 109.069435119629, 111.90364074707, 169.79216003418,
       121.043930053711, 63.0136985778809, 44.6150207519531, 21.4926776885986, 7.93575811386108]
ms1data = [0.0, 110.0, 200.0, 270.0, 320.0, 350.0, 360.0, 350.0, 320.0, 270.0, 200.0]
data1 = (data1 - mean(data1) ) / std(data1)
data2 = (data2 - mean(data2) ) / std(data2)
ms1data = (ms1data - mean(ms1data) ) / std(ms1data)
xcorrmatrix_0_0 = correlate(data1, data1, "same") / len(data1)
xcorrmatrix_0_1 = correlate(data1, data2, "same") / len(data1)

max_el0 = max(enumerate(xcorrmatrix_0_0), key= lambda x: x[1])
max_el1 = max(enumerate(xcorrmatrix_0_1), key= lambda x: x[1])

xcorr_deltas = [0, abs(max_el0[0] - max_el1[0]), 0]
xcorr_max = [1, max_el1[1], 1]

mean(xcorr_deltas) + std(xcorr_deltas, ddof=1) # coelution score
# 2.7320508075688772
mean(xcorr_max) # shape score
# 0.13232774079239637

# MS1 level

xcorrvector_1 = correlate(ms1data, data1, "same") / len(data1)
xcorrvector_2 = correlate(ms1data, data2, "same") / len(data2)
max_el0 = max(enumerate(xcorrvector_1), key= lambda x: x[1])
max_el1 = max(enumerate(xcorrvector_2), key= lambda x: x[1])
xcorr_deltas = [0, abs(max_el0[0] - max_el1[0])]
xcorr_max = [max_el0[1], max_el1[1]]

mean(xcorr_deltas) + std(xcorr_deltas, ddof=1) # coelution score
# 1.8213672050459184
mean(xcorr_max) # shape score
# 0.54120799790227003

 *
*/

//START_SECTION((void initializeXCorrMatrix(MRMFeature& mrmfeature, MRMTransitionGroup< SpectrumType, PeakType >& transition_group, bool normalize)))
    BOOST_AUTO_TEST_CASE(initializeXCorrMatrix)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);

          //initialize the XCorr Matrix
          mrmscore.initializeXCorrMatrix(imrmfeature, native_ids);
          delete imrmfeature;

          TEST_EQUAL(mrmscore.getXCorrMatrix().rows(), 2)
          TEST_EQUAL(mrmscore.getXCorrMatrix().cols(), 2)
          TEST_EQUAL(mrmscore.getXCorrMatrix()(0, 0).data.size(), 23)

          // test auto-correlation = xcorrmatrix_0_0
          const OpenSwath::Scoring::XCorrArrayType auto_correlation =
                  mrmscore.getXCorrMatrix()(0, 0);

          TEST_EQUAL(auto_correlation.data[11].first, 0)
          TEST_EQUAL(auto_correlation.data[12].first, 1)
          TEST_EQUAL(auto_correlation.data[10].first, -1)
          TEST_EQUAL(auto_correlation.data[13].first, 2)
          TEST_EQUAL(auto_correlation.data[ 9].first, -2)

          TEST_REAL_SIMILAR(auto_correlation.data[11].second , 1)                   // find(0)->second,
          TEST_REAL_SIMILAR(auto_correlation.data[12].second , -0.227352707759245)  // find(1)->second,
          TEST_REAL_SIMILAR(auto_correlation.data[10].second ,  -0.227352707759245) // find(-1)->second,
          TEST_REAL_SIMILAR(auto_correlation.data[13].second , -0.07501116)         // find(2)->second,
          TEST_REAL_SIMILAR(auto_correlation.data[ 9].second ,  -0.07501116)        // find(-2)->second,

          // test cross-correlation = xcorrmatrix_0_1
          const OpenSwath::Scoring::XCorrArrayType cross_correlation =
                  mrmscore.getXCorrMatrix()(0, 1);

          TEST_REAL_SIMILAR(cross_correlation.data[13].second, -0.31165141)   // find(2)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[12].second, -0.35036919)   // find(1)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[11].second, 0.03129565)    // find(0)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[10].second,  0.30204049)   // find(-1)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[ 9].second,  0.13012441)   // find(-2)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[ 8].second,  0.39698322)   // find(-3)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[ 7].second,  0.16608774)   // find(-4)->second,
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(initializeXCorrPrecursorContrastMatrix)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeXCorrPrecursorContrastMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_EQUAL(mrmscore.getXCorrPrecursorContrastMatrix().rows(), 3)
          TEST_EQUAL(mrmscore.getXCorrPrecursorContrastMatrix().cols(), 2)

          std::vector<double> sum_matrix;

          const auto& cm = mrmscore.getXCorrPrecursorContrastMatrix();
          // Note: the original code depens on col vs. row order and
          // the old code: for (auto e : mrmscore.getXCorrPrecursorContrastMatrix()) fails with different data
          for (size_t r = 0; r != cm.rows(); ++r) 
            for (size_t c = 0; c != cm.cols(); ++c) 
            {
              double sum{0};
              for (size_t i = 0; i < cm(r,c).data.size(); ++i)
              {
                sum += abs(cm(r,c).data[i].second);
              }
              sum_matrix.push_back(sum);
            }
/*
          for (auto e : mrmscore.getXCorrPrecursorContrastMatrix())
          {
          }
*/
          TEST_REAL_SIMILAR(sum_matrix[0], 3.40949220)
          TEST_REAL_SIMILAR(sum_matrix[1], 6.19794611)
          TEST_REAL_SIMILAR(sum_matrix[2], 3.68912454)
          TEST_REAL_SIMILAR(sum_matrix[3], 6.60757921)
          TEST_REAL_SIMILAR(sum_matrix[4], 0)
          TEST_REAL_SIMILAR(sum_matrix[5], 0)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(initializeXCorrPrecursorCombinedMatrix)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeXCorrPrecursorCombinedMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_EQUAL(mrmscore.getXCorrPrecursorCombinedMatrix().rows(), 5)
          TEST_EQUAL(mrmscore.getXCorrPrecursorCombinedMatrix().cols(), 5)

          std::vector<double> sum_matrix;

          const auto& cm = mrmscore.getXCorrPrecursorCombinedMatrix();
          // Note: the original code depens on col vs. row order and
          // the old code: for (auto e : mrmscore.getXCorrPrecursorCombinedMatrix()) fails with different data
          for (size_t r = 0; r != cm.rows(); ++r) 
            for (size_t c = 0; c != cm.cols(); ++c) 
            {
              double sum{0};
              for (size_t i = 0; i < cm(r,c).data.size(); ++i)
              {
                sum += abs(cm(r,c).data[i].second);
              }
              sum_matrix.push_back(sum);
            }

          // Check upper triangular matrix 
          TEST_REAL_SIMILAR(sum_matrix[0], 5.86440677)
          TEST_REAL_SIMILAR(sum_matrix[1], 6.05410398)
          TEST_REAL_SIMILAR(sum_matrix[2], 0)
          TEST_REAL_SIMILAR(sum_matrix[3], 3.40949220)
          TEST_REAL_SIMILAR(sum_matrix[4], 6.19794611)
          TEST_REAL_SIMILAR(sum_matrix[6], 6.30751744)
          TEST_REAL_SIMILAR(sum_matrix[7], 0)
          TEST_REAL_SIMILAR(sum_matrix[8], 3.68912454)
          TEST_REAL_SIMILAR(sum_matrix[9], 6.60757921)
          TEST_REAL_SIMILAR(sum_matrix[12], 0)
          TEST_REAL_SIMILAR(sum_matrix[13], 0)
          TEST_REAL_SIMILAR(sum_matrix[14], 0)
          TEST_REAL_SIMILAR(sum_matrix[18], 3.13711983)
          TEST_REAL_SIMILAR(sum_matrix[19], 3.57832717)
          TEST_REAL_SIMILAR(sum_matrix[24], 6.78303987)
        }
    END_SECTION

/*BOOST_AUTO_TEST_CASE(initializeXCorrPrecursorMatrix)
{
  MockMRMFeature * imrmfeature = new MockMRMFeature();
  MRMScoring mrmscore;

  std::vector<std::string> precursor_ids;
  fill_mock_objects(imrmfeature, precursor_ids);

  //initialize the XCorr vector
  mrmscore.initializeXCorrPrecursorMatrix(imrmfeature, precursor_ids);
  delete imrmfeature;

  TEST_EQUAL(mrmscore.getXCorrPrecursorMatrix().rows(), 3)
  TEST_EQUAL(mrmscore.getXCorrPrecurso().cols(), 2)
}
END_SECTION*/

    BOOST_AUTO_TEST_CASE(initializeXCorrContrastMatrix)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);

          //initialize the XCorr Matrix
          mrmscore.initializeXCorrContrastMatrix(imrmfeature, native_ids, native_ids);
          delete imrmfeature;

          TEST_EQUAL(mrmscore.getXCorrContrastMatrix().rows(), 2)
          TEST_EQUAL(mrmscore.getXCorrContrastMatrix().cols(), 2)
          TEST_EQUAL(mrmscore.getXCorrContrastMatrix()(0, 0).data.size(), 23)

          // test auto-correlation = xcorrmatrix_0_0
          const OpenSwath::Scoring::XCorrArrayType auto_correlation =
                  mrmscore.getXCorrContrastMatrix()(0, 0);
          TEST_REAL_SIMILAR(auto_correlation.data[11].second, 1)                     // find(0)->second,
          TEST_REAL_SIMILAR(auto_correlation.data[12].second, -0.227352707759245)    // find(1)->second,
          TEST_REAL_SIMILAR(auto_correlation.data[10].second,  -0.227352707759245)   // find(-1)->second,
          TEST_REAL_SIMILAR(auto_correlation.data[13].second, -0.07501116)           // find(2)->second,
          TEST_REAL_SIMILAR(auto_correlation.data[ 9].second,  -0.07501116)          // find(-2)->second,

          // // test cross-correlation = xcorrmatrix_0_1
          const OpenSwath::Scoring::XCorrArrayType cross_correlation =
                  mrmscore.getXCorrContrastMatrix()(0, 1);
          TEST_REAL_SIMILAR(cross_correlation.data[13].second, -0.31165141)   // find(2)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[12].second, -0.35036919)   // find(1)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[11].second, 0.03129565)    // find(0)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[10].second,  0.30204049)   // find(-1)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[ 9].second,  0.13012441)   // find(-2)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[ 8].second,  0.39698322)   // find(-3)->second,
          TEST_REAL_SIMILAR(cross_correlation.data[ 7].second,  0.16608774)   // find(-4)->second,
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcXcorrCoelutionScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;
          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);
          mrmscore.initializeXCorrMatrix(imrmfeature, native_ids);
          delete imrmfeature;
          TEST_REAL_SIMILAR(mrmscore.calcXcorrCoelutionScore(), 1 + std::sqrt(3.0)) // mean + std deviation
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcSeparateXcorrContrastCoelutionScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;
          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);
          mrmscore.initializeXCorrContrastMatrix(imrmfeature, native_ids, native_ids);
          delete imrmfeature;
          TEST_REAL_SIMILAR(mrmscore.calcSeparateXcorrContrastCoelutionScore()[0], 1.5)
          TEST_REAL_SIMILAR(mrmscore.calcSeparateXcorrContrastCoelutionScore()[1], 1.5)

        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcXcorrCoelutionWeightedScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;
          static const double weights_[] = { 0.5, 0.5 };
          std::vector<double> weights (weights_, weights_ + sizeof(weights_) / sizeof(weights_[0]) );
          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);
          mrmscore.initializeXCorrMatrix(imrmfeature, native_ids);
          delete imrmfeature;

          // xcorr_deltas = [0, 3, 0] * array([0.25, 2*0.5*0.5,0.25])
          // sum(xcorr_deltas)
          TEST_REAL_SIMILAR(mrmscore.calcXcorrCoelutionWeightedScore(weights), 1.5)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcXcorrShapeScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;
          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);
          mrmscore.initializeXCorrMatrix(imrmfeature, native_ids);
          delete imrmfeature;
          TEST_REAL_SIMILAR(mrmscore.calcXcorrShapeScore(), (1.0 + 0.3969832 + 1.0)/3.0) // mean + std deviation
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcSeparateXcorrContrastShapeScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;
          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);
          mrmscore.initializeXCorrContrastMatrix(imrmfeature, native_ids, native_ids);
          delete imrmfeature;
          TEST_REAL_SIMILAR(mrmscore.calcSeparateXcorrContrastShapeScore()[0], 0.698492)
          TEST_REAL_SIMILAR(mrmscore.calcSeparateXcorrContrastShapeScore()[1], 0.698492)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcXcorrShapeWeightedScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;
          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);
          static const double weights_[] = { 0.5, 0.5 };
          std::vector<double> weights (weights_, weights_ + sizeof(weights_) / sizeof(weights_[0]) );
          mrmscore.initializeXCorrMatrix(imrmfeature, native_ids);
          delete imrmfeature;

          // xcorr_deltas = [1, 0.3969832, 1] * array([0.25, 2*0.5*0.5,0.25])
          // sum(xcorr_deltas)
          TEST_REAL_SIMILAR(mrmscore.calcXcorrShapeWeightedScore(weights), 0.6984916)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(calcXcorrPrecursorContrastCoelutionScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeXCorrPrecursorContrastMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_REAL_SIMILAR(mrmscore.calcXcorrPrecursorContrastCoelutionScore(), 9.5741984 )
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(calcXcorrPrecursorCombinedCoelutionScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeXCorrPrecursorCombinedMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_REAL_SIMILAR(mrmscore.calcXcorrPrecursorCombinedCoelutionScore(), 9.2444789 )
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(calcXcorrPrecursorContrastShapeScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeXCorrPrecursorContrastMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_REAL_SIMILAR(mrmscore.calcXcorrPrecursorContrastShapeScore(), 0.3772868 )
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(calcXcorrPrecursorCombinedShapeScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeXCorrPrecursorCombinedMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_REAL_SIMILAR(mrmscore.calcXcorrPrecursorCombinedShapeScore(), 0.5079334 )
        }
    END_SECTION

//START_SECTION((virtual void test_calcLibraryScore()))
    BOOST_AUTO_TEST_CASE(test_Library_score)
        {
          /*
           * Validation in Python of the different library correlation scores
           *

        from numpy import *
        data1 = array([1,10000,2000])
        data2 = array([782.380737304688, 58.3845062255859, 58.3845062255859])
        ndata1 = (data1 / (sum(data1) *1.0) )
        ndata2 = (data2 / (sum(data2) *1.0) )

        dotprod = sum([ (a*b) for (a,b) in zip(ndata1, ndata2) ])
        lenx = sqrt(sum([ (a*a) for (a,b) in zip(ndata1, ndata2) ]))
        leny = sqrt(sum([ (b*b) for (a,b) in zip(ndata1, ndata2) ]))

        math.acos(dotprod/(lenx*leny))
        # 1.483262002242929
        res = [ (a-b)*(a-b) for (a,b) in zip(ndata1, ndata2) ]
        sqrt(sum(res)/len(data1))
        # 0.67272266738875497

        import scipy.stats.stats
        scipy.stats.stats.pearsonr(ndata1, ndata2)
        # (-0.65459131605877441, 0.54568145960752545)

        deltas = [ abs(a-b) for (a,b) in zip(ndata1, ndata2) ]
        sum(deltas) / len(data1)
        #0.5800337593857342

        sqrtdata1 = sqrt(data1)
        sqrtdata2 = sqrt(data2)
        norm1 = sqrtdata1 / sqrt( sum([s*s for s in sqrtdata1]) )
        norm2 = sqrtdata2 / sqrt( sum([s*s for s in sqrtdata2]) )
        sum([ (a*b) for (a,b) in zip(norm1, norm2) ])
        # 0.34514800971521764

        ndata1 = (data1 / (sum(data1) *1.0) )
        ndata2 = (data2 / (sum(data2) *1.0) )

        nsqrtdata1 = (sqrtdata1 / (sum(sqrtdata1) *1.0) )
        nsqrtdata2 = (sqrtsdata2 / (sum(sqrtdata2) *1.0) )
        sum([ abs(a-b) for (a,b) in zip(nsqrtdata1, nsqrtdata2) ])
        # 1.2796447146892949

          */

          MockMRMFeature * imrmfeature = new MockMRMFeature();

          // create mrmfeature, add "experimental" intensities
          boost::shared_ptr<MockFeature> f1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
          boost::shared_ptr<MockFeature> f2_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
          boost::shared_ptr<MockFeature> f3_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
          f1_ptr->m_intensity = (float)782.38073;
          f2_ptr->m_intensity = (float)58.384506;
          f3_ptr->m_intensity = (float)58.384506;
          std::map<std::string, boost::shared_ptr<MockFeature> > features;
          features["group1"] = f1_ptr;
          features["group2"] = f2_ptr;
          features["group3"] = f2_ptr;
          imrmfeature->m_features = features;

          // create transitions, e.g. library intensity
          std::vector<OpenSwath::LightTransition> transitions;
          { OpenSwath::LightTransition t; t.library_intensity = 1;     t.transition_name = "group1"; transitions.push_back(t); }
          { OpenSwath::LightTransition t; t.library_intensity = 10000; t.transition_name = "group2"; transitions.push_back(t); }
          { OpenSwath::LightTransition t; t.library_intensity = 2000;  t.transition_name = "group3"; transitions.push_back(t); }

          MRMScoring mrmscore;
          double manhatten, dotproduct;
          double spectral_angle, rmsd;
          double library_corr, library_rmsd;
          mrmscore.calcLibraryScore(imrmfeature, transitions, library_corr, library_rmsd, manhatten, dotproduct, spectral_angle, rmsd);
          TEST_REAL_SIMILAR(library_corr, -0.654591316)
          TEST_REAL_SIMILAR(library_rmsd, 0.5800337593)

          TEST_REAL_SIMILAR(manhatten, 1.279644714)
          TEST_REAL_SIMILAR(dotproduct, 0.34514801)

          TEST_REAL_SIMILAR(spectral_angle, 1.483262)
          TEST_REAL_SIMILAR(rmsd, 0.6727226674)

          delete imrmfeature;
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_RT_score)
        {
          MRMScoring mrmscore;
          OpenSwath::LightCompound pep;
          pep.rt = 100;
          TEST_REAL_SIMILAR(mrmscore.calcRTScore(pep, 100), 0)
          TEST_REAL_SIMILAR(mrmscore.calcRTScore(pep, 0), 100)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_SN_score)
        {
          MRMScoring mrmscore;
          std::vector<OpenSwath::ISignalToNoisePtr> sn_estimators;
          boost::shared_ptr<MockSignalToNoise> sn1 = boost::shared_ptr<MockSignalToNoise>(new MockSignalToNoise());
          sn1->m_sn_value = 500;
          boost::shared_ptr<MockSignalToNoise> sn2 = boost::shared_ptr<MockSignalToNoise>(new MockSignalToNoise());
          sn2->m_sn_value = 1500;
          sn_estimators.push_back(sn1);
          sn_estimators.push_back(sn2);

          MockMRMFeature imrmfeature;
          boost::shared_ptr<MockFeature> f1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
          boost::shared_ptr<MockFeature> f2_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
          f1_ptr->m_rt = 1200;
          f2_ptr->m_rt = 1200;
          std::map<std::string, boost::shared_ptr<MockFeature> > features;
          features["group1"] = f1_ptr;
          features["group2"] = f2_ptr;
          imrmfeature.m_features = features;

          TEST_REAL_SIMILAR(mrmscore.calcSNScore(&imrmfeature, sn_estimators), 1000.0)
          TEST_REAL_SIMILAR(mrmscore.calcSeparateSNScore(&imrmfeature, sn_estimators)[0], 6.21461)
          TEST_REAL_SIMILAR(mrmscore.calcSeparateSNScore(&imrmfeature, sn_estimators)[1], 7.31322)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(initializeMIMatrix)
        {
/*
* Requires Octave with installed MIToolbox

y = [5.97543668746948 4.2749171257019 3.3301842212677 4.08597040176392 5.50307035446167 5.24326848983765 8.40812492370605 2.83419919013977 6.94378805160522 7.69957494735718 4.08597040176392]';
x = [15.8951349258423 41.5446395874023 76.0746307373047 109.069435119629 111.90364074707 169.79216003418 121.043930053711 63.0136985778809 44.6150207519531 21.4926776885986 7.93575811386108]';

[~, ~, y_ranking] = unique(y);
[~, ~, x_ranking] = unique(x);

% test_calcMIScore matrices
m1 = [mi(x_ranking,y_ranking) mi(y_ranking,y_ranking) mi(x_ranking,x_ranking)]
mean(m1)

% test_calcSeparateMIContrastScore
m2 = zeros(2,2)
m2(1,1) = mi(x_ranking,y_ranking)
m2(2,1) = mi(y_ranking,y_ranking)
m2(1,2) = mi(x_ranking,x_ranking)
m2(2,2) = mi(y_ranking,x_ranking)
mean(m2)

% test_calcMIWeightedScore
m3 = [mi(x_ranking,y_ranking)*0.5*0.5 mi(y_ranking,y_ranking)*0.5*0.5*2 mi(x_ranking,x_ranking)*0.5*0.5]
sum(m3)

% test_calcMIPrecursorContrastScore
ms1 = [0.0 110.0 200.0 270.0 320.0 350.0 360.0 350.0 320.0 270.0 200.0]'
[~, ~, ms1_ranking] = unique(ms1);

m4 = [mi(x_ranking,ms1_ranking) mi(y_ranking,ms1_ranking)]
mean(m4)

*/

          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);

          //initialize the MI Matrix
          mrmscore.initializeMIMatrix(imrmfeature, native_ids);
          delete imrmfeature;

          TEST_EQUAL(mrmscore.getMIMatrix().rows(), 2)
          TEST_EQUAL(mrmscore.getMIMatrix().cols(), 2)

          TEST_REAL_SIMILAR(mrmscore.getMIMatrix()(0, 0), 3.2776)
          TEST_REAL_SIMILAR(mrmscore.getMIMatrix()(0, 1), 3.2776)
          TEST_REAL_SIMILAR(mrmscore.getMIMatrix()(1, 1), 3.4594)
          TEST_REAL_SIMILAR(mrmscore.getMIMatrix()(1, 0), 0) // value not initialized for lower diagonal half of matrix
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(initializeMIPrecursorContrastMatrix)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeMIPrecursorContrastMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_EQUAL(mrmscore.getMIPrecursorContrastMatrix().rows(), 3)
          TEST_EQUAL(mrmscore.getMIPrecursorContrastMatrix().cols(), 2)
          double sum = mrmscore.getMIPrecursorContrastMatrix().getEigenMatrix().sum();
          TEST_REAL_SIMILAR(sum, 12.01954465)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(initializeMIPrecursorCombinedMatrix)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeMIPrecursorCombinedMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_EQUAL(mrmscore.getMIPrecursorCombinedMatrix().rows(), 5)
          TEST_EQUAL(mrmscore.getMIPrecursorCombinedMatrix().cols(), 5)

          double sum = mrmscore.getMIPrecursorCombinedMatrix().getEigenMatrix().sum();
          TEST_REAL_SIMILAR(sum, 48.98726953)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(initializeMIContrastMatrix)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> native_ids1, native_ids2;
          fill_mock_objects(imrmfeature, native_ids1);
          for (int i=native_ids1.size()-1; i>=0; i--)
          {
            native_ids2.push_back(native_ids1[i]);
          }   

          //initialize the XCorr Matrix
          mrmscore.initializeMIContrastMatrix(imrmfeature, native_ids1, native_ids2);
          delete imrmfeature;

          TEST_REAL_SIMILAR(mrmscore.getMIContrastMatrix()(0, 0), 3.2776)
          TEST_REAL_SIMILAR(mrmscore.getMIContrastMatrix()(0, 1), 3.2776)
          TEST_REAL_SIMILAR(mrmscore.getMIContrastMatrix()(1, 1), 3.2776)
          TEST_REAL_SIMILAR(mrmscore.getMIContrastMatrix()(1, 0), 3.4594)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcMIScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;
          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);
          mrmscore.initializeMIMatrix(imrmfeature, native_ids);
          delete imrmfeature;
          TEST_REAL_SIMILAR(mrmscore.calcMIScore(), 3.3382) // mean + std deviation
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcSeparateMIContrastScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;
          std::vector<std::string> native_ids1, native_ids2;
          fill_mock_objects(imrmfeature, native_ids1);
          for (int i=native_ids1.size()-1; i>=0; i--)
          {
            native_ids2.push_back(native_ids1[i]);
          }   
          mrmscore.initializeMIContrastMatrix(imrmfeature, native_ids1, native_ids2);
          delete imrmfeature;
          TEST_REAL_SIMILAR(mrmscore.calcSeparateMIContrastScore()[0], 3.27761)
          TEST_REAL_SIMILAR(mrmscore.calcSeparateMIContrastScore()[1], 3.36852)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcMIWeightedScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;
          std::vector<std::string> native_ids;
          fill_mock_objects(imrmfeature, native_ids);
          static const double weights_[] = { 0.5, 0.5 };
          std::vector<double> weights (weights_, weights_ + sizeof(weights_) / sizeof(weights_[0]) );
          mrmscore.initializeMIMatrix(imrmfeature, native_ids);
          delete imrmfeature;

          // xcorr_deltas = [1, 0.3969832, 1] * array([0.25, 2*0.5*0.5,0.25])
          // sum(xcorr_deltas)
          TEST_REAL_SIMILAR(mrmscore.calcMIWeightedScore(weights), 3.3231)
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcMIPrecursorContrastScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeMIPrecursorContrastMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_REAL_SIMILAR(mrmscore.calcMIPrecursorContrastScore(), 2.003257 )
        }
    END_SECTION

    BOOST_AUTO_TEST_CASE(test_calcMIPrecursorCombinedScore)
        {
          MockMRMFeature * imrmfeature = new MockMRMFeature();
          MRMScoring mrmscore;

          std::vector<std::string> precursor_ids;
          std::vector<std::string> native_ids;
          fill_mock_objects2(imrmfeature, precursor_ids, native_ids);

          //initialize the XCorr vector
          mrmscore.initializeMIPrecursorCombinedMatrix(imrmfeature, precursor_ids, native_ids);
          delete imrmfeature;

          TEST_REAL_SIMILAR(mrmscore.calcMIPrecursorCombinedScore(), 1.959490 )
        }
    END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
