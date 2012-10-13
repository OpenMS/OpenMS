// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/ALGO/MRMScoring.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/MockObjects.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h"

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

void fill_mock_objects(MockMRMFeature * imrmfeature, MockTransitionGroup * itransition_group)
{
  std::vector<std::string> native_ids;
  native_ids.push_back("group1");
  native_ids.push_back("group2");
  itransition_group->m_native_ids = native_ids;

  itransition_group->m_size = 2;
  static const double arr1[] = {
    5.97543668746948, 4.2749171257019, 3.3301842212677, 4.08597040176392, 5.50307035446167, 5.24326848983765,
    8.40812492370605, 2.83419919013977, 6.94378805160522, 7.69957494735718, 4.08597040176392 
  };
  std::vector<double> intensity1 (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = {
    15.8951349258423, 41.5446395874023, 76.0746307373047, 109.069435119629, 111.90364074707, 169.79216003418,
    121.043930053711, 63.0136985778809, 44.6150207519531, 21.4926776885986, 7.93575811386108
  };
  std::vector<double> intensity2 (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  boost::shared_ptr<MockFeature> f1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> f2_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  f1_ptr->m_intensity_vec = intensity1;
  f2_ptr->m_intensity_vec = intensity2;
  std::map<std::string, boost::shared_ptr<MockFeature> > features;
  features["group1"] = f1_ptr;
  features["group2"] = f2_ptr;
  imrmfeature->m_features = features;

}

///////////////////////////

START_TEST(MRMScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#ifndef USE_BOOST_UNIT_TEST
{
MRMScoring* ptr = 0;
MRMScoring* nullPointer = 0;

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

//START_SECTION((void initializeXCorrMatrix(MRMFeature& mrmfeature, MRMTransitionGroup< SpectrumType, PeakType >& transition_group, bool normalize)))
BOOST_AUTO_TEST_CASE(initializeXCorrMatrix)
{
  // test cross-correlation
  //
  // numpy: 
  //
  // data1 = [5.97543668746948, 4.2749171257019, 3.3301842212677, 4.08597040176392, 5.50307035446167, 5.24326848983765,
  //          8.40812492370605, 2.83419919013977, 6.94378805160522, 7.69957494735718, 4.08597040176392 ]
  // data2 = [15.8951349258423, 41.5446395874023, 76.0746307373047, 109.069435119629, 111.90364074707, 169.79216003418,
  //          121.043930053711, 63.0136985778809, 44.6150207519531, 21.4926776885986, 7.93575811386108]
  // data1 = (data1 - mean(data1) ) / std(data1)
  // data2 = (data2 - mean(data2) ) / std(data2)
  // xcorrmatrix_0_0 = correlate(data1, data1, "same") / len(data1)
  // xcorrmatrix_0_1 = correlate(data1, data2, "same") / len(data1)
  // xcorr_deltas = [0, 3, 0]
  // xcorr_coeltion = mean(deltas) + std(deltas, ddof=1)

  MockMRMFeature * imrmfeature = new MockMRMFeature();
  MockTransitionGroup * itransition_group = new MockTransitionGroup();
  MRMScoring mrmscore;

  fill_mock_objects(imrmfeature, itransition_group);

  //initialize the XCorr Matrix
  mrmscore.initializeXCorrMatrix(imrmfeature, itransition_group, true);

  TEST_EQUAL(mrmscore.getXCorrMatrix().size(), 2)
  TEST_EQUAL(mrmscore.getXCorrMatrix()[0].size(), 2)
  TEST_EQUAL(mrmscore.getXCorrMatrix()[0][0].size(), 23)

  // test auto-correlation
  const std::map<int, double> auto_correlation =
      mrmscore.getXCorrMatrix()[0][0];
  TEST_REAL_SIMILAR(auto_correlation.find(0)->second, 1)
  TEST_REAL_SIMILAR(auto_correlation.find(1)->second, -0.227352707759245)
  TEST_REAL_SIMILAR(auto_correlation.find(-1)->second, -0.227352707759245)
  TEST_REAL_SIMILAR(auto_correlation.find(2)->second, -0.07501116)
  TEST_REAL_SIMILAR(auto_correlation.find(-2)->second, -0.07501116)

  const std::map<int, double> cross_correlation =
      mrmscore.getXCorrMatrix()[0][1];
  TEST_REAL_SIMILAR(cross_correlation.find(2)->second, -0.31165141)
  TEST_REAL_SIMILAR(cross_correlation.find(1)->second, -0.35036919)
  TEST_REAL_SIMILAR(cross_correlation.find(0)->second, 0.03129565)
  TEST_REAL_SIMILAR(cross_correlation.find(-1)->second, 0.30204049)
  TEST_REAL_SIMILAR(cross_correlation.find(-2)->second, 0.13012441)
  TEST_REAL_SIMILAR(cross_correlation.find(-3)->second, 0.39698322)
  TEST_REAL_SIMILAR(cross_correlation.find(-4)->second, 0.16608774)
}
END_SECTION

BOOST_AUTO_TEST_CASE(test_calcXcorrCoelutionScore)
{
  MockMRMFeature * imrmfeature = new MockMRMFeature();
  MockTransitionGroup * itransition_group = new MockTransitionGroup();
  MRMScoring mrmscore;

  fill_mock_objects(imrmfeature, itransition_group);

  //initialize the XCorr Matrix
  mrmscore.initializeXCorrMatrix(imrmfeature, itransition_group, true);
  
  // xcorr_deltas = [0, 3, 0]
  // xcorr_coeltion = mean(xcorr_deltas) + std(xcorr_deltas, ddof=1)
  TEST_REAL_SIMILAR(mrmscore.calcXcorrCoelutionScore(), 1 + std::sqrt(3.0)) // mean + std deviation
}
END_SECTION

BOOST_AUTO_TEST_CASE(test_calcXcorrCoelutionScore_weighted)
{
  MockMRMFeature * imrmfeature = new MockMRMFeature();
  MockTransitionGroup * itransition_group = new MockTransitionGroup();
  MRMScoring mrmscore;

  static const double weights_[] = { 0.5, 0.5 };
  std::vector<double> weights (weights_, weights_ + sizeof(weights_) / sizeof(weights_[0]) );

  fill_mock_objects(imrmfeature, itransition_group);

  //initialize the XCorr Matrix
  mrmscore.initializeXCorrMatrix(imrmfeature, itransition_group, true);
  
  // xcorr_deltas = [0, 3, 0] * array([0.25, 2*0.5*0.5,0.25])
  // sum(xcorr_deltas) 
  TEST_REAL_SIMILAR(mrmscore.calcXcorrCoelutionScore_weighted(weights), 1.5) 
}
END_SECTION

BOOST_AUTO_TEST_CASE(test_calcXcorrShape_score)
{
  MockMRMFeature * imrmfeature = new MockMRMFeature();
  MockTransitionGroup * itransition_group = new MockTransitionGroup();
  MRMScoring mrmscore;

  fill_mock_objects(imrmfeature, itransition_group);

  //initialize the XCorr Matrix
  mrmscore.initializeXCorrMatrix(imrmfeature, itransition_group, true);
  
  // xcorr_deltas = [1, 0.3969832, 1]
  // xcorr_coeltion = mean(xcorr_deltas)
  TEST_REAL_SIMILAR(mrmscore.calcXcorrShape_score(), (1+0.3969832+1)/3.0) // mean + std deviation
}
END_SECTION

BOOST_AUTO_TEST_CASE(test_calcXcorrShape_score_weighted)
{
  MockMRMFeature * imrmfeature = new MockMRMFeature();
  MockTransitionGroup * itransition_group = new MockTransitionGroup();
  MRMScoring mrmscore;

  fill_mock_objects(imrmfeature, itransition_group);

  static const double weights_[] = { 0.5, 0.5 };
  std::vector<double> weights (weights_, weights_ + sizeof(weights_) / sizeof(weights_[0]) );

  //initialize the XCorr Matrix
  mrmscore.initializeXCorrMatrix(imrmfeature, itransition_group, true);
  
  // xcorr_deltas = [1, 0.3969832, 1] * array([0.25, 2*0.5*0.5,0.25])
  // sum(xcorr_deltas)
  TEST_REAL_SIMILAR(mrmscore.calcXcorrShape_score_weighted(weights), 0.6984916)
}
END_SECTION

//START_SECTION((virtual void test_calcLibraryScore()))
BOOST_AUTO_TEST_CASE(test_Library_score)
{
  // numpy
  // data1 = array([1,10000,2000])
  // data2 = array([782.380737304688, 58.3845062255859, 58.3845062255859])
  double library_corr, library_rmsd; 
  MockMRMFeature * imrmfeature = new MockMRMFeature();

  // create mrmfeature, add "experimental" intensities
  boost::shared_ptr<MockFeature> f1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> f2_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> f3_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  f1_ptr->m_intensity = 782.38073;
  f2_ptr->m_intensity = 58.384506;
  f3_ptr->m_intensity = 58.384506;
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
  double d1, d2;
  mrmscore.calcLibraryScore(imrmfeature, transitions, library_corr, library_rmsd, d1, d2);
  TEST_REAL_SIMILAR(library_corr, -0.654591316)
  TEST_REAL_SIMILAR(library_rmsd, 0.5800337593)

  delete imrmfeature;
}
END_SECTION

BOOST_AUTO_TEST_CASE(test_RT_score)
{
  MRMScoring mrmscore;
  OpenSwath::LightPeptide pep;
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
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
