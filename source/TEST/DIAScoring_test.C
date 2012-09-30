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

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenMSHelper.h>

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/MockObjects.h"

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;

void getMRMFeatureTest(MockMRMFeature * imrmfeature_test)
{
  boost::shared_ptr<MockFeature> f1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> f2_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  f1_ptr->m_intensity = 0.3;
  f2_ptr->m_intensity = 0.7;
  std::map<std::string, boost::shared_ptr<MockFeature> > features;
  features["group1"] = f1_ptr;
  features["group2"] = f2_ptr;
  imrmfeature_test->m_features = features;
  imrmfeature_test->m_intensity = 1.0;
}

START_TEST(DIAScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DIAScoring* ptr = 0;
DIAScoring* nullPointer = 0;

START_SECTION(DIAScoring())
{
  ptr = new DIAScoring();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~DIAScoring())
{
  delete ptr;
}
END_SECTION

START_SECTION((void MRMFeatureScoring::standardize_data(std::vector<double>& data)))
{
// see seperate test
NOT_TESTABLE
}
END_SECTION

START_SECTION((void MRMFeatureScoring::getBYSeries(AASequence& a, int charge, std::vector<double>& bseries, std::vector<double>& yseries)))
{

  OpenSwath::DIAScoring diascoring;
  String sequence = "SYVAWDR";
  std::vector<double> bseries, yseries;
  OpenMS::AASequence a = OpenMS::AASequence(sequence);
  OpenSwath::getBYSeries(a,  bseries, yseries, 1);

  TEST_EQUAL(bseries.size(), 5)
  TEST_EQUAL(yseries.size(), 6)

  //TEST_REAL_SIMILAR (bseries[0],  88.03990 );
  TEST_REAL_SIMILAR (bseries[0], 251.10323 );
  TEST_REAL_SIMILAR (bseries[1], 350.17164 );
  TEST_REAL_SIMILAR (bseries[2], 421.20875 );
  TEST_REAL_SIMILAR (bseries[3], 607.28807 );
  TEST_REAL_SIMILAR (bseries[4], 722.31501 );
  //TEST_REAL_SIMILAR (bseries[5], 878.41612 );

  TEST_REAL_SIMILAR (yseries[0], 175.11955  );
  TEST_REAL_SIMILAR (yseries[1], 290.14649  );
  TEST_REAL_SIMILAR (yseries[2], 476.22580  );
  TEST_REAL_SIMILAR (yseries[3], 547.26291  );
  TEST_REAL_SIMILAR (yseries[4], 646.33133  );
  TEST_REAL_SIMILAR (yseries[5], 809.39466  );
  //TEST_REAL_SIMILAR (yseries[6], 896.42668  );

  // now add a modification to the sequence
  bseries.clear();
  yseries.clear();
  a.setModification(1, "Phospho" ); // modify the Y
  OpenSwath::getBYSeries(a,  bseries, yseries ,1);

  TEST_EQUAL(bseries.size(), 5)
  TEST_EQUAL(yseries.size(), 6)

  //TEST_REAL_SIMILAR (bseries[0],  88.03990 );
  TEST_REAL_SIMILAR (bseries[0], 251.10323 + 79.9657 );
  TEST_REAL_SIMILAR (bseries[1], 350.17164 + 79.9657 );
  TEST_REAL_SIMILAR (bseries[2], 421.20875 + 79.9657 );
  TEST_REAL_SIMILAR (bseries[3], 607.28807 + 79.9657 );
  TEST_REAL_SIMILAR (bseries[4], 722.31501 + 79.9657 );
  //TEST_REAL_SIMILAR (bseries[5], 878.41612 );

  TEST_REAL_SIMILAR (yseries[0], 175.11955  );
  TEST_REAL_SIMILAR (yseries[1], 290.14649  );
  TEST_REAL_SIMILAR (yseries[2], 476.22580  );
  TEST_REAL_SIMILAR (yseries[3], 547.26291  );
  TEST_REAL_SIMILAR (yseries[4], 646.33133  );
  TEST_REAL_SIMILAR (yseries[5], 809.39466 + 79.9657);
  //TEST_REAL_SIMILAR (yseries[6], 896.42668  );
    
}
END_SECTION

OpenSwath::LightTransition mock_tr1; 
mock_tr1.product_mz = 500;
mock_tr1.charge = 1;
mock_tr1.transition_name = "group1";

OpenSwath::LightTransition mock_tr2; 
mock_tr2.product_mz = 600;
mock_tr2.charge = 1;
mock_tr2.transition_name = "group2";

START_SECTION ( forward
void dia_isotope_scores(const std::vector<TransitionType> & transitions,
                       SpectrumType  spectrum, OpenSwath::IMRMFeature * mrmfeature, int putative_fragment_charge,
                       double & isotope_corr, double & isotope_overlap))
{
  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2(new OpenSwath::BinaryDataArray);

  static const double arr1[] = {
    /*
    10, 20, 50, 100, 50, 20, 10, // peak at 499
    3, 7, 15, 30, 15, 7, 3,      // peak at 500
    1, 3, 9, 15, 9, 3, 1,        // peak at 501
    3, 9, 3,                     // peak at 502
    */

    10, 20, 50, 100, 50, 20, 10, // peak at 600
    3, 7, 15, 30, 15, 7, 3,      // peak at 601
    1, 3, 9, 15, 9, 3, 1,        // peak at 602
    3, 9, 3                      // peak at 603
  };
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = {
    /*
    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03,
    499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03,
    500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
    501.99, 502.0, 502.01,
    */
    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
    600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  data1->data = mz;
  data2->data = intensity;
  binaryDataArrayPtrs.push_back(data1);
  binaryDataArrayPtrs.push_back(data2);
  sptr->binaryDataArrayPtrs = binaryDataArrayPtrs;

  MockMRMFeature * imrmfeature_test = new MockMRMFeature();
  getMRMFeatureTest(imrmfeature_test);
  imrmfeature_test->m_intensity = 0.7;
  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr2);

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use 50 ppm and a cutoff of 30 in intensity 
  double isotope_corr = 0, isotope_overlap = 0;
  diascoring.dia_isotope_scores(transitions, sptr, imrmfeature_test, isotope_corr, isotope_overlap);
  // >>> exp = [240, 74, 39, 15, 0]
  // >>> theo = [1, 0.325757771553019, 0.0678711748364005, 0.0105918703087134, 0.00134955223787482]
  // >>> from scipy.stats.stats import pearsonr
  // >>> pearsonr(exp, theo)
  // (0.99463189043051314, 0.00047175434098498532)
  //
  TEST_REAL_SIMILAR(isotope_corr, 0.995361286111832)
  TEST_REAL_SIMILAR(isotope_overlap, 0.0)
}
END_SECTION

START_SECTION ( backward
void dia_isotope_scores(const std::vector<TransitionType> & transitions,
                       SpectrumType  spectrum, OpenSwath::IMRMFeature * mrmfeature, int putative_fragment_charge,
                       double & isotope_corr, double & isotope_overlap))
{
  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);

  static const double arr1[] = {
    10, 20, 50, 100, 50, 20, 10, // peak at 499
    3, 7, 15, 30, 15, 7, 3,      // peak at 500
    1, 3, 9, 15, 9, 3, 1,        // peak at 501
    3, 9, 3                      // peak at 502

    /*
    10, 20, 50, 100, 50, 20, 10, // peak at 600
    3, 7, 15, 30, 15, 7, 3,      // peak at 601
    1, 3, 9, 15, 9, 3, 1,        // peak at 602
    3, 9, 3                      // peak at 603
    */
  };
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = {
    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03,
    499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03,
    500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
    501.99, 502.0, 502.01

    /*
    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
    600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
    */
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  data1->data = mz;
  data2->data = intensity;
  binaryDataArrayPtrs.push_back(data1);
  binaryDataArrayPtrs.push_back(data2);
  sptr->binaryDataArrayPtrs = binaryDataArrayPtrs;

  MockMRMFeature * imrmfeature_test = new MockMRMFeature();
  getMRMFeatureTest(imrmfeature_test);
  imrmfeature_test->m_intensity = 0.3;
  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr1);

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use 50 ppm and a cutoff of 30 in intensity 
  double isotope_corr = 0, isotope_overlap = 0;
  diascoring.dia_isotope_scores(transitions, sptr, imrmfeature_test, isotope_corr, isotope_overlap);

  // >>> exp = [74, 39, 15, 0, 0]
  // >>> theo = [1, 0.266799519434277, 0.0486475002325161, 0.0066525896497495, 0.000747236543377621]
  // >>> from scipy.stats.stats import pearsonr
  // >>> pearsonr(exp, theo)
  // (0.959570883150479, 0.0096989307464742554)
  // there is one peak (this one) which has an overlap in isotopes

  TEST_REAL_SIMILAR(isotope_corr, 0.959570883150479)
  TEST_REAL_SIMILAR(isotope_overlap, 1.0)

}
END_SECTION

START_SECTION (
void dia_isotope_scores(const std::vector<TransitionType> & transitions,
                       SpectrumType  spectrum, OpenSwath::IMRMFeature * mrmfeature, int putative_fragment_charge,
                       double & isotope_corr, double & isotope_overlap))
{
  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);

  static const double arr1[] = {
    10, 20, 50, 100, 50, 20, 10, // peak at 499
    3, 7, 15, 30, 15, 7, 3,      // peak at 500
    1, 3, 9, 15, 9, 3, 1,        // peak at 501
    3, 9, 3,                     // peak at 502

    10, 20, 50, 100, 50, 20, 10, // peak at 600
    3, 7, 15, 30, 15, 7, 3,      // peak at 601
    1, 3, 9, 15, 9, 3, 1,        // peak at 602
    3, 9, 3                      // peak at 603
  };
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = {
    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03,
    499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03,
    500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
    501.99, 502.0, 502.01,

    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
    600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  data1->data = mz;
  data2->data = intensity;
  binaryDataArrayPtrs.push_back(data1);
  binaryDataArrayPtrs.push_back(data2);
  sptr->binaryDataArrayPtrs = binaryDataArrayPtrs;

  MockMRMFeature * imrmfeature_test = new MockMRMFeature();
  getMRMFeatureTest(imrmfeature_test);

  // create transitions, e.g. library intensity
  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr1);
  transitions.push_back(mock_tr2);

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use 50 ppm and a cutoff of 30 in intensity 
  double isotope_corr = 0, isotope_overlap = 0;
  diascoring.dia_isotope_scores(transitions, sptr, imrmfeature_test, isotope_corr, isotope_overlap);

  // see above for the the two individual numbers (forward and backward)
  TEST_REAL_SIMILAR(isotope_corr, 0.984624164796771)
  TEST_REAL_SIMILAR(isotope_overlap, 1.0 * 0.3)

}
END_SECTION

START_SECTION (
void dia_massdiff_score(const std::vector<TransitionType> & transitions,
                       SpectrumType  spectrum, const std::vector<double> & normalized_library_intensity,
                       double & ppm_score, double & ppm_score_weighted))
{ 
  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);

  static const double arr1[] = {
    10, 20, 50, 100, 50, 20, 10, // peak at 499
    3, 7, 15, 30, 15, 7, 3,      // peak at 500
    1, 3, 9, 15, 9, 3, 1,        // peak at 501
    3, 9, 3,                     // peak at 502

    10, 20, 50, 100, 50, 20, 10, // peak at 600
    3, 7, 15, 30, 15, 7, 3,      // peak at 601
    1, 3, 9, 15, 9, 3, 1,        // peak at 602
    3, 9, 3                      // peak at 603
  };
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = {
    498.97, 498.98, 498.99, 499.0, 499.01, 499.02, 499.03,
    499.97, 499.98, 499.99, 500.0, 500.01, 500.02, 500.03,
    500.97, 500.98, 500.99, 501.0, 501.01, 501.02, 501.03,
    501.99, 502.0, 502.01,

    599.97, 599.98, 599.99, 600.0, 600.01, 600.02, 600.03,
    600.97, 600.98, 600.99, 601.0, 601.01, 601.02, 601.03,
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  // shift the peaks by a fixed amount in ppm
  for (std::size_t i = 0; i < mz.size() / 2.0; i++)
  {
    mz[i] +=  mz[i] / 1000000 * 15; // shift first peak by 15 ppm
    //std::cout << " new mz " << mz[i] << std::endl;
  }
  for (std::size_t i = mz.size() / 2.0; i < mz.size(); i++)
  {
    mz[i] +=  mz[i] / 1000000 * 10; // shift second peak by 10 ppm
    //std::cout << " new mz " << mz[i] << std::endl;
  }
  data1->data = mz;
  data2->data = intensity;
  binaryDataArrayPtrs.push_back(data1);
  binaryDataArrayPtrs.push_back(data2);
  sptr->binaryDataArrayPtrs = binaryDataArrayPtrs;

  MockMRMFeature * imrmfeature_test = new MockMRMFeature();
  getMRMFeatureTest(imrmfeature_test);

  // create transitions, e.g. library intensity
  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr1);
  transitions.push_back(mock_tr2);

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.5, false, 30, 50, 4, 4); // here we use a large enough window so that none of our peaks falls out
  double ppm_score = 0, ppm_score_weighted = 0;
  std::vector<double> normalized_library_intensity;
  normalized_library_intensity.push_back(0.7);
  normalized_library_intensity.push_back(0.3);
  diascoring.dia_massdiff_score(transitions, sptr, normalized_library_intensity, ppm_score, ppm_score_weighted); 

  TEST_REAL_SIMILAR(ppm_score, 15 + 10); // 15 ppm and 10 ppm 
  TEST_REAL_SIMILAR(ppm_score_weighted, 15 * 0.7 + 10* 0.3); // weighted
}
END_SECTION

START_SECTION ( void dia_by_ion_score(SpectrumType & spectrum, AASequence & sequence, int charge, double & bseries_score, double & yseries_score))
{

  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);

  static const double arr1[] = {
    100, 100, 100, 100,
    100, 100, 100
  };
  std::vector<double> intensity (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = {
    // four of the naked b/y ions 
    // as well as one of the modified b and y ions ion each
    350.17164, // b
    421.20875, // b
    421.20875 + 79.9657, // b + P
    547.26291, // y
    646.33133, // y
    809.39466 + 79.9657 // y + P
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  data1->data = mz;
  data2->data = intensity;
  binaryDataArrayPtrs.push_back(data1);
  binaryDataArrayPtrs.push_back(data2);
  sptr->binaryDataArrayPtrs = binaryDataArrayPtrs;

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use a large enough window so that none of our peaks falls out
  String sequence = "SYVAWDR";
  std::vector<double> bseries, yseries;
  AASequence a = AASequence(sequence);

  double bseries_score = 0, yseries_score = 0;
  diascoring.dia_by_ion_score(sptr, a, 1, bseries_score, yseries_score); 

  TEST_REAL_SIMILAR (bseries_score, 2);
  TEST_REAL_SIMILAR (yseries_score, 2);

  // now add a modification to the sequence
  a.setModification(1, "Phospho" ); // modify the Y
  bseries_score = 0, yseries_score = 0;
  diascoring.dia_by_ion_score(sptr, a, 1, bseries_score, yseries_score); 

  TEST_REAL_SIMILAR (bseries_score, 1);
  TEST_REAL_SIMILAR (yseries_score, 3);
}
END_SECTION

START_SECTION (void integrateWindows(const SpectrumType spectrum, const double & mz_start, const double & mz_end, double & mz, double & intensity, bool centroided))
{
  // TODO is tested above, maybe seperate test here?
}
END_SECTION

START_SECTION((void set_dia_parameters(double dia_extract_window, double dia_centroided, double dia_byseries_intensity_min, double dia_byseries_ppm_diff, double dia_nr_isotopes, double dia_nr_charges)))
{
  NOT_TESTABLE
}
END_SECTION
START_SECTION((void score_with_isotopes(SpectrumType spectrum, const std::vector<TransitionType> & transitions, double & dotprod, double & manhattan);))
// TODO (wolski): write tests
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

