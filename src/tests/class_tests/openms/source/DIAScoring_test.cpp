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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DIAHelper.h>

#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h"
#include "OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/MockObjects.h"

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;

void getMRMFeatureTest(MockMRMFeature * imrmfeature_test)
{
  boost::shared_ptr<MockFeature> f1_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  boost::shared_ptr<MockFeature> f2_ptr = boost::shared_ptr<MockFeature>(new MockFeature());
  f1_ptr->m_intensity = 0.3f;
  f2_ptr->m_intensity = 0.7f;
  std::map<std::string, boost::shared_ptr<MockFeature> > features;
  features["group1"] = f1_ptr;
  features["group2"] = f2_ptr;
  imrmfeature_test->m_features = features;
  imrmfeature_test->m_intensity = 1.0;
}

OpenSwath::SpectrumPtr prepareSpectrum()
{
  OpenSwath::SpectrumPtr sptr = (OpenSwath::SpectrumPtr)(new OpenSwath::Spectrum);
  std::vector<OpenSwath::BinaryDataArrayPtr> binaryDataArrayPtrs;
  OpenSwath::BinaryDataArrayPtr data1 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr data2 = (OpenSwath::BinaryDataArrayPtr)(new OpenSwath::BinaryDataArray);

  static const double arr1[] = {
    10, 20, 50, 100, 50, 20, 10, // peak at 499 -> 260-20 = 240 intensity within 0.05 Th
    3, 7, 15, 30, 15, 7, 3,      // peak at 500 -> 80-6 = 74 intensity within 0.05 Th
    1, 3, 9, 15, 9, 3, 1,        // peak at 501 -> 41-2 = 39 intensity within 0.05 Th
    3, 9, 3,                     // peak at 502 -> 15 intensity within 0.05 Th

    10, 20, 50, 100, 50, 20, 10, // peak at 600 -> 260-20 = 240 intensity within 0.05 Th
    3, 7, 15, 30, 15, 7, 3,      // peak at 601 -> 80-6 = 74 intensity within 0.05 Th
    1, 3, 9, 15, 9, 3, 1,        // peak at 602 -> sum([ 9, 15, 9, 3, 1]) = 37 intensity within 0.05 Th
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
    // note that this peak at 602 is special since it is integrated from 
    // [(600+2*1.0033548) - 0.025, (600+2*1.0033548)  + 0.025] = [601.9817096 to 602.0317096]
    601.97, 601.98, 601.99, 602.0, 602.01, 602.02, 602.03,
    602.99, 603.0, 603.01
  };
  std::vector<double> mz (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );
  data1->data = mz;
  data2->data = intensity;

  sptr->setMZArray(data1);
  sptr->setIntensityArray( data2 );
  return sptr;
}

OpenSwath::SpectrumPtr prepareShiftedSpectrum()
{
  OpenSwath::SpectrumPtr sptr = prepareSpectrum();
  // shift the peaks by a fixed amount in ppm
  for (std::size_t i = 0; i < sptr->getMZArray()->data.size() / 2.0; i++)
  {
    sptr->getMZArray()->data[i] +=  sptr->getMZArray()->data[i] / 1000000 * 15; // shift first peak by 15 ppm
  }
  for (std::size_t i = sptr->getMZArray()->data.size() / 2.0; i < sptr->getMZArray()->data.size(); i++)
  {
    sptr->getMZArray()->data[i] +=  sptr->getMZArray()->data[i] / 1000000 * 10; // shift second peak by 10 ppm
  }
  return sptr;
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

START_SECTION(([EXTRA] void MRMFeatureScoring::getBYSeries(AASequence& a, int charge, std::vector<double>& bseries, std::vector<double>& yseries)))
{
  OpenMS::DIAScoring diascoring;
  String sequence = "SYVAWDR";
  std::vector<double> bseries, yseries;
  OpenMS::AASequence a = OpenMS::AASequence::fromString(sequence);
  OpenMS::DIAHelpers::getBYSeries(a,  bseries, yseries, 1);

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
  OpenMS::DIAHelpers::getBYSeries(a,  bseries, yseries ,1);

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
mock_tr1.product_charge = 1;
mock_tr1.transition_name = "group1";

OpenSwath::LightTransition mock_tr2;
mock_tr2.product_mz = 600;
mock_tr2.product_charge = 1;
mock_tr2.transition_name = "group2";

START_SECTION([EXTRA] forward void dia_isotope_scores(const std::vector<TransitionType> & transitions, SpectrumType  spectrum, OpenSwath::IMRMFeature * mrmfeature, int putative_fragment_charge, double & isotope_corr, double & isotope_overlap))
{
  OpenSwath::SpectrumPtr sptr = prepareSpectrum();
  MockMRMFeature * imrmfeature_test = new MockMRMFeature();
  getMRMFeatureTest(imrmfeature_test);
  imrmfeature_test->m_intensity = 0.7f;
  std::vector<OpenSwath::LightTransition> transitions;

  // Try with transition at 600 m/z 
  transitions.push_back(mock_tr2);

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use 50 ppm and a cutoff of 30 in intensity
  double isotope_corr = 0, isotope_overlap = 0;
  diascoring.dia_isotope_scores(transitions, sptr, imrmfeature_test, isotope_corr, isotope_overlap);

  // >> exp = [240, 74, 37, 15, 0]
  // >> theo = [1, 0.325757771553019, 0.0678711748364005, 0.0105918703087134, 0.00134955223787482]
  // >> from scipy.stats.stats import pearsonr
  // >> pearsonr(exp, theo)
  // (0.99536128611183172, 0.00037899006151919545)
  //
  TEST_REAL_SIMILAR(isotope_corr, 0.995361286111832)
  TEST_REAL_SIMILAR(isotope_overlap, 0.0)
}
END_SECTION

START_SECTION([EXTRA] backward void dia_isotope_scores(const std::vector<TransitionType> & transitions, SpectrumType  spectrum, OpenSwath::IMRMFeature * mrmfeature, int putative_fragment_charge, double & isotope_corr, double & isotope_overlap))
{
  OpenSwath::SpectrumPtr sptr = prepareSpectrum();
  MockMRMFeature * imrmfeature_test = new MockMRMFeature();
  getMRMFeatureTest(imrmfeature_test);
  imrmfeature_test->m_intensity = 0.3f;
  std::vector<OpenSwath::LightTransition> transitions;

  // Try with transition at 500 m/z 
  // This peak is not monoisotopic (e.g. at 499 there is another, more intense, peak)
  transitions.push_back(mock_tr1);

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use 50 ppm and a cutoff of 30 in intensity
  double isotope_corr = 0, isotope_overlap = 0;
  diascoring.dia_isotope_scores(transitions, sptr, imrmfeature_test, isotope_corr, isotope_overlap);

  // >> exp = [74, 39, 15, 0, 0]
  // >> theo = [1, 0.266799519434277, 0.0486475002325161, 0.0066525896497495, 0.000747236543377621]
  // >> from scipy.stats.stats import pearsonr
  // >> pearsonr(exp, theo)
  // (0.959570883150479, 0.0096989307464742554)
  TEST_REAL_SIMILAR(isotope_corr, 0.959570883150479)
  TEST_REAL_SIMILAR(isotope_overlap, 1.0)
}
END_SECTION

START_SECTION ( void dia_isotope_scores(const std::vector< TransitionType > &transitions, SpectrumType spectrum, OpenSwath::IMRMFeature *mrmfeature, double &isotope_corr, double &isotope_overlap) )
{
  OpenSwath::SpectrumPtr sptr = prepareSpectrum();

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

  // see above for the two individual numbers (forward and backward)
  TEST_REAL_SIMILAR(isotope_corr, 0.995361286111832 * 0.7 + 0.959570883150479 * 0.3)
  TEST_REAL_SIMILAR(isotope_overlap, 0.0 * 0.7 + 1.0 * 0.3)

}
END_SECTION

START_SECTION(void dia_ms1_isotope_scores(double precursor_mz, SpectrumPtrType spectrum, size_t charge_state, 
                                double& isotope_corr, double& isotope_overlap))
{
  OpenSwath::SpectrumPtr sptr = prepareSpectrum();

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use 50 ppm and a cutoff of 30 in intensity

  // Check for charge 1+ and m/z at 500
  {
    size_t precursor_charge_state = 1;
    double precursor_mz = 500;

    double isotope_corr = 0, isotope_overlap = 0;
    diascoring.dia_ms1_isotope_scores(precursor_mz, sptr, precursor_charge_state, isotope_corr, isotope_overlap);

    // see above for the two individual numbers (forward and backward)
    TEST_REAL_SIMILAR(isotope_corr, 0.959570883150479)
    TEST_REAL_SIMILAR(isotope_overlap, 240/74.0)
  }

  // Check if charge state is assumed 2+
  {
    size_t precursor_charge_state = 2;
    double precursor_mz = 500;

    double isotope_corr = 0, isotope_overlap = 0;
    diascoring.dia_ms1_isotope_scores(precursor_mz, sptr, precursor_charge_state, isotope_corr, isotope_overlap);

    // >>> theo = [0.57277789564886, 0.305415548811564, 0.0952064968352544, 0.0218253361702587, 0.00404081869309618]
    // >>> exp = [74, 0, 39, 0, 15]
    // >>> pearsonr(exp, theo)
    // (0.68135233883093205, 0.20528953804781694)
    TEST_REAL_SIMILAR(isotope_corr, 0.681352338830933)
    TEST_REAL_SIMILAR(isotope_overlap, 240/74.0)
  }

  // Check and confirm that monoisotopic is at m/z 499
  {
    size_t precursor_charge_state = 1;
    double precursor_mz = 499;

    double isotope_corr = 0, isotope_overlap = 0;
    diascoring.dia_ms1_isotope_scores(precursor_mz, sptr, precursor_charge_state, isotope_corr, isotope_overlap);

    // >> exp = [240, 74, 39, 15, 0]
    // >> theo = [0.755900817146293, 0.201673974754608, 0.0367726851778834, 0.00502869795238462, 0.000564836713740715]
    // >> from scipy.stats.stats import pearsonr
    // >> pearsonr(exp, theo)
    // (0.99463189043051314, 0.00047175434098498532)
    TEST_REAL_SIMILAR(isotope_corr, 0.995485552148335)
    TEST_REAL_SIMILAR(isotope_overlap, 0.0) // monoisotopic
  }
}
END_SECTION

START_SECTION ( void dia_massdiff_score(const std::vector< TransitionType > &transitions, SpectrumType spectrum, const std::vector< double > &normalized_library_intensity, double &ppm_score, double &ppm_score_weighted) )
{
  OpenSwath::SpectrumPtr sptr = prepareShiftedSpectrum();

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

START_SECTION ( bool DIAScoring::dia_ms1_massdiff_score(double precursor_mz, transitions, SpectrumType spectrum, double& ppm_score) )
{ 
  OpenSwath::SpectrumPtr sptr = prepareShiftedSpectrum();
  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.5, false, 30, 50, 4, 4); // here we use a large enough window so that none of our peaks falls out
  double ppm_score = 0;

  TEST_EQUAL(diascoring.dia_ms1_massdiff_score(500.0, sptr, ppm_score), true);
  TEST_REAL_SIMILAR(ppm_score, 15); // 15 ppm shifted

  TEST_EQUAL(diascoring.dia_ms1_massdiff_score(600.0, sptr, ppm_score), true);
  TEST_REAL_SIMILAR(ppm_score, 10); // 10 ppm shifted

  TEST_EQUAL(diascoring.dia_ms1_massdiff_score(100.0, sptr, ppm_score), false);
  TEST_REAL_SIMILAR(ppm_score, 0.5 * 1000000 / 100.0); // not present
}
END_SECTION

START_SECTION ( void dia_by_ion_score(SpectrumType spectrum, AASequence &sequence, int charge, double &bseries_score, double &yseries_score) )
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

  sptr->setMZArray(data1);
  sptr->setIntensityArray( data2 );

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use a large enough window so that none of our peaks falls out
  String sequence = "SYVAWDR";
  std::vector<double> bseries, yseries;
  AASequence a = AASequence::fromString(sequence);

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

START_SECTION((void set_dia_parameters(double dia_extract_window, double dia_centroided, double dia_byseries_intensity_min, double dia_byseries_ppm_diff, double dia_nr_isotopes, double dia_nr_charges)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void score_with_isotopes(SpectrumType spectrum, const std::vector< TransitionType > &transitions, double &dotprod, double &manhattan))
{
  OpenSwath::LightTransition mock_tr1;
  mock_tr1.product_mz = 500.;
  mock_tr1.product_charge = 1;
  mock_tr1.transition_name = "group1";
  mock_tr1.library_intensity = 5.;

  OpenSwath::LightTransition mock_tr2;
  mock_tr2.product_mz = 600.;
  mock_tr2.product_charge = 1;
  mock_tr2.transition_name = "group2";
  mock_tr2.library_intensity = 5.;

  OpenSwath::SpectrumPtr sptr = prepareSpectrum();

  std::vector<OpenSwath::LightTransition> transitions;
  transitions.push_back(mock_tr1);
  transitions.push_back(mock_tr2);

  DIAScoring diascoring;
  diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use a large enough window so that none of our peaks falls out
  double dotprod, manhattan;
  diascoring.score_with_isotopes(sptr,transitions,dotprod,manhattan);
  TEST_REAL_SIMILAR (dotprod, 0.729377020214466);
  TEST_REAL_SIMILAR (manhattan, 0.645047527406156);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

