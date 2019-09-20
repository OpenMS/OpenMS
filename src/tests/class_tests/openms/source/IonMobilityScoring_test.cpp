// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/IonMobilityScoring.h>
///////////////////////////

using namespace std;
using namespace OpenMS;
using namespace OpenSwath;

START_TEST(IonMobilityScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IonMobilityScoring* ptr = nullptr;
IonMobilityScoring* nullPointer = nullptr;

START_SECTION(IonMobilityScoring())
{
  ptr = new IonMobilityScoring();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~IonMobilityScoring())
{
  delete ptr;
}
END_SECTION

OpenSwath::LightTransition mock_tr1;
mock_tr1.product_mz = 500.2;
mock_tr1.precursor_mz = 700.2;
mock_tr1.fragment_charge = 1;
mock_tr1.transition_name = "group1";

OpenSwath::LightTransition mock_tr2;
mock_tr2.product_mz = 600.5;
mock_tr2.precursor_mz = 700.2;
mock_tr2.fragment_charge = 1;
mock_tr2.transition_name = "group2";

// create transitions, e.g. library intensity
std::vector<OpenSwath::LightTransition> transitions;
transitions.push_back(mock_tr1);
transitions.push_back(mock_tr2);

OpenSwath::SpectrumPtr spec(new OpenSwath::Spectrum());
{
  OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);

  mass->data.push_back(500.2);
  mass->data.push_back(500.2);
  mass->data.push_back(500.2);
  mass->data.push_back(500.2);
  mass->data.push_back(500.2);
  mass->data.push_back(500.2);
  mass->data.push_back(500.2);
  mass->data.push_back(500.2);

  mass->data.push_back(500.3);
  mass->data.push_back(500.3);
  mass->data.push_back(500.3);
  mass->data.push_back(500.3);
  mass->data.push_back(500.3);
  mass->data.push_back(500.3);
  mass->data.push_back(500.3);
  mass->data.push_back(500.3);

  mass->data.push_back(600.2);
  mass->data.push_back(600.3);
  mass->data.push_back(600.4);
  mass->data.push_back(600.5);
  mass->data.push_back(600.6);
  mass->data.push_back(600.7);
  mass->data.push_back(600.8);
  mass->data.push_back(600.9);

  OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
  intensity->data.push_back(10);
  intensity->data.push_back(20);
  intensity->data.push_back(30);
  intensity->data.push_back(40);
  intensity->data.push_back(40);
  intensity->data.push_back(30);
  intensity->data.push_back(20);
  intensity->data.push_back(10);

  intensity->data.push_back(10);
  intensity->data.push_back(20);
  intensity->data.push_back(30);
  intensity->data.push_back(40);
  intensity->data.push_back(40);
  intensity->data.push_back(30);
  intensity->data.push_back(20);
  intensity->data.push_back(10);

  intensity->data.push_back(20);
  intensity->data.push_back(20);
  intensity->data.push_back(20);
  intensity->data.push_back(20);
  intensity->data.push_back(20);
  intensity->data.push_back(20);
  intensity->data.push_back(20);
  intensity->data.push_back(20);

  OpenSwath::BinaryDataArrayPtr ion_mobility(new OpenSwath::BinaryDataArray);
  ion_mobility->data.push_back(0.1);
  ion_mobility->data.push_back(0.2);
  ion_mobility->data.push_back(0.3);
  ion_mobility->data.push_back(0.4); // 40
  ion_mobility->data.push_back(0.5); // 40
  ion_mobility->data.push_back(0.6);
  ion_mobility->data.push_back(0.7);
  ion_mobility->data.push_back(0.8);

  ion_mobility->data.push_back(0.5);
  ion_mobility->data.push_back(0.6);
  ion_mobility->data.push_back(0.7);
  ion_mobility->data.push_back(0.8); // 40
  ion_mobility->data.push_back(0.9); // 40
  ion_mobility->data.push_back(1.1);
  ion_mobility->data.push_back(1.2);
  ion_mobility->data.push_back(1.3);

  ion_mobility->data.push_back(0.1);
  ion_mobility->data.push_back(0.2);
  ion_mobility->data.push_back(0.3);
  ion_mobility->data.push_back(0.4);
  ion_mobility->data.push_back(0.5);
  ion_mobility->data.push_back(0.6);
  ion_mobility->data.push_back(0.7);
  ion_mobility->data.push_back(0.8);

  ion_mobility->description = "Ion Mobility";

  spec->setMZArray( mass);
  spec->setIntensityArray( intensity);
  spec->getDataArrays().push_back( ion_mobility );
}

OpenSwath::SpectrumPtr ms1spec(new OpenSwath::Spectrum());
{
  OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);

  mass->data.push_back(700.2);
  mass->data.push_back(700.3);
  mass->data.push_back(700.4);
  mass->data.push_back(700.5);
  mass->data.push_back(700.6);
  mass->data.push_back(700.7);
  mass->data.push_back(700.8);
  mass->data.push_back(700.9);

  OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
  intensity->data.push_back(10);
  intensity->data.push_back(20);
  intensity->data.push_back(30);
  intensity->data.push_back(40);
  intensity->data.push_back(40);
  intensity->data.push_back(30);
  intensity->data.push_back(20);
  intensity->data.push_back(10);

  OpenSwath::BinaryDataArrayPtr ion_mobility(new OpenSwath::BinaryDataArray);
  ion_mobility->data.push_back(0.2); // shifted by one
  ion_mobility->data.push_back(0.3);
  ion_mobility->data.push_back(0.4); // 40
  ion_mobility->data.push_back(0.5); // 40
  ion_mobility->data.push_back(0.6);
  ion_mobility->data.push_back(0.7);
  ion_mobility->data.push_back(0.8);
  ion_mobility->data.push_back(0.9);

  ion_mobility->description = "Ion Mobility";

  ms1spec->setMZArray( mass);
  ms1spec->setIntensityArray( intensity);
  ms1spec->getDataArrays().push_back( ion_mobility );
}

START_SECTION(([EXTRA]
    static void driftScoring(OpenSwath::SpectrumPtr spectrum,
                             const std::vector<TransitionType> & transitions,
                             OpenSwath_Scores & scores,
                             const double drift_lower,
                             const double drift_upper,
                             const double drift_target,
                             const double dia_extraction_window_,
                             const bool dia_extraction_ppm_,
                             const bool use_spline,
                             const double drift_extra) ))
{
  OpenSwath_Scores scores;

  double drift_lower = 0.5;
  double drift_upper = 1.5;
  double drift_target = 1.0;
  double im_drift_extra_pcnt_ = 0.25;

  double dia_extract_window_ = 0.3;
  bool dia_extraction_ppm_ = false;
  OpenSwath::SpectrumPtr drift_spectrum(new OpenSwath::Spectrum());
  OpenSwath::BinaryDataArrayPtr ion_mobility(new OpenSwath::BinaryDataArray);
  drift_spectrum->getDataArrays().push_back( ion_mobility );

  IonMobilityScoring::driftScoring(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
  drift_spectrum->setMZArray( mass);
  drift_spectrum->setIntensityArray( intensity);

  IonMobilityScoring::driftScoring(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  drift_spectrum = spec;

  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), 24)
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getIntensityArray()->data.size())
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getDriftTimeArray()->data.size())

  IonMobilityScoring::driftScoring(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, (0.705405 + 0.4)/2.0 )
  TEST_REAL_SIMILAR(scores.im_drift_weighted, 0.662790697674419)
  TEST_REAL_SIMILAR(scores.im_delta_score, (0.294595 + 0.6)/2.0 )

  TEST_REAL_SIMILAR(scores.im_xcorr_shape_score, 0.892124778448826)
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 2.73205080756888)

  dia_extract_window_ = 0.1;
  IonMobilityScoring::driftScoring(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, (0.5 + 0.4)/2.0 )
  TEST_REAL_SIMILAR(scores.im_drift_weighted, 0.489473684210526)
  TEST_REAL_SIMILAR(scores.im_delta_score, (0.5 + 0.6)/2.0 )

  TEST_REAL_SIMILAR(scores.im_xcorr_shape_score, 0.833513903989399)
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 0.910683602522959)

  // deal with exactly one entry in mobilogram
  dia_extract_window_ = 0.3;
  drift_lower = 1.0;
  drift_upper = 1.1;
  drift_target = 1.05;
  IonMobilityScoring::driftScoring(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, 1.1)
  TEST_REAL_SIMILAR(scores.im_drift_weighted, 1.1)
  TEST_REAL_SIMILAR(scores.im_delta_score, 0.05)

  TEST_REAL_SIMILAR(scores.im_xcorr_shape_score, 0.0) // higher is better
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 1.) // lower is better

  // deal with one zero transitions
  dia_extract_window_ = 0.3;
  drift_lower = 1.0;
  drift_upper = 1.3;
  drift_target = 1.1;
  IonMobilityScoring::driftScoring(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, 1.16666666666667)
  TEST_REAL_SIMILAR(scores.im_drift_weighted, 1.16666666666667)
  TEST_REAL_SIMILAR(scores.im_delta_score, 0.0666666666666667)

  TEST_REAL_SIMILAR(scores.im_xcorr_shape_score, 1.0/3)
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 3.73205080756888)

  // deal with all-zero transitions
  drift_lower = 2.5;
  drift_upper = 3.5;
  drift_target = 3.0;
  IonMobilityScoring::driftScoring(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, -1)
  TEST_REAL_SIMILAR(scores.im_drift_weighted, -1)
  TEST_EQUAL(std::isnan(scores.im_delta_score), true)

  TEST_EQUAL(std::isnan(scores.im_xcorr_shape_score), true)
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 0)
}
END_SECTION

START_SECTION([EXTRA]
    static void driftScoringMS1(OpenSwath::SpectrumPtr spectrum,
                                const std::vector<TransitionType> & transitions,
                                OpenSwath_Scores & scores,
                                const double drift_lower,
                                const double drift_upper,
                                const double drift_target,
                                const double dia_extract_window_,
                                const bool dia_extraction_ppm_,
                                const bool use_spline,
                                const double drift_extra))
{
  OpenSwath_Scores scores;

  double drift_lower = 0.5;
  double drift_upper = 1.5;
  double drift_target = 1.0;
  double im_drift_extra_pcnt_ = 0.25;

  double dia_extract_window_ = 0.3;
  bool dia_extraction_ppm_ = false;
  OpenSwath::SpectrumPtr drift_spectrum(new OpenSwath::Spectrum());
  OpenSwath::BinaryDataArrayPtr ion_mobility(new OpenSwath::BinaryDataArray);
  drift_spectrum->getDataArrays().push_back( ion_mobility );

  IonMobilityScoring::driftScoringMS1(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
  drift_spectrum->setMZArray( mass);
  drift_spectrum->setIntensityArray( intensity);

  IonMobilityScoring::driftScoringMS1(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  drift_spectrum = ms1spec;

  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), 8)
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getIntensityArray()->data.size())
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getDriftTimeArray()->data.size())

  IonMobilityScoring::driftScoringMS1(drift_spectrum, transitions, scores,
                                   drift_lower, drift_upper, drift_target,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_delta_score, 0.7)
}
END_SECTION

START_SECTION(([EXTRA]
    static void driftScoringMS1Contrast(OpenSwath::SpectrumPtr spectrum, OpenSwath::SpectrumPtr ms1spectrum, 
                             const std::vector<TransitionType> & transitions,
                             OpenSwath_Scores & scores,
                             const double drift_lower,
                             const double drift_upper,
                             const double drift_target,
                             const double dia_extraction_window_,
                             const bool dia_extraction_ppm_,
                             const bool use_spline,
                             const double drift_extra) ))
{
  OpenSwath_Scores scores;

  double drift_lower = 0.5;
  double drift_upper = 1.5;
  double im_drift_extra_pcnt_ = 0.25;

  double dia_extract_window_ = 0.3;
  bool dia_extraction_ppm_ = false;

  OpenSwath::SpectrumPtr drift_spectrum(new OpenSwath::Spectrum());
  OpenSwath::SpectrumPtr drift_spectrum_ms1(new OpenSwath::Spectrum());
  OpenSwath::BinaryDataArrayPtr ion_mobility(new OpenSwath::BinaryDataArray);

  drift_spectrum_ms1->getDataArrays().push_back( ion_mobility );
  drift_spectrum->getDataArrays().push_back( ion_mobility );

  IonMobilityScoring::driftScoringMS1Contrast(drift_spectrum, drift_spectrum_ms1, transitions, scores,
                                   drift_lower, drift_upper,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
  drift_spectrum->setMZArray( mass);
  drift_spectrum->setIntensityArray( intensity);
  drift_spectrum_ms1->setMZArray( mass);
  drift_spectrum_ms1->setIntensityArray( intensity);

  IonMobilityScoring::driftScoringMS1Contrast(drift_spectrum, drift_spectrum_ms1, transitions, scores,
                                   drift_lower, drift_upper,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  drift_spectrum = spec;
  drift_spectrum_ms1 = ms1spec;

  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), 24)
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getIntensityArray()->data.size())
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getDriftTimeArray()->data.size())

  IonMobilityScoring::driftScoringMS1Contrast(drift_spectrum, drift_spectrum_ms1, transitions, scores,
                                   drift_lower, drift_upper,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_contrast_coelution, 5.62132034355964)
  TEST_REAL_SIMILAR(scores.im_ms1_contrast_shape, 0.50991093654836)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_coelution, 2)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_shape, 0.56486260935015)

  dia_extract_window_ = 0.1;
  IonMobilityScoring::driftScoringMS1Contrast(drift_spectrum, drift_spectrum_ms1, transitions, scores,
                                   drift_lower, drift_upper,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_contrast_coelution, 6)
  TEST_REAL_SIMILAR(scores.im_ms1_contrast_shape, 0)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_coelution, 6)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_shape, 0)

  // deal with exactly one entry in mobilogram
  dia_extract_window_ = 0.3;
  drift_lower = 1.0;
  drift_upper = 1.1;
  IonMobilityScoring::driftScoringMS1Contrast(drift_spectrum, drift_spectrum_ms1, transitions, scores,
                                   drift_lower, drift_upper,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_contrast_coelution, 1)
  TEST_REAL_SIMILAR(scores.im_ms1_contrast_shape, 0)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_coelution, 1)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_shape, 0)

  // deal with one zero transitions
  dia_extract_window_ = 0.3;
  drift_lower = 1.0;
  drift_upper = 1.3;
  IonMobilityScoring::driftScoringMS1Contrast(drift_spectrum, drift_spectrum_ms1, transitions, scores,
                                   drift_lower, drift_upper,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_contrast_coelution, 3)
  TEST_REAL_SIMILAR(scores.im_ms1_contrast_shape, 0)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_coelution, 3)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_shape, 0)

  // deal with all-zero transitions
  drift_lower = 2.5;
  drift_upper = 3.5;
  IonMobilityScoring::driftScoringMS1Contrast(drift_spectrum, drift_spectrum_ms1, transitions, scores,
                                   drift_lower, drift_upper,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_contrast_coelution, 0)
  TEST_EQUAL(std::isnan(scores.im_ms1_contrast_shape), true)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_coelution, 0)
  TEST_EQUAL(std::isnan(scores.im_ms1_sum_contrast_shape), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

