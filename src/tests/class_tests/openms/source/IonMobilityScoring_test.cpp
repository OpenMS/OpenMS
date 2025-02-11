// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
    static void driftScoring(const SpectrumSequence& spectrum,
                             const std::vector<TransitionType> & transitions,
                             OpenSwath_Scores & scores,
                             const double drift_target,
                             const RangeMobility im_range,
                             const double dia_extraction_window_,
                             const bool dia_extraction_ppm_,
                             const bool use_spline,
                             const double drift_extra) ))
{
  OpenSwath_Scores scores;

  double drift_target = 1.0;
  OpenMS::RangeMobility im_range_1(1.0);
  im_range_1.minSpanIfSingular(1.);
  double im_drift_extra_pcnt_ = 0.25;

  double dia_extract_window_ = 0.3;
  bool dia_extraction_ppm_ = false;

  // Test #1: Empty Spectrum
  OpenSwath::SpectrumPtr drift_spectrum(new OpenSwath::Spectrum());
  OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr ion_mobility(new OpenSwath::BinaryDataArray);

  drift_spectrum->setMZArray( mass);
  drift_spectrum->setIntensityArray( intensity);
  drift_spectrum->getDataArrays().push_back( ion_mobility );

  std::vector<OpenSwath::SpectrumPtr> sptrArr;
  sptrArr.push_back(drift_spectrum);

  IonMobilityScoring::driftScoring(sptrArr, transitions, scores,
                                   drift_target, im_range_1,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, 0);
  TEST_REAL_SIMILAR(scores.im_drift_weighted, 0);
  TEST_REAL_SIMILAR(scores.im_delta_score, 0);
  TEST_REAL_SIMILAR(scores.im_xcorr_shape_score, 0)
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 0)

  // Test #2: IM Scores (Condition 1/2)
  drift_spectrum = spec;
  std::vector<OpenSwath::SpectrumPtr> sptrArr2;
  sptrArr2.push_back(spec);

  // Test integrity of spectrum
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), 24)
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getIntensityArray()->data.size())
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getDriftTimeArray()->data.size())


  IonMobilityScoring::driftScoring(sptrArr2, transitions, scores,
                                   drift_target, im_range_1,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, (0.705405 + 0.4)/2.0 )
  TEST_REAL_SIMILAR(scores.im_drift_weighted, 0.662790697674419)
  TEST_REAL_SIMILAR(scores.im_delta_score, (0.294595 + 0.6)/2.0 )

  TEST_REAL_SIMILAR(scores.im_xcorr_shape_score, 0.892124778448826)
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 2.73205080756888)

  // Test #3: IM Scores (Condition 2/2)
  dia_extract_window_ = 0.1;
  IonMobilityScoring::driftScoring(sptrArr2, transitions, scores,
                                   drift_target, im_range_1,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, (0.5 + 0.4)/2.0 )
  TEST_REAL_SIMILAR(scores.im_drift_weighted, 0.489473684210526)
  TEST_REAL_SIMILAR(scores.im_delta_score, (0.5 + 0.6)/2.0 )

  TEST_REAL_SIMILAR(scores.im_xcorr_shape_score, 0.833513903989399)
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 0.910683602522959)

  // Test #4: deal with exactly one entry in mobilogram
  dia_extract_window_ = 0.3;
  drift_target = 1.05;
  OpenMS::RangeMobility im_range_2(drift_target);
  im_range_2.minSpanIfSingular(0.1);
  IonMobilityScoring::driftScoring(sptrArr2, transitions, scores,
                                   drift_target, im_range_2,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, 1.1)
  TEST_REAL_SIMILAR(scores.im_drift_weighted, 1.1)
  TEST_REAL_SIMILAR(scores.im_delta_score, 0.05)

  TEST_REAL_SIMILAR(scores.im_xcorr_shape_score, 0.0) // higher is better
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 1.) // lower is better

  // Test #5: deal with one zero transitions
  dia_extract_window_ = 0.3;
  OpenMS::RangeMobility im_range_3;
  im_range_3.setMin(1.0);
  im_range_3.setMax(1.3);
  drift_target = 1.1;
  IonMobilityScoring::driftScoring(sptrArr2, transitions, scores,
                                   drift_target, im_range_3,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, 1.16666666666667)
  TEST_REAL_SIMILAR(scores.im_drift_weighted, 1.16666666666667)
  TEST_REAL_SIMILAR(scores.im_delta_score, 0.0666666666666667)

  TEST_REAL_SIMILAR(scores.im_xcorr_shape_score, 1.0/3)
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 3.73205080756888)

  // Test #6: deal with all-zero transitions
  // IM range from 2.5 to 3.5
  OpenMS::RangeMobility im_range_4(3.0);
  im_range_4.minSpanIfSingular(1.);

  IonMobilityScoring::driftScoring(sptrArr2, transitions, scores,
                                   drift_target, im_range_4,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_drift, -1)
  TEST_REAL_SIMILAR(scores.im_drift_weighted, -1)
  TEST_REAL_SIMILAR(scores.im_delta_score, -1)

  TEST_EQUAL(std::isnan(scores.im_xcorr_shape_score), true)
  TEST_REAL_SIMILAR(scores.im_xcorr_coelution_score, 0)
}
END_SECTION

START_SECTION([EXTRA]
    static void driftScoringMS1(const SpectrumSequence& spectrum&,
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

  // IM range from 0.5 to 1.5
  double drift_target = 1.0;
  OpenMS::RangeMobility im_range(drift_target);
  im_range.minSpanIfSingular(1.);
  double im_drift_extra_pcnt_ = 0.25;

  double dia_extract_window_ = 0.3;
  bool dia_extraction_ppm_ = false;
  OpenSwath::SpectrumPtr drift_spectrum(new OpenSwath::Spectrum());
  OpenSwath::BinaryDataArrayPtr ion_mobility(new OpenSwath::BinaryDataArray);
  drift_spectrum->getDataArrays().push_back( ion_mobility );

  std::vector<OpenSwath::SpectrumPtr> sptrArr;
  sptrArr.push_back(drift_spectrum);

  IonMobilityScoring::driftScoringMS1(sptrArr, transitions, scores,
                                   drift_target, im_range,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
  drift_spectrum->setMZArray( mass);
  drift_spectrum->setIntensityArray( intensity);

  std::vector<OpenSwath::SpectrumPtr> sptrArr2;
  sptrArr2.push_back(drift_spectrum);

  IonMobilityScoring::driftScoringMS1(sptrArr2, transitions, scores,
                                   drift_target, im_range,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  drift_spectrum = ms1spec;

  std::vector<OpenSwath::SpectrumPtr> sptrArr3;
  sptrArr3.push_back(drift_spectrum);

  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), 8)
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getIntensityArray()->data.size())
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getDriftTimeArray()->data.size())

  IonMobilityScoring::driftScoringMS1(sptrArr3, transitions, scores,
                                   drift_target, im_range,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   false, im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_delta_score, 0.7)
}
END_SECTION

START_SECTION(([EXTRA]
    static void driftScoringMS1Contrast(std::vector<OpenSwath::SpectrumPtr> spectrum, OpenSwath::SpectrumPtr ms1spectrum,
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

  // IM from 0.5 to 1.5
  OpenMS::RangeMobility im_range_1(1);
  im_range_1.minSpanIfSingular(1.);
  double im_drift_extra_pcnt_ = 0.25;

  double dia_extract_window_ = 0.3;
  bool dia_extraction_ppm_ = false;

  OpenSwath::SpectrumPtr drift_spectrum(new OpenSwath::Spectrum());
  OpenSwath::SpectrumPtr drift_spectrum_ms1(new OpenSwath::Spectrum());
  OpenSwath::BinaryDataArrayPtr ion_mobility(new OpenSwath::BinaryDataArray);

  drift_spectrum_ms1->getDataArrays().push_back( ion_mobility );
  drift_spectrum->getDataArrays().push_back( ion_mobility );

  std::vector<OpenSwath::SpectrumPtr> sptrArr;
  std::vector<OpenSwath::SpectrumPtr> sptrArrMS1;

  sptrArr.push_back(drift_spectrum);
  sptrArrMS1.push_back(drift_spectrum_ms1);

  IonMobilityScoring::driftScoringMS1Contrast(sptrArr, sptrArrMS1, transitions, scores,
                                   im_range_1,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  OpenSwath::BinaryDataArrayPtr mass(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr intensity(new OpenSwath::BinaryDataArray);
  drift_spectrum->setMZArray( mass);
  drift_spectrum->setIntensityArray( intensity);
  drift_spectrum_ms1->setMZArray( mass);
  drift_spectrum_ms1->setIntensityArray( intensity);


  std::vector<OpenSwath::SpectrumPtr> sptrArr_2;
  std::vector<OpenSwath::SpectrumPtr> sptrArrMS1_2;

  sptrArr_2.push_back(drift_spectrum);
  sptrArrMS1_2.push_back(drift_spectrum_ms1);

  IonMobilityScoring::driftScoringMS1Contrast(sptrArr_2, sptrArrMS1, transitions, scores,
                                   im_range_1,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);


  std::vector<OpenSwath::SpectrumPtr> sptrArr_3;
  std::vector<OpenSwath::SpectrumPtr> sptrArrMS1_3;

  drift_spectrum = spec;
  drift_spectrum_ms1 = ms1spec;

  sptrArr_3.push_back(drift_spectrum);
  sptrArrMS1_3.push_back(drift_spectrum_ms1);

  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), 24)
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getIntensityArray()->data.size())
  TEST_EQUAL(drift_spectrum->getMZArray()->data.size(), drift_spectrum->getDriftTimeArray()->data.size())

  IonMobilityScoring::driftScoringMS1Contrast(sptrArr_3, sptrArrMS1_3, transitions, scores,
                                   im_range_1,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_contrast_coelution, 5.62132034355964)
  TEST_REAL_SIMILAR(scores.im_ms1_contrast_shape, 0.50991093654836)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_coelution, 2)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_shape, 0.56486260935015)

  dia_extract_window_ = 0.1;
  IonMobilityScoring::driftScoringMS1Contrast(sptrArr_3, sptrArrMS1_3, transitions, scores,
                                   im_range_1,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_contrast_coelution, 6)
  TEST_REAL_SIMILAR(scores.im_ms1_contrast_shape, 0)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_coelution, 6)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_shape, 0)

  // deal with exactly one entry in mobilogram
  dia_extract_window_ = 0.3;

  // IM Span from 1.0 to 1.1
  OpenMS::RangeMobility im_range_2(1.05);
  im_range_2.minSpanIfSingular(0.1);

  IonMobilityScoring::driftScoringMS1Contrast(sptrArr_3, sptrArrMS1_3, transitions, scores,
                                   im_range_2,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_contrast_coelution, 1)
  TEST_REAL_SIMILAR(scores.im_ms1_contrast_shape, 0)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_coelution, 1)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_shape, 0)

  // deal with one zero transitions
  dia_extract_window_ = 0.3;
  //Im Span from 1.0 to 1.3
  OpenMS::RangeMobility im_range_3(1.15);
  im_range_3.minSpanIfSingular(0.3);
  IonMobilityScoring::driftScoringMS1Contrast(sptrArr_3, sptrArrMS1_3, transitions, scores,
                                   im_range_3,
                                   dia_extract_window_, dia_extraction_ppm_,
                                   im_drift_extra_pcnt_);

  TEST_REAL_SIMILAR(scores.im_ms1_contrast_coelution, 3)
  TEST_REAL_SIMILAR(scores.im_ms1_contrast_shape, 0)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_coelution, 3)
  TEST_REAL_SIMILAR(scores.im_ms1_sum_contrast_shape, 0)

  // deal with all-zero transitions
  // IM span from 2.5 to 3.5
  OpenMS::RangeMobility im_range_4(3.);
  im_range_4.minSpanIfSingular(1.);

  IonMobilityScoring::driftScoringMS1Contrast(sptrArr_3, sptrArrMS1_3, transitions, scores,
                                   im_range_4,
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
