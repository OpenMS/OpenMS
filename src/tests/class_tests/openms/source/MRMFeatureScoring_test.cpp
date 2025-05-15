// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>

#include "OpenSwathTestHelper.h"

#include <OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>

#include <OpenMS/ANALYSIS/OPENSWATH/DIAScoring.h>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMScoring.h>
#include <OpenMS/OPENSWATHALGO/ALGO/Scoring.h>

///////////////////////////


///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MRMScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OpenSwath::MRMScoring* ptr = nullptr;
OpenSwath::MRMScoring* nullPointer = nullptr;

START_SECTION(MRMScoring())
{
  ptr = new OpenSwath::MRMScoring();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMScoring())
{
  delete ptr;
}
END_SECTION

///////////////////////////////////////////////////////////////////////////
// testing the individual scores that are produced
// calcXcorrCoelutionScore
// calcXcorrCoelutionWeightedScore
// calcXcorrShapeScore
// calcXcorrShapeWeightedScore
// calcLibraryScore
START_SECTION([EXTRA] test_scores())
{
  // load the mock objects
  MRMFeature mrmfeature = OpenSWATH_Test::createMockFeature();
  OpenSWATH_Test::MRMTransitionGroupType transition_group = OpenSWATH_Test::createMockTransitionGroup();

  // create the Interface objects
  OpenSwath::IMRMFeature * imrmfeature;
  imrmfeature = new MRMFeatureOpenMS(mrmfeature);

  //initialize the XCorr Matrix
  OpenSwath::MRMScoring mrmscore;
  std::vector<std::string> native_ids;
  for (Size i = 0; i < transition_group.getTransitions().size(); i++) {native_ids.push_back(transition_group.getTransitions()[i].getNativeID());}
  mrmscore.initializeXCorrMatrix(imrmfeature, native_ids);

  static const double arr_lib[] = {0.5,1,0.5};
  std::vector<double> normalized_library_intensity (arr_lib, arr_lib + sizeof(arr_lib) / sizeof(arr_lib[0]) );
  //mrmscore.standardize_data(normalized_library_intensity);
  double sumx = std::accumulate( normalized_library_intensity.begin(), normalized_library_intensity.end(), 0.0 );
  for(Size m =0; m<normalized_library_intensity.size();m++) { normalized_library_intensity[m] /= sumx;}

  TEST_REAL_SIMILAR(mrmscore.calcXcorrCoelutionScore(), 2.26491106406735)
  TEST_REAL_SIMILAR(mrmscore.calcXcorrCoelutionWeightedScore(normalized_library_intensity), 1.375)
  TEST_REAL_SIMILAR(mrmscore.calcXcorrShapeScore(), 0.757687954406132)
  TEST_REAL_SIMILAR(mrmscore.calcXcorrShapeWeightedScore(normalized_library_intensity), 0.7130856895)

  // numpy
  double library_corr, library_rmsd;
  double manhatten, dotproduct;
  double spectral_angle, rmsd;
  mrmscore.calcLibraryScore(imrmfeature, transition_group.getTransitions(), library_corr, library_rmsd, manhatten, dotproduct, spectral_angle, rmsd);
  TEST_REAL_SIMILAR(library_corr, -0.654591316)
  TEST_REAL_SIMILAR(library_rmsd, 0.5800337593)

  TEST_REAL_SIMILAR(manhatten, 1.279644714)
  TEST_REAL_SIMILAR(dotproduct, 0.34514801)

  TEST_REAL_SIMILAR(spectral_angle, 1.483262)
  TEST_REAL_SIMILAR(rmsd, 0.6727226674)

  delete imrmfeature;
}
END_SECTION

// testing the individual DIA (data independent / SWATH) scores that are produced
// dia_isotope_scores
// dia_massdiff_score
// dia_by_ion_score
START_SECTION((virtual void test_dia_scores()))
{
  OpenSWATH_Test::MRMTransitionGroupType transition_group;
  transition_group = OpenSWATH_Test::createMockTransitionGroup();

  PeakMap swath_map;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("ChromatogramExtractor_input.mzML"), swath_map);

  MRMFeature mrmfeature = OpenSWATH_Test::createMockFeature();

  int by_charge_state = 1;
  RangeMobility empty_im_range;

  // find spectrum that is closest to the apex of the peak (set to 3120) using binary search
  MSSpectrum OpenMSspectrum = (*swath_map.RTBegin( 3120 ));

  OpenSwath::BinaryDataArrayPtr intensity_array(new OpenSwath::BinaryDataArray);
  OpenSwath::BinaryDataArrayPtr mz_array(new OpenSwath::BinaryDataArray);
  for(MSSpectrum::iterator it = OpenMSspectrum.begin(); it != OpenMSspectrum.end(); it++)
  {
    mz_array->data.push_back(it->getMZ());
    intensity_array->data.push_back(it->getIntensity());
  }
  OpenSwath::SpectrumPtr sptr(new OpenSwath::Spectrum);
  sptr->setMZArray( mz_array );
  sptr->setIntensityArray( intensity_array);

  OpenSwath::MRMScoring mrmscore;
  DIAScoring diascoring;
  // diascoring.set_dia_parameters(0.05, false, 30, 50, 4, 4); // here we use 50 ppm and a cutoff of 30 in intensity -- because our peptide does not match with the testdata :-)
  Param p_dia = diascoring.getDefaults();
  p_dia.setValue("dia_extraction_window", 0.05);
  p_dia.setValue("dia_extraction_unit", "Th");
  p_dia.setValue("dia_centroided", "false");
  p_dia.setValue("dia_byseries_intensity_min", 30.0);
  p_dia.setValue("dia_byseries_ppm_diff", 50.0);
  p_dia.setValue("dia_nr_isotopes", 4);
  p_dia.setValue("dia_nr_charges", 4);
  diascoring.setParameters(p_dia);

  // calculate the normalized library intensity (expected value of the intensities)
  // Numpy
  // arr1 = [ 0,1,3,5,2,0 ];
  // arr2 = [ 1,3,5,2,0,0 ];
  // (arr1 - mean(arr1) ) / std(arr1)
  // (arr2 - mean(arr2) ) / std(arr2)
  static const double arr_lib[] = {1.0,0.5,0.5};
  std::vector<double> normalized_library_intensity (arr_lib, arr_lib + sizeof(arr_lib) / sizeof(arr_lib[0]) );
  //mrmscore.standardize_data(normalized_library_intensity);
  double sumx = std::accumulate( normalized_library_intensity.begin(), normalized_library_intensity.end(), 0.0 );
  for(Size m =0; m<normalized_library_intensity.size();m++) { normalized_library_intensity[m] /= sumx;}

  // Isotope correlation / overlap score: Is this peak part of an
  // isotopic pattern or is it the monoisotopic peak in an isotopic
  // pattern?
  OpenSwath::IMRMFeature * imrmfeature;
  imrmfeature = new MRMFeatureOpenMS(mrmfeature);
  // We have to reorder the transitions to make the tests work
  std::vector<OpenSWATH_Test::TransitionType> transitions = transition_group.getTransitions();
  double isotope_corr = 0, isotope_overlap = 0;

  std::vector<OpenSwath::SpectrumPtr> sptrArr;
  sptrArr.push_back(sptr);

  diascoring.dia_isotope_scores(transitions, sptrArr, imrmfeature, empty_im_range, isotope_corr, isotope_overlap);

  delete imrmfeature;

  // Mass deviation score
  double ppm_score = 0, ppm_score_weighted = 0;
  std::vector<double> ppm_errors;
  diascoring.dia_massdiff_score(transition_group.getTransitions(),
    sptrArr, normalized_library_intensity, empty_im_range, ppm_score, ppm_score_weighted, ppm_errors);

  // Presence of b/y series score
  double bseries_score = 0, yseries_score = 0;
  String sequence = "SYVAWDR";
  OpenMS::AASequence aas = AASequence::fromString(sequence);
  diascoring.dia_by_ion_score(sptrArr, aas, by_charge_state, empty_im_range, bseries_score, yseries_score);

  TEST_REAL_SIMILAR(isotope_corr, 0.2866618 * transition_group.getTransitions().size() )
  TEST_REAL_SIMILAR(isotope_corr, 0.85998565339479)
  TEST_REAL_SIMILAR(isotope_overlap, 0.0599970892071724)

  TEST_REAL_SIMILAR(ppm_score, 1.76388919944981 / 3)
  TEST_REAL_SIMILAR(ppm_score_weighted, 0.484116946070573)

  double ppm_expected[] = {0.17257858483247876, 0.79565530730866774, 0.79565530730866774};
  for (size_t i = 0; i < ppm_errors.size(); ++i)
  {
    TEST_REAL_SIMILAR(ppm_errors[i], ppm_expected[i]);
  }

  TEST_EQUAL(bseries_score, 0)
  TEST_EQUAL(yseries_score, 1)

  // b/y series score with modifications
  bseries_score = 0, yseries_score = 0;
  aas.setModification(1, "Phospho" ); // modify the Y
  diascoring.dia_by_ion_score(sptrArr, aas, by_charge_state, empty_im_range, bseries_score, yseries_score);
  TEST_EQUAL(bseries_score, 0)
  TEST_EQUAL(yseries_score, 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

