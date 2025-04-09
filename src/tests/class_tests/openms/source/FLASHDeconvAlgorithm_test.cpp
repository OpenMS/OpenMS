// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim$
// $Authors: Jihyung Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvAlgorithm.h>
#include <OpenMS/FORMAT/MzMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FLASHDeconvAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FLASHDeconvAlgorithm* ptr = 0;
FLASHDeconvAlgorithm* null_ptr = 0;
START_SECTION(FLASHDeconvAlgorithm())
{
  ptr = new FLASHDeconvAlgorithm();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FLASHDeconvAlgorithm())
{
  delete ptr;
}
END_SECTION


/// < public methods without tests >
/// - default constructors and operators are not used (copy, move, assignment)
/// - setTargetMasses : only private member (which can not be accessed) is affected
/// - getDecoyDeconvolvedSpectrum, isDecoy, addPreviouslyDeconvolvedMonoMass, clearPreviouslyDeconvolvedMonoMasses: under development
/// - getAvgPPMError

FLASHDeconvAlgorithm fd_algo = FLASHDeconvAlgorithm();
Param fd_param;
fd_param.setValue("min_charge", 5);
fd_param.setValue("max_charge", 20);

START_SECTION((static int getNominalMass(const double mass)))
{
  double tmp_mass1 = 10000;
  double tmp_mass2 = 25000;

  TEST_EQUAL(fd_algo.getNominalMass(tmp_mass1), 9995);
  TEST_EQUAL(fd_algo.getNominalMass(tmp_mass2), 24987);
}
END_SECTION

START_SECTION((void calculateAveragine(const bool use_RNA_averagine)))
{
  fd_param.setValue("max_mass", 2000.);
  fd_algo.setParameters(fd_param);

  FLASHDeconvAlgorithm tmp_algo = FLASHDeconvAlgorithm();
  fd_param.setValue("max_mass", 100.);
  tmp_algo.setParameters(fd_param);

  fd_algo.calculateAveragine(false);
  tmp_algo.calculateAveragine(true);
  const auto &precalculated_avg = fd_algo.getAveragine();
  const auto &precalculated_avg_tmp = tmp_algo.getAveragine();

  TEST_EQUAL(precalculated_avg.getMaxIsotopeIndex(), 8);
  TEST_EQUAL(precalculated_avg.getApexIndex(50), 0);
  TOLERANCE_ABSOLUTE(0.1)
  TEST_REAL_SIMILAR(precalculated_avg.getAverageMassDelta(50), 0.0296591659229435);

  TEST_EQUAL(precalculated_avg_tmp.getMaxIsotopeIndex(), 3);
  TEST_EQUAL(precalculated_avg_tmp.getApexIndex(50), 0);
  TEST_REAL_SIMILAR(precalculated_avg_tmp.getAverageMassDelta(50), 0.025145817950033234);
}
END_SECTION

START_SECTION((PrecalculatedAveragine& getAveragine()))
{
  const auto &precalculated_avg = fd_algo.getAveragine();

  TEST_EQUAL(precalculated_avg.getMaxIsotopeIndex(), 8);
  TEST_EQUAL(precalculated_avg.getApexIndex(50), 0);
  TEST_REAL_SIMILAR(precalculated_avg.getAverageMassDelta(50), 0.0296591659229435);
}
END_SECTION

// load test data
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FLASHDeconv_sample_input1.mzML"), input);

// resetting fd_algo based on the test data
fd_param.setValue("max_mass", 50000.);
fd_algo.setParameters(fd_param);
fd_algo.calculateAveragine(false);

START_SECTION(DeconvolvedSpectrum& getDeconvolvedSpectrum())
{
  std::vector<DeconvolvedSpectrum> survey_specs;
  const std::map<int, std::vector<std::vector<float>>> null_map;

  fd_algo.performSpectrumDeconvolution(input[3], survey_specs, 4, null_map);

  DeconvolvedSpectrum d_ms1_spec = fd_algo.getDeconvolvedSpectrum();
  TEST_EQUAL(d_ms1_spec.size(), 1);
}
END_SECTION

START_SECTION((DeconvolvedSpectrum& performSpectrumDeconvolution(const MSSpectrum &spec, const std::vector< DeconvolvedSpectrum > &survey_scans, const int scan_number, const bool write_detail, const std::map< int, std::vector< std::vector< double >>> &precursor_map_for_FLASHIda)))
{
  std::vector<DeconvolvedSpectrum> survey_specs;
  const std::map<int, std::vector<std::vector<float>>> null_map;

  fd_algo.performSpectrumDeconvolution(input[3], survey_specs, 4, null_map);
  DeconvolvedSpectrum d_ms1_spec = fd_algo.getDeconvolvedSpectrum();
  survey_specs.push_back(d_ms1_spec);
  fd_algo.performSpectrumDeconvolution(input[5], survey_specs, 6, null_map);
  DeconvolvedSpectrum d_ms2_spec = fd_algo.getDeconvolvedSpectrum();
  TEST_EQUAL(d_ms1_spec.getScanNumber(), 4);
  TEST_EQUAL(d_ms1_spec.size(), 1);
  Precursor precursor = d_ms2_spec.getPrecursor();
  TOLERANCE_ABSOLUTE(1);
  TEST_EQUAL(d_ms1_spec.getPrecursorPeakGroup().size(), 0);
  TEST_EQUAL(d_ms2_spec.getPrecursorPeakGroup().size(), 0);
  TEST_EQUAL(precursor.getCharge(), 9);
  TOLERANCE_ABSOLUTE(100);
  TEST_REAL_SIMILAR(precursor.getIntensity(), 12293.4);
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

