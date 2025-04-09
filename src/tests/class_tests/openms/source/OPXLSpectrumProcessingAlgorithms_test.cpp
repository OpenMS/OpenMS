// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/XLMS/OPXLSpectrumProcessingAlgorithms.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>

using namespace OpenMS;

START_TEST(OPXLSpectrumProcessingAlgorithms, "$Id$")
TheoreticalSpectrumGeneratorXLMS specGen;
Param param = specGen.getParameters();
param.setValue("add_isotopes", "false");
param.setValue("add_metainfo", "true");
param.setValue("add_first_prefix_ion", "false");
param.setValue("y_intensity", 10.0, "Intensity of the y-ions");
param.setValue("add_a_ions", "false");
param.setValue("add_losses", "false");
param.setValue("add_precursor_peaks", "false");
param.setValue("add_k_linked_ions", "false");
specGen.setParameters(param);

PeakSpectrum theo_spec_1, theo_spec_2, exp_spec_1, exp_spec_2;
AASequence peptide = AASequence::fromString("PEPTIDE");
AASequence peptedi = AASequence::fromString("PEPTEDI");
specGen.getLinearIonSpectrum(exp_spec_1, peptide, 2, true, 3);
specGen.getLinearIonSpectrum(exp_spec_2, peptedi, 3, true, 3);

specGen.getLinearIonSpectrum(theo_spec_1, peptide, 3, true, 3);
specGen.getLinearIonSpectrum(theo_spec_2, peptedi, 4, true, 3);

START_SECTION(static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum))

  PeakSpectrum merged_spec = OPXLSpectrumProcessingAlgorithms::mergeAnnotatedSpectra(theo_spec_1, theo_spec_2);

  TEST_EQUAL(merged_spec.size(), 36)
  TEST_EQUAL(merged_spec.getIntegerDataArrays().size(), 1)
  TEST_EQUAL(merged_spec.getIntegerDataArrays()[0].size(), 36)
  TEST_EQUAL(merged_spec.getStringDataArrays()[0].size(), 36)
  TEST_EQUAL(merged_spec.getIntegerDataArrays()[0][10], 3)
  TEST_EQUAL(merged_spec.getStringDataArrays()[0][10], "[alpha|ci$y2]")
  TEST_EQUAL(merged_spec.getIntegerDataArrays()[0][20], 2)
  TEST_EQUAL(merged_spec.getStringDataArrays()[0][20], "[alpha|ci$y2]")
  TEST_REAL_SIMILAR(merged_spec[10].getMZ(), 83.04780)
  TEST_REAL_SIMILAR(merged_spec[20].getMZ(), 132.04732)

  for (Size i = 0; i < merged_spec.size()-1; ++i)
  {
    TEST_EQUAL(merged_spec[i].getMZ() <= merged_spec[i+1].getMZ(), true)
  }

END_SECTION

START_SECTION(static void getSpectrumAlignmentFastCharge(
      std::vector<std::pair<Size, Size> > & alignment, double fragment_mass_tolerance,
      bool fragment_mass_tolerance_unit_ppm,
      const PeakSpectrum& theo_spectrum,
      const PeakSpectrum& exp_spectrum,
      const DataArrays::IntegerDataArray& theo_charges,
      const DataArrays::IntegerDataArray& exp_charges,
      DataArrays::FloatDataArray& ppm_error_array,
      double intensity_cutoff = 0.0))

  std::vector <std::pair <Size, Size> > alignment1;
  std::vector <std::pair <Size, Size> > alignment2;

  theo_spec_1.sortByPosition();

  // slightly shift one of the exp spectra to get non-zero ppm error values
  PeakSpectrum exp_spec_3 = exp_spec_2;
  for (Peak1D p : exp_spec_3)
  {
    p.setMZ(p.getMZ() + 0.00001);
  }

  auto theo_2_it = getDataArrayByName(theo_spec_2.getIntegerDataArrays(), "charge");
  DataArrays::IntegerDataArray theo_2_charges = *theo_2_it;
  auto exp_3_it = getDataArrayByName(exp_spec_3.getIntegerDataArrays(), "charge");
  DataArrays::IntegerDataArray exp_3_charges = *exp_3_it;
  DataArrays::IntegerDataArray dummy_charges;

  DataArrays::FloatDataArray dummy_array;
  DataArrays::FloatDataArray ppm_error_array;
  OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment1, 50, true, theo_spec_1, exp_spec_1, dummy_charges, dummy_charges, dummy_array);
  OPXLSpectrumProcessingAlgorithms::getSpectrumAlignmentFastCharge(alignment2, 50, true, theo_spec_2, exp_spec_3, theo_2_charges, exp_3_charges, ppm_error_array);

  TEST_EQUAL(alignment1.size(), 15)
  TEST_EQUAL(alignment2.size(), 15)
  for (Size i = 0; i < alignment2.size(); i++)
  {
    TEST_REAL_SIMILAR(theo_spec_2[alignment2[i].first].getMZ(), exp_spec_3[alignment2[i].second].getMZ())
    TEST_REAL_SIMILAR((theo_spec_2[alignment2[i].first].getMZ() - exp_spec_3[alignment2[i].second].getMZ()) / theo_spec_2[alignment2[i].first].getMZ() / 1e-6, ppm_error_array[i])
  }

END_SECTION

END_TEST
