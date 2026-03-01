// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>
///////////////////////////

#include <OpenMS/FORMAT/MSPGenericFile.h>
#include <OpenMS/CONCEPT/Constants.h>

using namespace OpenMS;
using namespace std;

START_TEST(MetaboliteSpectralMatching, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MetaboliteSpectralMatching* ptr = nullptr;
MetaboliteSpectralMatching* null_ptr = nullptr;
START_SECTION(MetaboliteSpectralMatching())
{
	ptr = new MetaboliteSpectralMatching();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MetaboliteSpectralMatching())
{
	delete ptr;
}
END_SECTION

START_SECTION((virtual ~MetaboliteSpectralMatching()))
{
  // TODO
}
END_SECTION

START_SECTION((double computeHyperScore(MSSpectrum, MSSpectrum, const double &, const double &)))
{
  // TODO
}
END_SECTION

START_SECTION((void run(PeakMap &, MzTab &)))
{
  // TODO
}
END_SECTION

/////////////////////////////////////////////////////////////
// SpectralMatch getter/setter tests for CCS/IM
/////////////////////////////////////////////////////////////

START_SECTION(([SpectralMatch] double getObservedPrecursorDriftTime()))
{
  SpectralMatch sm;
  TEST_REAL_SIMILAR(sm.getObservedPrecursorDriftTime(), 0.0)
}
END_SECTION

START_SECTION(([SpectralMatch] void setObservedPrecursorDriftTime(const double&)))
{
  SpectralMatch sm;
  sm.setObservedPrecursorDriftTime(1.234);
  TEST_REAL_SIMILAR(sm.getObservedPrecursorDriftTime(), 1.234)
}
END_SECTION

START_SECTION(([SpectralMatch] double getFoundPrecursorCCS()))
{
  SpectralMatch sm;
  TEST_REAL_SIMILAR(sm.getFoundPrecursorCCS(), 0.0)
}
END_SECTION

START_SECTION(([SpectralMatch] void setFoundPrecursorCCS(const double&)))
{
  SpectralMatch sm;
  sm.setFoundPrecursorCCS(156.78);
  TEST_REAL_SIMILAR(sm.getFoundPrecursorCCS(), 156.78)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Integration test: MSP library with CCS parsed correctly
/////////////////////////////////////////////////////////////

START_SECTION(([Integration] MSPGenericFile parses CCS from MSP library))
{
  MSExperiment library;
  MSPGenericFile msp_parser;
  msp_parser.load(OPENMS_GET_TEST_DATA_PATH("MetaboliteSpectralMatching_test_library.msp"), library);

  // We should have 2 spectra
  TEST_EQUAL(library.size(), 2)

  // First spectrum: Caffeine with CCS = 141.2
  TEST_EQUAL(library[0].metaValueExists(Constants::UserParam::MSM_CCS), true)
  TEST_STRING_EQUAL(library[0].getMetaValue(Constants::UserParam::MSM_CCS).toString(), "141.2")

  // Second spectrum: L-Tryptophan with CCS = 156.78
  TEST_EQUAL(library[1].metaValueExists(Constants::UserParam::MSM_CCS), true)
  TEST_STRING_EQUAL(library[1].getMetaValue(Constants::UserParam::MSM_CCS).toString(), "156.78")

  // Verify other MSM meta values are also parsed
  TEST_EQUAL(library[0].metaValueExists(Constants::UserParam::MSM_SUM_FORMULA), true)
  TEST_STRING_EQUAL(library[0].getMetaValue(Constants::UserParam::MSM_SUM_FORMULA).toString(), "C8H10N4O2")
  TEST_EQUAL(library[0].metaValueExists(Constants::UserParam::MSM_PRECURSOR_ADDUCT), true)
  TEST_STRING_EQUAL(library[0].getMetaValue(Constants::UserParam::MSM_PRECURSOR_ADDUCT).toString(), "[M+H]+")
}
END_SECTION

/////////////////////////////////////////////////////////////
// Integration test: run() extracts CCS and drift time
/////////////////////////////////////////////////////////////

START_SECTION(([Integration] run() populates CCS and drift time in results))
{
  // Load the MSP library
  MSExperiment spec_db;
  MSPGenericFile msp_parser;
  msp_parser.load(OPENMS_GET_TEST_DATA_PATH("MetaboliteSpectralMatching_test_library.msp"), spec_db);

  // Create two synthetic experimental spectra (clustering requires >= 2)
  PeakMap msexp;

  // Spectrum 1: matching Caffeine (195.0877 m/z)
  {
    MSSpectrum exp_spectrum;
    exp_spectrum.setMSLevel(2);
    exp_spectrum.setRT(60.0);
    exp_spectrum.setDriftTime(5.5);

    Precursor prec;
    prec.setMZ(195.0877);
    prec.setCharge(1);
    exp_spectrum.setPrecursors({prec});

    Peak1D p1; p1.setMZ(42.0); p1.setIntensity(100);
    Peak1D p2; p2.setMZ(69.0); p2.setIntensity(200);
    Peak1D p3; p3.setMZ(110.0); p3.setIntensity(360);
    Peak1D p4; p4.setMZ(138.0); p4.setIntensity(800);
    Peak1D p5; p5.setMZ(195.0); p5.setIntensity(999);
    exp_spectrum.push_back(p1);
    exp_spectrum.push_back(p2);
    exp_spectrum.push_back(p3);
    exp_spectrum.push_back(p4);
    exp_spectrum.push_back(p5);

    msexp.addSpectrum(exp_spectrum);
  }

  // Spectrum 2: matching L-Tryptophan (205.0972 m/z)
  {
    MSSpectrum exp_spectrum;
    exp_spectrum.setMSLevel(2);
    exp_spectrum.setRT(120.0);
    exp_spectrum.setDriftTime(7.2);

    Precursor prec;
    prec.setMZ(205.0972);
    prec.setCharge(1);
    exp_spectrum.setPrecursors({prec});

    Peak1D p1; p1.setMZ(118.0); p1.setIntensity(300);
    Peak1D p2; p2.setMZ(146.0); p2.setIntensity(500);
    Peak1D p3; p3.setMZ(159.0); p3.setIntensity(700);
    Peak1D p4; p4.setMZ(188.0); p4.setIntensity(999);
    exp_spectrum.push_back(p1);
    exp_spectrum.push_back(p2);
    exp_spectrum.push_back(p3);
    exp_spectrum.push_back(p4);

    msexp.addSpectrum(exp_spectrum);
  }

  // Run the matching
  MzTab mztab_output;
  MetaboliteSpectralMatching msm;
  Param params;
  params.setValue("prec_mass_error_value", 10.0);
  params.setValue("frag_mass_error_value", 500.0);
  params.setValue("mass_error_unit", "ppm");
  params.setValue("report_mode", "top3");
  msm.setParameters(params);

  String dummy_db_name;
  msm.run(msexp, spec_db, mztab_output, dummy_db_name);

  // Check that results were produced
  const MzTabSmallMoleculeSectionRows& rows = mztab_output.getSmallMoleculeSectionRows();
  TEST_EQUAL(rows.empty(), false)

  // Verify results contain our drift time and CCS in optional columns
  if (!rows.empty())
  {
    bool found_dt = false;
    bool found_ccs = false;
    for (const auto& row : rows)
    {
      for (const auto& opt : row.opt_)
      {
        if (opt.first == "opt_observed_drift_time" && !found_dt)
        {
          found_dt = true;
        }
        if (opt.first == "opt_library_ccs" && !found_ccs)
        {
          found_ccs = true;
        }
      }
    }
    TEST_EQUAL(found_dt, true)
    TEST_EQUAL(found_ccs, true)
  }
}
END_SECTION




/////////////////////////////////////////////////////////////
// Test: CCS tolerance filtering rejects matches outside tolerance
/////////////////////////////////////////////////////////////

START_SECTION(([Integration] CCS tolerance filtering rejects matches outside tolerance))
{
  // Load the MSP library
  MSExperiment spec_db;
  MSPGenericFile msp_parser;
  msp_parser.load(OPENMS_GET_TEST_DATA_PATH("MetaboliteSpectralMatching_test_library.msp"), spec_db);

  // Create two experimental spectra with drift times far from library CCS
  PeakMap msexp;

  // Spectrum 1: matching Caffeine (195.0877 m/z), library CCS = 141.2
  // Set drift time to 5.5 (very far from 141.2 -> will be filtered out)
  {
    MSSpectrum exp_spectrum;
    exp_spectrum.setMSLevel(2);
    exp_spectrum.setRT(60.0);
    exp_spectrum.setDriftTime(5.5);

    Precursor prec;
    prec.setMZ(195.0877);
    prec.setCharge(1);
    exp_spectrum.setPrecursors({prec});

    Peak1D p1; p1.setMZ(42.0); p1.setIntensity(100);
    Peak1D p2; p2.setMZ(69.0); p2.setIntensity(200);
    Peak1D p3; p3.setMZ(110.0); p3.setIntensity(360);
    Peak1D p4; p4.setMZ(138.0); p4.setIntensity(800);
    Peak1D p5; p5.setMZ(195.0); p5.setIntensity(999);
    exp_spectrum.push_back(p1);
    exp_spectrum.push_back(p2);
    exp_spectrum.push_back(p3);
    exp_spectrum.push_back(p4);
    exp_spectrum.push_back(p5);

    msexp.addSpectrum(exp_spectrum);
  }

  // Spectrum 2: matching L-Tryptophan (205.0972 m/z), library CCS = 156.78
  {
    MSSpectrum exp_spectrum;
    exp_spectrum.setMSLevel(2);
    exp_spectrum.setRT(120.0);
    exp_spectrum.setDriftTime(7.2);

    Precursor prec;
    prec.setMZ(205.0972);
    prec.setCharge(1);
    exp_spectrum.setPrecursors({prec});

    Peak1D p1; p1.setMZ(118.0); p1.setIntensity(300);
    Peak1D p2; p2.setMZ(146.0); p2.setIntensity(500);
    Peak1D p3; p3.setMZ(159.0); p3.setIntensity(700);
    Peak1D p4; p4.setMZ(188.0); p4.setIntensity(999);
    exp_spectrum.push_back(p1);
    exp_spectrum.push_back(p2);
    exp_spectrum.push_back(p3);
    exp_spectrum.push_back(p4);

    msexp.addSpectrum(exp_spectrum);
  }

  // Run with tight CCS tolerance (1%) - should reject all since
  // drift time (5.5, 7.2) is very far from library CCS (141.2, 156.78)
  MzTab mztab_filtered;
  MetaboliteSpectralMatching msm_filtered;
  Param params_filtered;
  params_filtered.setValue("prec_mass_error_value", 10.0);
  params_filtered.setValue("frag_mass_error_value", 500.0);
  params_filtered.setValue("mass_error_unit", "ppm");
  params_filtered.setValue("report_mode", "top3");
  params_filtered.setValue("ccs_error_value", 1.0);  // 1% tolerance
  msm_filtered.setParameters(params_filtered);

  String dummy_db;
  msm_filtered.run(msexp, spec_db, mztab_filtered, dummy_db);

  // All matches should be filtered out because rel CCS error >> 1%
  const MzTabSmallMoleculeSectionRows& filtered_rows = mztab_filtered.getSmallMoleculeSectionRows();
  TEST_EQUAL(filtered_rows.empty(), true)

  // Now run WITHOUT CCS filtering (default 0 = disabled)
  MzTab mztab_unfiltered;
  MetaboliteSpectralMatching msm_unfiltered;
  Param params_unfiltered;
  params_unfiltered.setValue("prec_mass_error_value", 10.0);
  params_unfiltered.setValue("frag_mass_error_value", 500.0);
  params_unfiltered.setValue("mass_error_unit", "ppm");
  params_unfiltered.setValue("report_mode", "top3");
  // ccs_error_value defaults to 0.0 (disabled)
  msm_unfiltered.setParameters(params_unfiltered);

  String dummy_db2;
  msm_unfiltered.run(msexp, spec_db, mztab_unfiltered, dummy_db2);

  // Matches should NOT be filtered out when CCS filtering is disabled
  const MzTabSmallMoleculeSectionRows& unfiltered_rows = mztab_unfiltered.getSmallMoleculeSectionRows();
  TEST_EQUAL(unfiltered_rows.empty(), false)
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
