// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/FLASHDeconvSpectrumFile.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <sstream>

///////////////////////////

using namespace OpenMS;
using namespace std;

/////////////////////////////////////////////////////////////
// Helper function to create a test PeakGroup
/////////////////////////////////////////////////////////////

PeakGroup createTestPeakGroup(double mono_mass, int min_charge, int max_charge, bool is_positive = true)
{
  PeakGroup pg(min_charge, max_charge, is_positive);
  pg.setMonoisotopicMass(mono_mass);
  pg.setIsotopeCosine(0.95f);
  pg.setChargeScore(0.9f);
  pg.setSNR(50.0f);
  pg.setQscore(0.85);
  pg.setQscore2D(0.80);
  pg.setRepAbsCharge((min_charge + max_charge) / 2);
  pg.setScanNumber(100);
  pg.setTargetDecoyType(PeakGroup::TargetDecoyType::target);

  // Add some LogMzPeaks
  for (int z = min_charge; z <= max_charge; ++z)
  {
    FLASHHelperClasses::LogMzPeak lp;
    lp.mz = mono_mass / z + (is_positive ? Constants::PROTON_MASS_U : -Constants::PROTON_MASS_U);
    lp.intensity = 1000.0f * z;
    lp.abs_charge = z;
    lp.is_positive = is_positive;
    lp.isotopeIndex = 0;
    lp.mass = mono_mass;
    pg.push_back(lp);
  }

  return pg;
}

/////////////////////////////////////////////////////////////
// Helper function to create a test DeconvolvedSpectrum
/////////////////////////////////////////////////////////////

DeconvolvedSpectrum createTestDeconvolvedSpectrum(int scan_number, uint ms_level, double rt)
{
  DeconvolvedSpectrum dspec(scan_number);

  // Create original spectrum
  MSSpectrum spec;
  spec.setRT(rt);
  spec.setMSLevel(ms_level);

  // Add some peaks to original spectrum
  Peak1D p1, p2, p3;
  p1.setMZ(500.0);
  p1.setIntensity(10000.0f);
  p2.setMZ(600.0);
  p2.setIntensity(15000.0f);
  p3.setMZ(700.0);
  p3.setIntensity(8000.0f);
  spec.push_back(p1);
  spec.push_back(p2);
  spec.push_back(p3);

  dspec.setOriginalSpectrum(spec);

  // Add at least 3 peak groups (required by topFD_min_peak_count_ = 3)
  dspec.push_back(createTestPeakGroup(10000.0, 5, 15));
  dspec.push_back(createTestPeakGroup(15000.0, 8, 20));
  dspec.push_back(createTestPeakGroup(20000.0, 10, 25));

  // For MS2, set precursor information
  if (ms_level > 1)
  {
    Precursor prec;
    prec.setMZ(800.0);
    prec.setIntensity(50000.0f);
    prec.setCharge(10);
    dspec.setPrecursor(prec);
    dspec.setPrecursorScanNumber(scan_number - 1);
    dspec.setActivationMethod(Precursor::ActivationMethod::HCD);

    // Set precursor peak group
    PeakGroup precPg = createTestPeakGroup(8000.0, 8, 12);
    precPg.setFeatureIndex(1);
    dspec.setPrecursorPeakGroup(precPg);
  }

  return dspec;
}

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_TEST(FLASHDeconvSpectrumFile, "$Id$")

/////////////////////////////////////////////////////////////
// Test writeDeconvolvedMassesHeader
/////////////////////////////////////////////////////////////

START_SECTION(static void writeDeconvolvedMassesHeader(std::ostream& os, uint ms_level, bool detail, bool report_decoy))
{
  // Test MS1 header without detail, without decoy
  {
    ostringstream oss;
    FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(oss, 1, false, false);
    String header = oss.str();

    TEST_EQUAL(header.hasSubstring("Index"), true)
    TEST_EQUAL(header.hasSubstring("FileName"), true)
    TEST_EQUAL(header.hasSubstring("ScanNum"), true)
    TEST_EQUAL(header.hasSubstring("RetentionTime"), true)
    TEST_EQUAL(header.hasSubstring("MonoisotopicMass"), true)
    TEST_EQUAL(header.hasSubstring("AverageMass"), true)
    TEST_EQUAL(header.hasSubstring("SumIntensity"), true)
    TEST_EQUAL(header.hasSubstring("IsotopeCosine"), true)
    TEST_EQUAL(header.hasSubstring("Qscore"), true)
    // Should NOT contain precursor info for MS1
    TEST_EQUAL(header.hasSubstring("PrecursorScanNum"), false)
    // Should NOT contain detail columns
    TEST_EQUAL(header.hasSubstring("PeakMZs"), false)
    // Should NOT contain decoy column
    TEST_EQUAL(header.hasSubstring("TargetDecoyType"), false)
  }

  // Test MS1 header with detail
  {
    ostringstream oss;
    FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(oss, 1, true, false);
    String header = oss.str();

    // Should contain detail columns
    TEST_EQUAL(header.hasSubstring("PeakMZs"), true)
    TEST_EQUAL(header.hasSubstring("PeakIntensities"), true)
    TEST_EQUAL(header.hasSubstring("PeakCharges"), true)
    TEST_EQUAL(header.hasSubstring("PerChargeIntensity"), true)
    TEST_EQUAL(header.hasSubstring("PerIsotopeIntensity"), true)
  }

  // Test MS1 header with decoy reporting
  {
    ostringstream oss;
    FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(oss, 1, false, true);
    String header = oss.str();

    TEST_EQUAL(header.hasSubstring("TargetDecoyType"), true)
    TEST_EQUAL(header.hasSubstring("Qvalue"), true)
  }

  // Test MS2 header (should contain precursor info)
  {
    ostringstream oss;
    FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(oss, 2, false, false);
    String header = oss.str();

    TEST_EQUAL(header.hasSubstring("PrecursorScanNum"), true)
    TEST_EQUAL(header.hasSubstring("PrecursorMz"), true)
    TEST_EQUAL(header.hasSubstring("PrecursorCharge"), true)
    TEST_EQUAL(header.hasSubstring("PrecursorMonoisotopicMass"), true)
  }

  // Test MS2 header with detail and decoy
  {
    ostringstream oss;
    FLASHDeconvSpectrumFile::writeDeconvolvedMassesHeader(oss, 2, true, true);
    String header = oss.str();

    TEST_EQUAL(header.hasSubstring("PrecursorScanNum"), true)
    TEST_EQUAL(header.hasSubstring("PeakMZs"), true)
    TEST_EQUAL(header.hasSubstring("TargetDecoyType"), true)
    TEST_EQUAL(header.hasSubstring("PrecursorQvalue"), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Test writeTopFDHeader
/////////////////////////////////////////////////////////////

START_SECTION(static void writeTopFDHeader(std::ostream& os, const Param& param))
{
  ostringstream oss;
  Param params;
  params.setValue("test_param", 42);
  params.setValue("another_param", "value");

  FLASHDeconvSpectrumFile::writeTopFDHeader(oss, params);
  String header = oss.str();

  TEST_EQUAL(header.hasSubstring("#FLASHDeconv generated msalign file"), true)
  TEST_EQUAL(header.hasSubstring("Parameters"), true)
  TEST_EQUAL(header.hasSubstring("test_param"), true)
  TEST_EQUAL(header.hasSubstring("42"), true)
  TEST_EQUAL(header.hasSubstring("another_param"), true)
  TEST_EQUAL(header.hasSubstring("value"), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Test writeIsobaricQuantification
/////////////////////////////////////////////////////////////

START_SECTION(static void writeIsobaricQuantification(std::ostream& os, std::vector<DeconvolvedSpectrum>& deconvolved_spectra))
{
  // Test with empty spectra - should produce header only (channel_count == 0)
  {
    ostringstream oss;
    std::vector<DeconvolvedSpectrum> empty_spectra;

    FLASHDeconvSpectrumFile::writeIsobaricQuantification(oss, empty_spectra);
    String output = oss.str();

    // Header is always written, even when no spectra have isobaric quantities
    TEST_EQUAL(output.hasPrefix("Scan\tPrecursorScan\tPrecursorMZ\t"), true)
  }

  // Test with MS2 spectra without isobaric quantities - should produce header only
  {
    ostringstream oss;
    std::vector<DeconvolvedSpectrum> spectra;
    spectra.push_back(createTestDeconvolvedSpectrum(100, 2, 600.0));
    spectra.push_back(createTestDeconvolvedSpectrum(102, 2, 610.0));

    FLASHDeconvSpectrumFile::writeIsobaricQuantification(oss, spectra);
    String output = oss.str();

    // Header is always written, even when channel_count == 0 (no isobaric quantities present)
    TEST_EQUAL(output.hasPrefix("Scan\tPrecursorScan\tPrecursorMZ\t"), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Test writeMzML
/////////////////////////////////////////////////////////////

START_SECTION(static void writeMzML(const MSExperiment& map, std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const String& deconvolved_mzML_file, const String& annotated_mzML_file, int mzml_charge, DoubleList tols))
{
  // Test with empty file names (should return early without error)
  {
    MSExperiment exp;
    std::vector<DeconvolvedSpectrum> empty_spectra;
    DoubleList tols = {10.0, 10.0};

    // Should not throw when both file names are empty
    FLASHDeconvSpectrumFile::writeMzML(exp, empty_spectra, "", "", 0, tols);
    TEST_EQUAL(true, true)  // If we get here, no exception was thrown
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Test writeTopFD with synthesized data
/////////////////////////////////////////////////////////////

START_SECTION(static void writeTopFD(DeconvolvedSpectrum& dspec, std::ostream& os, const String& filename, const double qval_threshold, const uint min_ms_level, bool randomize_precursor_mass, bool randomize_fragment_mass))
{
  // Test with MS1 spectrum
  {
    ostringstream oss;
    DeconvolvedSpectrum dspec = createTestDeconvolvedSpectrum(100, 1, 600.0);

    FLASHDeconvSpectrumFile::writeTopFD(dspec, oss, "test.mzML", 1.0, 1, false, false);
    String output = oss.str();

    // Check msalign format markers
    TEST_EQUAL(output.hasSubstring("BEGIN IONS"), true)
    TEST_EQUAL(output.hasSubstring("END IONS"), true)
    TEST_EQUAL(output.hasSubstring("FILE_NAME=test.mzML"), true)
    TEST_EQUAL(output.hasSubstring("SCANS=100"), true)
    TEST_EQUAL(output.hasSubstring("LEVEL=1"), true)
    // MS1 should not have precursor info
    TEST_EQUAL(output.hasSubstring("PRECURSOR_MZ"), false)
  }

  // Test with MS2 spectrum
  {
    ostringstream oss;
    DeconvolvedSpectrum dspec = createTestDeconvolvedSpectrum(101, 2, 605.0);

    FLASHDeconvSpectrumFile::writeTopFD(dspec, oss, "test.mzML", 1.0, 1, false, false);
    String output = oss.str();

    // Check msalign format for MS2
    TEST_EQUAL(output.hasSubstring("BEGIN IONS"), true)
    TEST_EQUAL(output.hasSubstring("END IONS"), true)
    TEST_EQUAL(output.hasSubstring("LEVEL=2"), true)
    // MS2 should have precursor info
    TEST_EQUAL(output.hasSubstring("PRECURSOR_MZ"), true)
    TEST_EQUAL(output.hasSubstring("PRECURSOR_CHARGE"), true)
    TEST_EQUAL(output.hasSubstring("PRECURSOR_MASS"), true)
    TEST_EQUAL(output.hasSubstring("MS_ONE_SCAN"), true)
    TEST_EQUAL(output.hasSubstring("ACTIVATION"), true)
  }

  // Test empty spectrum (should produce no output due to min peak count)
  {
    ostringstream oss;
    DeconvolvedSpectrum empty_dspec(50);
    MSSpectrum spec;
    spec.setRT(300.0);
    spec.setMSLevel(1);
    empty_dspec.setOriginalSpectrum(spec);

    FLASHDeconvSpectrumFile::writeTopFD(empty_dspec, oss, "test.mzML", 1.0, 1, false, false);
    // Empty spectrum should not produce output (below min peak count)
    TEST_EQUAL(oss.str().empty(), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
