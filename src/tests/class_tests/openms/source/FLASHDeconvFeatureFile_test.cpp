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

#include <OpenMS/FORMAT/FLASHDeconvFeatureFile.h>
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/PeakGroup.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <sstream>

///////////////////////////

START_TEST(FLASHDeconvFeatureFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

/////////////////////////////////////////////////////////////
// Helper function to create a test MassFeature
/////////////////////////////////////////////////////////////

FLASHHelperClasses::MassFeature createTestMassFeature(uint index, double mass, int ms_level = 1)
{
  FLASHHelperClasses::MassFeature mf;
  mf.index = index;
  mf.avg_mass = mass + 1.0;  // average mass slightly higher than mono
  mf.iso_offset = 0;
  mf.scan_number = 100;
  mf.min_scan_number = 95;
  mf.max_scan_number = 105;
  mf.rep_charge = 10;
  mf.min_charge = 5;
  mf.max_charge = 15;
  mf.charge_count = 11;
  mf.isotope_score = 0.95;
  mf.qscore = 0.85;
  mf.rep_mz = mass / 10.0 + 1.007276;  // approximate m/z
  mf.is_decoy = false;
  mf.ms_level = ms_level;

  // Initialize per charge/isotope intensities
  mf.per_charge_intensity.resize(20, 0.0f);
  mf.per_isotope_intensity.resize(10, 0.0f);
  for (int i = mf.min_charge; i <= mf.max_charge; ++i)
  {
    mf.per_charge_intensity[i] = 1000.0f * (i - mf.min_charge + 1);
  }
  for (size_t i = 0; i < 5; ++i)
  {
    mf.per_isotope_intensity[i] = 500.0f * (i + 1);
  }

  // Create a simple mass trace
  mf.mt = MassTrace();
  // Add some peaks to the mass trace
  Peak2D p1, p2, p3;
  p1.setRT(100.0);
  p1.setMZ(mass / 10.0);
  p1.setIntensity(10000.0f);
  p2.setRT(101.0);
  p2.setMZ(mass / 10.0);
  p2.setIntensity(15000.0f);
  p3.setRT(102.0);
  p3.setMZ(mass / 10.0);
  p3.setIntensity(8000.0f);

  std::vector<Peak2D> peaks = {p1, p2, p3};
  mf.mt = MassTrace(peaks);

  return mf;
}

/////////////////////////////////////////////////////////////
// Test writeHeader
/////////////////////////////////////////////////////////////

START_SECTION(static void writeHeader(std::ostream& os, bool report_decoy = false))
{
  // Test header without decoy reporting
  {
    ostringstream oss;
    FLASHDeconvFeatureFile::writeHeader(oss, false);
    String header = oss.str();

    TEST_EQUAL(header.hasSubstring("FeatureIndex"), true)
    TEST_EQUAL(header.hasSubstring("FileName"), true)
    TEST_EQUAL(header.hasSubstring("MSLevel"), true)
    TEST_EQUAL(header.hasSubstring("MonoisotopicMass"), true)
    TEST_EQUAL(header.hasSubstring("AverageMass"), true)
    TEST_EQUAL(header.hasSubstring("SumIntensity"), true)
    TEST_EQUAL(header.hasSubstring("PerChargeIntensity"), true)
    TEST_EQUAL(header.hasSubstring("PerIsotopeIntensity"), true)
    TEST_EQUAL(header.hasSubstring("IsDecoy"), false)  // Should NOT contain IsDecoy
    TEST_EQUAL(header.hasSubstring("\n"), true)  // Should end with newline
  }

  // Test header with decoy reporting
  {
    ostringstream oss;
    FLASHDeconvFeatureFile::writeHeader(oss, true);
    String header = oss.str();

    TEST_EQUAL(header.hasSubstring("IsDecoy"), true)  // Should contain IsDecoy
    TEST_EQUAL(header.hasSubstring("FeatureIndex"), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Test writeTopFDFeatureHeader
/////////////////////////////////////////////////////////////

START_SECTION(static void writeTopFDFeatureHeader(std::ostream& os, uint ms_level))
{
  // Test MS1 header
  {
    ostringstream oss;
    FLASHDeconvFeatureFile::writeTopFDFeatureHeader(oss, 1);
    String header = oss.str();

    TEST_EQUAL(header.hasSubstring("File_name"), true)
    TEST_EQUAL(header.hasSubstring("Feature_ID"), true)
    TEST_EQUAL(header.hasSubstring("Mass"), true)
    TEST_EQUAL(header.hasSubstring("Intensity"), true)
    TEST_EQUAL(header.hasSubstring("Min_time"), true)
    TEST_EQUAL(header.hasSubstring("Max_time"), true)
    TEST_EQUAL(header.hasSubstring("Apex_time"), true)
    TEST_EQUAL(header.hasSubstring("EC_score"), true)
    // MS1 header should NOT contain Spectrum_ID
    TEST_EQUAL(header.hasSubstring("Spectrum_ID"), false)
  }

  // Test MS2 header
  {
    ostringstream oss;
    FLASHDeconvFeatureFile::writeTopFDFeatureHeader(oss, 2);
    String header = oss.str();

    TEST_EQUAL(header.hasSubstring("File_name"), true)
    TEST_EQUAL(header.hasSubstring("Spectrum_ID"), true)  // MS2 has Spectrum_ID
    TEST_EQUAL(header.hasSubstring("Precursor_charge"), true)
    TEST_EQUAL(header.hasSubstring("Precursor_intensity"), true)
    TEST_EQUAL(header.hasSubstring("Fraction_feature_ID"), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Test writeFeatures with synthesized data
/////////////////////////////////////////////////////////////

START_SECTION(static void writeFeatures(const std::vector<FLASHHelperClasses::MassFeature>& mass_features, const String& file_name, std::ostream& os, bool report_decoy = false))
{
  // Test with empty features vector
  {
    ostringstream oss;
    std::vector<FLASHHelperClasses::MassFeature> empty_features;
    FLASHDeconvFeatureFile::writeFeatures(empty_features, "test.mzML", oss, false);
    TEST_EQUAL(oss.str().empty(), true)
  }

  // Test with synthesized features
  {
    ostringstream oss;
    std::vector<FLASHHelperClasses::MassFeature> features;
    features.push_back(createTestMassFeature(1, 10000.0));
    features.push_back(createTestMassFeature(2, 15000.0));

    FLASHDeconvFeatureFile::writeFeatures(features, "test_input.mzML", oss, false);
    String output = oss.str();

    // Check that output contains expected data
    TEST_EQUAL(output.hasSubstring("test_input.mzML"), true)
    TEST_EQUAL(output.hasSubstring("1\t"), true)  // Feature index 1
    TEST_EQUAL(output.hasSubstring("2\t"), true)  // Feature index 2

    // Should contain two lines (one per feature)
    Size line_count = std::count(output.begin(), output.end(), '\n');
    TEST_EQUAL(line_count, 2)
  }

  // Test with decoy reporting
  {
    ostringstream oss;
    std::vector<FLASHHelperClasses::MassFeature> features;
    auto mf = createTestMassFeature(1, 10000.0);
    mf.is_decoy = true;
    features.push_back(mf);

    FLASHDeconvFeatureFile::writeFeatures(features, "test.mzML", oss, true);
    String output = oss.str();

    // Should contain decoy indicator (1 for is_decoy=true)
    TEST_EQUAL(output.hasSubstring("\t1\t"), true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
// Test writeTopFDFeatures
/////////////////////////////////////////////////////////////

START_SECTION(static void writeTopFDFeatures(std::vector<DeconvolvedSpectrum>& deconvolved_spectra, const std::vector<FLASHHelperClasses::MassFeature>& mass_features, const std::map<int, double>& scan_rt_map, const String& file_name, std::ostream& os, uint ms_level))
{
  // Test with empty inputs
  {
    ostringstream oss;
    std::vector<DeconvolvedSpectrum> empty_spectra;
    std::vector<FLASHHelperClasses::MassFeature> empty_features;
    std::map<int, double> empty_rt_map;

    FLASHDeconvFeatureFile::writeTopFDFeatures(empty_spectra, empty_features, empty_rt_map, "test.mzML", oss, 1);
    TEST_EQUAL(oss.str().empty(), true)
  }

  // Test MS1 features with synthesized data
  {
    ostringstream oss;
    std::vector<DeconvolvedSpectrum> spectra;
    std::vector<FLASHHelperClasses::MassFeature> features;
    std::map<int, double> scan_rt_map;

    // Create test features
    features.push_back(createTestMassFeature(1, 10000.0, 1));
    features.push_back(createTestMassFeature(2, 15000.0, 1));

    // Create scan to RT mapping
    scan_rt_map[100] = 6000.0;  // 100 minutes in seconds

    FLASHDeconvFeatureFile::writeTopFDFeatures(spectra, features, scan_rt_map, "test_input.mzML", oss, 1);
    String output = oss.str();

    // Check output contains expected data for MS1
    TEST_EQUAL(output.hasSubstring("test_input.mzML"), true)
    // Should have feature entries
    TEST_EQUAL(output.empty(), false)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
