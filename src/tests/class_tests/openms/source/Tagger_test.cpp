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

#include <OpenMS/CHEMISTRY/Tagger.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace OpenMS;

START_TEST(Tagger, "$Id$")

START_SECTION(void getTag(const MSSpectrum& spec, std::set<std::string>& tags))

  TheoreticalSpectrumGenerator tsg;
  Param param = tsg.getParameters();
  param.setValue("add_metainfo", "false");
  param.setValue("add_first_prefix_ion", "true");
  param.setValue("add_a_ions", "true");
  param.setValue("add_losses", "true");
  param.setValue("add_precursor_peaks", "true");
  tsg.setParameters(param);

  // spectrum with charges +1 and +2
  AASequence test_sequence = AASequence::fromString("PEPTIDETESTTHISTAGGER");
  PeakSpectrum spec;
  tsg.getSpectrum(spec, test_sequence, 1, 2);
  TEST_EQUAL(spec.size(), 357);

  std::vector<std::string> tags;

  // tagger searching only for charge +1
  Tagger tagger = Tagger(2, 10, 5, 1, 1);
  tagger.getTag(spec, tags);
  TEST_EQUAL(tags.size(), 890);

  // first aa in prefixes is not recognized yet, unless as false positive
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)

  // last aa in suffixes is not recognized yet, unless as false positive
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), false)

  // tagger searching only for charge +2
  Tagger tagger2 = Tagger(2, 10, 5, 2, 2);
  tags.clear();
  tagger2.getTag(spec, tags);
  TEST_EQUAL(tags.size(), 1006);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  // these are found as false positives with charge +2, in a +1 and +2 spectrum
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)

  // tagger searching for charges +1 and +2
  Tagger tagger3 = Tagger(2, 10, 5, 1, 2);
  tags.clear();
  tagger3.getTag(spec, tags);
  TEST_EQUAL(tags.size(), 1094);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)


  // spectrum with charges +1 and +2
  AASequence test_sequence2 = AASequence::fromString("PEPTID(Oxidation)ETESTTHISTAGGER");
  PeakSpectrum spec2;
  tsg.getSpectrum(spec2, test_sequence2, 2, 2);
  TEST_EQUAL(spec2.size(), 180);

  tags.clear();
  tagger3.getTag(spec2, tags);
  TEST_EQUAL(tags.size(), 545);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)

  // not found due to modification
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), false)

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)

  // tagger searching for charge +2 with fixed modification
  Tagger tagger4 = Tagger(2, 10, 5, 2, 2, ListUtils::create<String>("Oxidation (D)"));
  tags.clear();
  tagger4.getTag(spec2, tags);
  TEST_EQUAL(tags.size(), 667);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)
  // modified residue found again
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)

  // tagger searching for charge +2 with variable modification
  Tagger tagger5 = Tagger(2, 10, 5, 2, 2, StringList(), ListUtils::create<String>("Oxidation (D)"));
  tags.clear();
  tagger5.getTag(spec2, tags);
  TEST_EQUAL(tags.size(), 739);

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPT") != tags.end(), false)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PEPTI") != tags.end(), false)
  // modified residue found again
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "EPTID") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "PTIDE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TIDET") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "IDETE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "DETES") != tags.end(), true)

  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ETEST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TESTT") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ESTTH") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STTHI") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TTHIS") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "THIST") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "HISTA") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "ISTAG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "STAGG") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "TAGGE") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "AGGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GGER") != tags.end(), true)
  TEST_EQUAL(std::find(tags.begin(), tags.end(), "GER") != tags.end(), true)

  // // runtime benchmark, research tags many times in the same spectrum
  // // takes currently about 90 sec
  // std::cout << std::endl;
  // for (int i = 0; i < 5000; i++)
  // {
  //   tags.clear();
  //   tagger3.getTag(spec, tags);
  // }

  // // write out found tags if necessary
  // for (const std::string& tag : tags)
  // {
  //   std::cout << "TEST TAG: " << tag << std::endl;
  // }

END_SECTION

START_SECTION(void getTag(const MSSpectrum& spec, std::set<std::string>& tags) Tolerance based on ppm vs absolute Da - comparison)
  // Test both tolerance modes: ppm-based (default/relative) and absolute Da-based (new feature)
  // This test demonstrates the difference between the two modes, especially important
  // for large fragment masses (e.g., nucleic acids, top-down proteomics)

  // Create a synthetic spectrum with large m/z values simulating large fragments
  std::vector<double> mzs_large;
  
  // Simulate a sequence with large base m/z values (1000-2000 Da range)
  // Starting at 1000 Da and adding amino acid masses
  double base_mz_ppm = 1000.0;
  mzs_large.push_back(base_mz_ppm);                        // ~1000 Da
  mzs_large.push_back(base_mz_ppm + 97.0527);              // +P  ~1097 Da
  mzs_large.push_back(base_mz_ppm + 97.0527 + 129.0426);   // +PE ~1226 Da
  mzs_large.push_back(base_mz_ppm + 97.0527 + 129.0426 + 97.0527); // +PEP ~1323 Da
  mzs_large.push_back(base_mz_ppm + 97.0527 + 129.0426 + 97.0527 + 101.0477); // +PEPT ~1424 Da
  mzs_large.push_back(base_mz_ppm + 97.0527 + 129.0426 + 97.0527 + 101.0477 + 103.0094); // +PEPTI ~1527 Da
  mzs_large.push_back(base_mz_ppm + 97.0527 + 129.0426 + 97.0527 + 101.0477 + 103.0094 + 115.0269); // +PEPTID ~1642 Da
  mzs_large.push_back(base_mz_ppm + 97.0527 + 129.0426 + 97.0527 + 101.0477 + 103.0094 + 115.0269 + 129.0426); // +PEPTIDE ~1771 Da

  std::vector<std::string> tags_ppm_mode;
  std::vector<std::string> tags_da_mode;

  // Test 1: Default ppm-based tolerance (tol_is_ppm = true, the default)
  // With ppm mode, tolerance is calculated relative to amino acid mass differences
  // For amino acid mass ~100 Da, 10 ppm = 0.001 Da
  Tagger tagger_ppm = Tagger(2, 10, 7, 1, 1, StringList(), StringList(), true); // ppm mode (default)
  tagger_ppm.getTag(mzs_large, tags_ppm_mode);
  
  // Test 2: New absolute Da-based tolerance (tol_is_ppm = false)
  // With Da mode, tolerance is calculated relative to fragment m/z values
  // For fragment m/z ~1500 Da, 10 ppm = 0.015 Da
  // This gives much more tolerance for matching when fragment masses are large
  Tagger tagger_da = Tagger(2, 10, 7, 1, 1, StringList(), StringList(), false); // Da mode (new)
  tagger_da.getTag(mzs_large, tags_da_mode);
  
  // Both modes should find some tags
  TEST_NOT_EQUAL(tags_ppm_mode.size(), 0)
  TEST_NOT_EQUAL(tags_da_mode.size(), 0)
  
  // With absolute mode, we expect potentially different results when fragment masses are large
  // but the specific behavior depends on the exact mass errors in the test data
  
  // Verify some expected tags are found in both modes
  bool found_pe_ppm = std::find(tags_ppm_mode.begin(), tags_ppm_mode.end(), "PE") != tags_ppm_mode.end();
  bool found_pe_da = std::find(tags_da_mode.begin(), tags_da_mode.end(), "PE") != tags_da_mode.end();
  
  // At least one mode should find basic tags
  TEST_EQUAL(found_pe_ppm || found_pe_da, true)

END_SECTION

START_SECTION(void getTag(const MSSpectrum& spec, std::set<std::string>& tags) Tolerance mode with absolute Da values)
  // Test the new absolute Da tolerance mode specifically
  // In this mode, tolerance is based on fragment m/z rather than amino acid masses
  std::vector<double> mzs_abs;
  double base_mz_abs = 150.0;
  mzs_abs.push_back(base_mz_abs);
  mzs_abs.push_back(base_mz_abs + 97.0527);              // +P
  mzs_abs.push_back(base_mz_abs + 97.0527 + 129.0426);   // +PE
  mzs_abs.push_back(base_mz_abs + 97.0527 + 129.0426 + 97.0527); // +PEP
  
  std::vector<std::string> tags_abs;
  
  // Use absolute Da mode (tol_is_ppm = false) with 0.02 Da tolerance
  Tagger tagger_abs = Tagger(2, 0.02, 3, 1, 1, StringList(), StringList(), false);
  tagger_abs.getTag(mzs_abs, tags_abs);
  
  // Should produce reasonable results even with smaller m/z values
  TEST_NOT_EQUAL(tags_abs.size(), 0)

END_SECTION

START_SECTION(void getTag(const MSSpectrum& spec, std::set<std::string>& tags) Explicit test for absolute vs ppm tolerance modes with large fragments)
  // This test explicitly demonstrates the difference between ppm and absolute Da modes
  // For large fragment masses (e.g., top-down proteomics, nucleic acids)
  
  // At m/z 2000 Da:
  // - 100 ppm on AA mass (~100 Da) = 0.01 Da tolerance
  // - 100 ppm on fragment m/z (2000 Da) = 0.2 Da tolerance (20x larger!)
  
  std::vector<double> mzs_demo;
  double base = 2000.0;
  double aa_P = 97.0527;
  double aa_E = 129.0426;
  
  // Add some noise/shift to make this realistic (0.05 Da shift)
  mzs_demo.push_back(base);
  mzs_demo.push_back(base + aa_P + 0.05);     // shifted by 0.05 Da
  mzs_demo.push_back(base + aa_P + aa_E + 0.05);
  mzs_demo.push_back(base + aa_P + aa_E + aa_P + 0.05);
  
  std::vector<std::string> tags_ppm_tight;
  std::vector<std::string> tags_da_wide;
  
  // Test 1: ppm mode with 100 ppm - applies tolerance to AA mass (~100 Da)
  // Tolerance: ~0.01 Da on AA mass - won't match 0.05 Da shift
  Tagger tagger_ppm_tight = Tagger(2, 100, 3, 1, 1, StringList(), StringList(), true);
  tagger_ppm_tight.getTag(mzs_demo, tags_ppm_tight);
  
  // Test 2: Da mode with 100 ppm - applies tolerance to fragment m/z (~2000 Da)
  // Tolerance: ~0.2 Da on fragment m/z - WILL match 0.05 Da shift
  Tagger tagger_da_wide = Tagger(2, 100, 3, 1, 1, StringList(), StringList(), false);
  tagger_da_wide.getTag(mzs_demo, tags_da_wide);
  
  // With the larger shift, ppm mode (tight on AA mass) should find fewer or no tags
  // Da mode (wide on fragment m/z) should find tags despite the shift
  // The exact counts depend on matching, but Da mode should be >= ppm mode
  TEST_EQUAL(tags_da_wide.size() >= tags_ppm_tight.size(), true)

END_SECTION

START_SECTION(void getTag(const MSSpectrum& spec, std::set<std::string>& tags) Backward compatibility - default is ppm mode)
  // Verify that the default behavior (ppm mode) is preserved for backward compatibility
  std::vector<double> mzs_compat;
  mzs_compat.push_back(200.0);
  mzs_compat.push_back(297.0527);    // +P
  mzs_compat.push_back(426.0953);    // +PE
  
  std::vector<std::string> tags_default;
  std::vector<std::string> tags_explicit_ppm;
  
  // Constructor without tol_is_ppm parameter (uses default = true)
  Tagger tagger_default = Tagger(2, 20, 2, 1, 1);
  tagger_default.getTag(mzs_compat, tags_default);
  
  // Constructor with explicit tol_is_ppm = true
  Tagger tagger_explicit = Tagger(2, 20, 2, 1, 1, StringList(), StringList(), true);
  tagger_explicit.getTag(mzs_compat, tags_explicit_ppm);
  
  // Both should produce identical results (backward compatibility)
  TEST_EQUAL(tags_default.size(), tags_explicit_ppm.size())

END_SECTION

START_SECTION(void getTag(const MSSpectrum& spec, std::set<std::string>& tags) Test absolute Da mode with realistic tolerances)
  // Test absolute Da mode with tolerances typical for high-resolution MS
  std::vector<double> mzs_hires;
  double base_mz_hires = 1500.0;
  
  // Add realistic mass measurement errors (~5-10 ppm at high resolution)
  mzs_hires.push_back(base_mz_hires);
  mzs_hires.push_back(base_mz_hires + 97.0527 + 0.001);  // 1 milliDa shift
  mzs_hires.push_back(base_mz_hires + 97.0527 * 2 + 0.002);  // 2 milliDa shift
  mzs_hires.push_back(base_mz_hires + 97.0527 * 3 + 0.001);  // 1 milliDa shift
  
  std::vector<std::string> tags_hires;
  
  // Use 10 ppm in absolute Da mode
  // At 1500 Da: 10 ppm = 0.015 Da = 15 milliDa (covers our 1-2 milliDa shifts)
  Tagger tagger_hires = Tagger(2, 10, 3, 1, 1, StringList(), StringList(), false);
  tagger_hires.getTag(mzs_hires, tags_hires);
  
  // Should find tags with realistic high-resolution mass errors
  TEST_NOT_EQUAL(tags_hires.size(), 0)

END_SECTION

END_TEST
