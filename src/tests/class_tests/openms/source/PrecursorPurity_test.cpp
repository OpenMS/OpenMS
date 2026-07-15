// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/ID/PrecursorPurity.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>


START_TEST(PrecursorPurity, "$Id$")

using namespace OpenMS;
using namespace std;

PeakMap spectra;
MzMLFile f;
PeakFileOptions options;
options.clearMSLevels();
options.addMSLevel(1);
options.addMSLevel(2);
f.getOptions() = options;

// the file is a copy of OPENMS_GET_TEST_DATA_PATH("IsobaricChannelExtractor_6.mzML")
// which contains two MS1 spectra and 5 MS2 spectra between them
f.load(OPENMS_GET_TEST_DATA_PATH("PrecursorPurity_input.mzML"), spectra);

START_SECTION(static PurityScores computePrecursorPurity(const PeakSpectrum& ms1, const Precursor& pre, const double precursor_mass_tolerance, const bool precursor_mass_tolerance_unit_ppm))

  TEST_EQUAL(spectra.size(), 7)
  // the MS1 spectra are soectra[0] and spectra[6]

  Precursor pre = spectra[2].getPrecursors()[0];
  PrecursorPurity::PurityScores score = PrecursorPurity::computePrecursorPurity(spectra[6], pre, 10, true);
  TEST_REAL_SIMILAR(score.total_intensity, 11557777.1875)
  TEST_REAL_SIMILAR(score.target_intensity, 8923915)
  TEST_REAL_SIMILAR(score.signal_proportion, 0.77211)
  TEST_EQUAL(score.target_peak_count, 1)
  TEST_EQUAL(score.interfering_peak_count, 3)

  pre = spectra[3].getPrecursors()[0];
  score = PrecursorPurity::computePrecursorPurity(spectra[0], pre, 0.2, false);
  TEST_REAL_SIMILAR(score.total_intensity, 9098343.89062)
  TEST_REAL_SIMILAR(score.target_intensity, 7057944)
  TEST_REAL_SIMILAR(score.signal_proportion, 0.77573)
  TEST_EQUAL(score.target_peak_count, 1)
  TEST_EQUAL(score.interfering_peak_count, 4)

END_SECTION

START_SECTION(static computePrecursorPurities(const PeakMap& spectra, double precursor_mass_tolerance, bool precursor_mass_tolerance_unit_ppm))

  unordered_map<std::string, PrecursorPurity::PurityScores> purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);

  TEST_EQUAL(purityscores.size(), 5)

  // using the ID of an MS1 spectrum, a new ID for the map, adds a new Score to the map, initialized to 0
  TEST_REAL_SIMILAR(purityscores[spectra[0].getNativeID()].total_intensity, 0)
  TEST_REAL_SIMILAR(purityscores[spectra[0].getNativeID()].target_intensity, 0)
  TEST_REAL_SIMILAR(purityscores[spectra[0].getNativeID()].signal_proportion, 0)
  TEST_EQUAL(purityscores[spectra[0].getNativeID()].target_peak_count, 0)
  TEST_EQUAL(purityscores[spectra[0].getNativeID()].interfering_peak_count, 0)

  // 5 MS2 spectra between two MS1 spectra
  TEST_REAL_SIMILAR(purityscores[spectra[1].getNativeID()].total_intensity, 5517171)
  TEST_REAL_SIMILAR(purityscores[spectra[1].getNativeID()].target_intensity, 5517171)
  TEST_REAL_SIMILAR(purityscores[spectra[1].getNativeID()].signal_proportion, 1)
  TEST_EQUAL(purityscores[spectra[1].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[1].getNativeID()].interfering_peak_count, 0)

  TEST_REAL_SIMILAR(purityscores[spectra[2].getNativeID()].total_intensity, 11287967.625)
  TEST_REAL_SIMILAR(purityscores[spectra[2].getNativeID()].target_intensity, 7390478.5)
  TEST_REAL_SIMILAR(purityscores[spectra[2].getNativeID()].signal_proportion, 0.65472)
  TEST_EQUAL(purityscores[spectra[2].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[2].getNativeID()].interfering_peak_count, 3)

  TEST_REAL_SIMILAR(purityscores[spectra[3].getNativeID()].total_intensity, 9098343.89062)
  TEST_REAL_SIMILAR(purityscores[spectra[3].getNativeID()].target_intensity, 7057944)
  TEST_REAL_SIMILAR(purityscores[spectra[3].getNativeID()].signal_proportion, 0.77573)
  TEST_EQUAL(purityscores[spectra[3].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[3].getNativeID()].interfering_peak_count, 4)

  TEST_REAL_SIMILAR(purityscores[spectra[4].getNativeID()].total_intensity, 9762418.03906)
  TEST_REAL_SIMILAR(purityscores[spectra[4].getNativeID()].target_intensity, 7029896.5)
  TEST_REAL_SIMILAR(purityscores[spectra[4].getNativeID()].signal_proportion, 0.72009)
  TEST_EQUAL(purityscores[spectra[4].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[4].getNativeID()].interfering_peak_count, 5)

  TEST_REAL_SIMILAR(purityscores[spectra[5].getNativeID()].total_intensity, 5465177)
  TEST_REAL_SIMILAR(purityscores[spectra[5].getNativeID()].target_intensity, 5465177)
  TEST_REAL_SIMILAR(purityscores[spectra[5].getNativeID()].signal_proportion, 1)
  TEST_EQUAL(purityscores[spectra[5].getNativeID()].target_peak_count, 1)
  TEST_EQUAL(purityscores[spectra[5].getNativeID()].interfering_peak_count, 0)

  // using the ID of an MS1 spectrum, a new ID for the map, adds a new Score to the map, initialized to 0
  TEST_REAL_SIMILAR(purityscores[spectra[6].getNativeID()].total_intensity, 0)
  TEST_REAL_SIMILAR(purityscores[spectra[6].getNativeID()].target_intensity, 0)
  TEST_REAL_SIMILAR(purityscores[spectra[6].getNativeID()].signal_proportion, 0)
  TEST_EQUAL(purityscores[spectra[6].getNativeID()].target_peak_count, 0)
  TEST_EQUAL(purityscores[spectra[6].getNativeID()].interfering_peak_count, 0)

  TEST_REAL_SIMILAR(purityscores["randomString"].total_intensity, 0)
  TEST_REAL_SIMILAR(purityscores["randomString"].target_intensity, 0)
  TEST_REAL_SIMILAR(purityscores["randomString"].signal_proportion, 0)
  TEST_EQUAL(purityscores["randomString"].target_peak_count, 0)
  TEST_EQUAL(purityscores["randomString"].interfering_peak_count, 0)

  PeakMap spectra_copy = spectra;

  // If the PeakMap does not start with an MS1 spectrum, PrecursorPurity is not applicable and returns an empty vector
  spectra[0].setMSLevel(2);
  purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);
  TEST_EQUAL(purityscores.size(), 0)

  spectra[5].setMSLevel(1);
  purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);
  TEST_EQUAL(purityscores.size(), 0)

  spectra = spectra_copy;

  // duplicate IDs, skip the computation
  spectra = spectra_copy;
  TEST_EQUAL(spectra[2].getNativeID(), "controllerType=0 controllerNumber=1 scan=4")
  spectra[2].setNativeID(spectra[4].getNativeID());
  TEST_EQUAL(spectra[2].getNativeID(), "controllerType=0 controllerNumber=1 scan=8")
  TEST_EQUAL(spectra[4].getNativeID(), "controllerType=0 controllerNumber=1 scan=8")
  purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);
  TEST_EQUAL(purityscores.size(), 0)

  // empty ID, skip the computation
  spectra = spectra_copy;
  spectra[3].setNativeID("");
  purityscores = PrecursorPurity::computePrecursorPurities(spectra, 0.1, false);
  TEST_EQUAL(purityscores.size(), 0)

START_SECTION(static computeInterpolatedPrecursorPurity edge cases with fallback to early scan purity)

  // Test case 1: next_ms1_spec_idx is out of range (too large)
  std::vector<double> early_purity = PrecursorPurity::computeSingleScanPrecursorPurities(1, 0, spectra, 20.0);
  std::vector<double> interpolated_invalid_idx = PrecursorPurity::computeInterpolatedPrecursorPurity(1, 0, 999, spectra, 20.0);
  
  TEST_EQUAL(early_purity.size(), interpolated_invalid_idx.size())
  for (Size i = 0; i < early_purity.size(); ++i)
  {
    TEST_REAL_SIMILAR(early_purity[i], interpolated_invalid_idx[i])
  }

  // Test case 2: next_ms1_spec_idx is negative
  std::vector<double> interpolated_negative_idx = PrecursorPurity::computeInterpolatedPrecursorPurity(1, 0, -1, spectra, 20.0);
  
  TEST_EQUAL(early_purity.size(), interpolated_negative_idx.size())
  for (Size i = 0; i < early_purity.size(); ++i)
  {
    TEST_REAL_SIMILAR(early_purity[i], interpolated_negative_idx[i])
  }

  // Test case 3: next_ms1_spec_idx points to an MS2 spectrum (not MS1)
  std::vector<double> interpolated_ms2_idx = PrecursorPurity::computeInterpolatedPrecursorPurity(1, 0, 2, spectra, 20.0);
  
  TEST_EQUAL(early_purity.size(), interpolated_ms2_idx.size())
  for (Size i = 0; i < early_purity.size(); ++i)
  {
    TEST_REAL_SIMILAR(early_purity[i], interpolated_ms2_idx[i])
  }

  // Test case 4: Valid interpolation case for comparison
  std::vector<double> interpolated_valid = PrecursorPurity::computeInterpolatedPrecursorPurity(1, 0, 6, spectra, 20.0);
  
  // The valid interpolation should differ from early_purity (unless by coincidence)
  // We just ensure it executes without error and returns reasonable values
  TEST_EQUAL(early_purity.size(), interpolated_valid.size())
  for (Size i = 0; i < interpolated_valid.size(); ++i)
  {
    TEST_EQUAL(interpolated_valid[i] >= 0.0 && interpolated_valid[i] <= 1.0, true)
  }

END_SECTION

END_SECTION



START_SECTION((static Size countSPSPrecursorsMatchingPeptideFragments(const std::vector<Precursor>& sps_precursors, const AASequence& peptide, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, int max_fragment_charge)))

  AASequence peptide = AASequence::fromString("PEPTIDER");

  // some singly charged b/y fragment ion m/z of PEPTIDER (computed independently)
  const double b4 = 425.20308;
  const double y3 = 419.18849;
  const double y6 = 730.37299;

  auto make_precursor = [](double mz)
  {
    Precursor p;
    p.setMZ(mz);
    return p;
  };

  // empty inputs return 0
  TEST_EQUAL(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments({}, peptide, 0.02, false, 1), 0)
  TEST_EQUAL(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments({make_precursor(b4)}, AASequence(), 0.02, false, 1), 0)

  // three matching fragments, all should be counted (Da tolerance)
  std::vector<Precursor> sps = {make_precursor(b4), make_precursor(y3), make_precursor(y6)};
  TEST_EQUAL(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments(sps, peptide, 0.02, false, 1), 3)

  // add a peak that does not correspond to any fragment ion -> still 3 matches
  sps.push_back(make_precursor(500.0));
  TEST_EQUAL(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments(sps, peptide, 0.02, false, 1), 3)

  // ppm tolerance: 10 ppm around b4 (~0.004 Th) still matches, but a peak 0.05 Th off does not
  std::vector<Precursor> sps_ppm = {make_precursor(b4 + 0.003), make_precursor(b4 + 0.05)};
  TEST_EQUAL(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments(sps_ppm, peptide, 20.0, true, 1), 1)

  // a charge < 1 is treated as charge 1
  TEST_EQUAL(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments({make_precursor(b4)}, peptide, 0.02, false, 0), 1)

  // self-consistency check: precursors set to exact theoretical fragment m/z (charge 1 and 2) must all match
  TheoreticalSpectrumGenerator tsg;
  std::vector<float> theo_mzs;
  tsg.getPrefixAndSuffixIonsMZ(theo_mzs, peptide, 2);
  std::vector<Precursor> sps_theo;
  for (float mz : theo_mzs)
  {
    sps_theo.push_back(make_precursor(mz));
  }
  TEST_EQUAL(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments(sps_theo, peptide, 0.001, false, 2), sps_theo.size())

  // a doubly charged fragment is not matched when only singly charged fragments are considered
  const double y7_2plus = 430.21143;
  TEST_EQUAL(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments({make_precursor(y7_2plus)}, peptide, 0.02, false, 1), 0)
  TEST_EQUAL(PrecursorPurity::countSPSPrecursorsMatchingPeptideFragments({make_precursor(y7_2plus)}, peptide, 0.02, false, 2), 1)

END_SECTION

END_TEST
