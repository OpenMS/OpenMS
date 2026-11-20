// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/test_config.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteringProfile.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <cmath>

using namespace OpenMS;

START_TEST(MultiplexFilteringProfile, "$Id$")

// read data
MSExperiment exp;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("MultiplexFiltering.mzML"), exp);
exp.updateRanges();

// pick data
PeakPickerHiRes picker;
Param param = picker.getParameters();
param.setValue("ms_levels", ListUtils::create<Int>("1"));
param.setValue("signal_to_noise", 0.0);
picker.setParameters(param);
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_s;
std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > boundaries_exp_c;
MSExperiment exp_picked;
picker.pickExperiment(exp, exp_picked, boundaries_exp_s, boundaries_exp_c);

// set parameters
int charge_min = 1;
int charge_max = 4;
int isotopes_per_peptide_min = 3;
int isotopes_per_peptide_max = 6;
double intensity_cutoff = 10.0;
double rt_band = 3;
double mz_tolerance = 40;
bool mz_tolerance_unit = true;    // ppm (true), Da (false)
double peptide_similarity = 0.8;
double averagine_similarity = 0.75;
double averagine_similarity_scaling = 0.75;
std::string averagine_type="peptide";

// construct list of peak patterns
MultiplexDeltaMasses shifts1;
shifts1.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
shifts1.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(8.0443702794, "Arg8"));
MultiplexDeltaMasses shifts2;
shifts2.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
MultiplexDeltaMasses::LabelSet label_set;
label_set.insert("Arg8");
label_set.insert("Arg8");
shifts2.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(2*8.0443702794, label_set));
std::vector<MultiplexIsotopicPeakPattern> patterns;
for (int c = charge_max; c >= charge_min; --c)
{
    MultiplexIsotopicPeakPattern pattern1(c, isotopes_per_peptide_max, shifts1, 0);
    patterns.push_back(pattern1);
    MultiplexIsotopicPeakPattern pattern2(c, isotopes_per_peptide_max, shifts2, 1);
    patterns.push_back(pattern2);
}

MultiplexFilteringProfile* nullPointer = nullptr;
MultiplexFilteringProfile* ptr;

START_SECTION(MultiplexFilteringProfile(MSExperiment& exp_profile, const MSExperiment& exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >& boundaries, const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, std::string averagine_type="peptide"))
  MultiplexFilteringProfile filtering(exp, exp_picked, boundaries_exp_s, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
  ptr = new MultiplexFilteringProfile(exp, exp_picked, boundaries_exp_s, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
  TEST_NOT_EQUAL(ptr, nullPointer);
  delete ptr;
END_SECTION

MultiplexFilteringProfile filtering(exp, exp_picked, boundaries_exp_s, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);

START_SECTION(std::vector<MultiplexFilterResult> filter())
    std::vector<MultiplexFilteredMSExperiment> results = filtering.filter();
    TEST_EQUAL(results[0].size(), 0);
    TEST_EQUAL(results[1].size(), 0);
    TEST_EQUAL(results[2].size(), 0);
    TEST_EQUAL(results[3].size(), 0);
    TEST_EQUAL(results[4].size(), 4);
    TEST_EQUAL(results[5].size(), 5);
    TEST_EQUAL(results[6].size(), 4);
    TEST_EQUAL(results[7].size(), 0);
END_SECTION

START_SECTION([EXTRA] zeroth peak check: a satellite at the M-1 of a heavier peptide is found after an earlier pattern blacklisted peaks)
{
  // The data of the same section in MultiplexFilteringCentroided_test as profile spectra (a Gaussian of FWHM m/z / 60000 for
  // every peak): in profile mode, the zeroth peak check runs on the white centroided spectra too. Peptide A (3+, m/z 400, six
  // peaks) is found by pattern 0 (3+ singlet), which blacklists its six peaks, so in the white spectra of pattern 1 (2+ doublet
  // B, light m/z 1500) every peak of B has an index six lower than in the experiment.
  const double c13 = Constants::C13C12_MASSDIFF_U;
  auto averagine = [](double mass)
  {
    std::vector<double> intensities;
    IsotopeDistribution distribution = CoarseIsotopePatternGenerator(6).estimateFromPeptideWeight(mass);
    for (Size i = 0; i < 6; ++i) { intensities.push_back(distribution[i].getIntensity()); }
    return intensities;
  };
  const std::vector<double> isotopes_a = averagine(3 * 400.0), isotopes_b = averagine(2 * 1500.0);
  // heavy / light, the mass shift of the heavy peptide, and whether a peak that is no satellite lies at the heavy M-1
  auto filter = [&](double ratio, double shift, bool interferer)
  {
    std::vector<Peak1D> peaks;
    for (int i = 0; i < 6; ++i) { peaks.emplace_back(400.0 + i * c13 / 3, 1e6 * isotopes_a[i]); }
    for (int i = 0; i < 6; ++i) { peaks.emplace_back(1500.0 + i * c13 / 2, 1e6 * isotopes_b[i]); }
    for (int i = 0; i < 6; ++i) { peaks.emplace_back(1500.0 + (shift + i * c13) / 2, 1e6 * ratio * isotopes_b[i]); }
    if (interferer) { peaks.emplace_back(1500.0 + (shift - c13) / 2, 0.8e6 * ratio * isotopes_b[0]); }
    MSExperiment experiment;
    for (int s = 0; s < 3; ++s)
    {
      MSSpectrum spectrum;
      spectrum.setRT(10.0 + s);
      spectrum.setMSLevel(1);
      spectrum.setType(SpectrumSettings::SpectrumType::PROFILE);
      for (const Peak1D& peak : peaks)
      {
        // 21 data points from -5 to +5 sigma
        const double sigma = peak.getMZ() / 60000 / 2.3548;
        for (int k = -10; k <= 10; ++k) { spectrum.push_back(Peak1D(peak.getMZ() + k * sigma / 2, float(peak.getIntensity() * std::exp(-k * k / 8.0)))); }
      }
      spectrum.sortByPosition();
      experiment.addSpectrum(spectrum);
    }
    experiment.updateRanges();
    MSExperiment experiment_picked;
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary>> boundaries_spectra, boundaries_chromatograms;
    picker.pickExperiment(experiment, experiment_picked, boundaries_spectra, boundaries_chromatograms);
    MultiplexDeltaMasses singlet;
    singlet.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
    MultiplexDeltaMasses doublet;
    doublet.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
    doublet.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(shift, "heavy"));
    std::vector<MultiplexIsotopicPeakPattern> two_patterns{MultiplexIsotopicPeakPattern(3, 6, singlet, 0), MultiplexIsotopicPeakPattern(2, 6, doublet, 1)};
    MultiplexFilteringProfile f(experiment, experiment_picked, boundaries_spectra, two_patterns, 3, 6, 1.0, 0.0, 10, true, 0.8, 0.75, 0.75, "peptide");
    f.setLogType(ProgressLogger::NONE);
    std::vector<MultiplexFilteredMSExperiment> results = f.filter();
    TEST_EQUAL(results[0].size(), 3) // A in every spectrum
    return results[1].size();
  };
  // heavy shifted by six 13C: its M-1 is the light M+5, a satellite. At heavy / light = 1 / 4 it is above the 0.3 threshold.
  TEST_EQUAL(filter(0.25, 6 * c13, false), 3)
  // At heavy / light = 1 the light M+5 is below the threshold either way.
  TEST_EQUAL(filter(1.0, 6 * c13, false), 3)
  // Lys8: a peak of 0.8 x the heavy M0 at the heavy M-1, no light isotope position, rejects B.
  TEST_EQUAL(filter(1.0, 8.0141988132, true), 0)
}
END_SECTION

END_TEST
