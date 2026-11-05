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
#include <OpenMS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>

#include <OpenMS/FORMAT/MzMLFile.h>

#include <algorithm>
#include <cmath>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

namespace
{
  // gives the test access to MultiplexFiltering::applyBlacklistEntries_(), to reserve peaks by hand, and to the number of owner
  // records it keeps
  class BlacklistOwnerTester : public MultiplexFiltering
  {
  public:
    using MultiplexFiltering::MultiplexFiltering;

    // the filtered peak peak_idx of pattern pattern_idx reserves peaks, given in order as (spectrum index, peak index, mass trace)
    void reserve(const std::vector<std::tuple<size_t, size_t, int>>& peaks, unsigned pattern_idx, size_t peak_idx)
    {
      std::vector<BlacklistEntry> entries;
      for (const auto& [rt_idx, mz_idx, mass_trace] : peaks)
      {
        entries.push_back({rt_idx, mz_idx, mass_trace});
      }
      applyBlacklistEntries_(entries, pattern_idx, peak_idx);
    }

    // one record per peak that more than one filtered peak reserved
    Size recordCount() const
    {
      return blacklist_owner_records_.size();
    }
  };
}

START_TEST(MultiplexFilteringCentroided, "$Id$")

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
std::vector<PeakPickerHiRes::PeakBoundary> boundaries;
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
shifts1.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0,"no_label"));
shifts1.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(8.0443702794,"Arg8"));
MultiplexDeltaMasses shifts2;
shifts2.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0,"no_label"));
MultiplexDeltaMasses::LabelSet label_set;
label_set.insert("Arg8");
label_set.insert("Arg8");
shifts2.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(2*8.0443702794,label_set));
std::vector<MultiplexIsotopicPeakPattern> patterns;
for (int c = charge_max; c >= charge_min; --c)
{
    MultiplexIsotopicPeakPattern pattern1(c, isotopes_per_peptide_max, shifts1, 0);
    patterns.push_back(pattern1);
    MultiplexIsotopicPeakPattern pattern2(c, isotopes_per_peptide_max, shifts2, 1);
    patterns.push_back(pattern2);
}

MultiplexFilteringCentroided* nullPointer = nullptr;
MultiplexFilteringCentroided* ptr;

START_SECTION(MultiplexFilteringCentroided(const MSExperiment& exp_picked, const std::vector<MultiplexIsotopicPeakPattern>& patterns, int isotopes_per_peptide_min, int isotopes_per_peptide_max, double intensity_cutoff, double rt_band, double mz_tolerance, bool mz_tolerance_unit, double peptide_similarity, double averagine_similarity, double averagine_similarity_scaling, std::string averagine_type="peptide"))
  MultiplexFilteringCentroided filtering(exp_picked, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
  ptr = new MultiplexFilteringCentroided(exp_picked, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
  TEST_NOT_EQUAL(ptr, nullPointer);
  delete ptr;
END_SECTION

MultiplexFilteringCentroided filtering(exp_picked, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);

START_SECTION(std::vector<MultiplexFilterResult> filter())
    std::vector<MultiplexFilteredMSExperiment> results = filtering.filter();
    TEST_EQUAL(results[0].size(), 0);
    TEST_EQUAL(results[1].size(), 0);
    TEST_EQUAL(results[2].size(), 0);
    TEST_EQUAL(results[3].size(), 0);
    TEST_EQUAL(results[4].size(), 4);
    TEST_EQUAL(results[5].size(), 4);
    TEST_EQUAL(results[6].size(), 4);
    TEST_EQUAL(results[7].size(), 0);
END_SECTION

START_SECTION([EXTRA] filter() gives the same result for any number of threads)
{
  // Dense data: in every spectrum the isotope patterns of neighbouring peptides interleave, so the peptides that
  // pass blacklist peaks that other peaks of the same spectrum need. The peaks of a spectrum are filtered in
  // parallel, and the result must still be that of filtering them one after another.
  MSExperiment dense;
  const double isotope_intensities[] = {1.0, 0.85, 0.45, 0.18, 0.06};
  for (int s = 0; s < 40; ++s)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setRT(100.0 + s);
    for (int k = 0; k < 400; ++k)
    {
      const int charge = 2 + k % 2;
      const double elution = std::exp(-0.5 * std::pow((s - k % 40) / 4.0, 2));
      if (elution < 0.01) continue;
      for (int i = 0; i < 5; ++i)
      {
        spectrum.push_back(Peak1D(400.0 + 0.61 * k + i * 1.003355 / charge, 1e6 * (1.0 + k % 7) * elution * isotope_intensities[i]));
      }
    }
    spectrum.sortByPosition();
    dense.addSpectrum(spectrum);
  }
  dense.updateRanges();

  MultiplexDeltaMasses unlabelled;
  unlabelled.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
  std::vector<MultiplexIsotopicPeakPattern> dense_patterns;
  for (int c = 3; c >= 2; --c)
  {
    dense_patterns.push_back(MultiplexIsotopicPeakPattern(c, isotopes_per_peptide_max, unlabelled, 0));
  }

  auto run = [&]([[maybe_unused]] int threads)
  {
#ifdef _OPENMP
    const int max_threads = omp_get_max_threads();
    omp_set_num_threads(threads);
#endif
    MultiplexFilteringCentroided dense_filtering(dense, dense_patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, 3.0, 10.0, true, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
    std::vector<MultiplexFilteredMSExperiment> results = dense_filtering.filter();
#ifdef _OPENMP
    omp_set_num_threads(max_threads);
#endif
    // every peak that passed, in the order it was accepted, with the number of its satellites
    std::vector<std::tuple<size_t, double, double, size_t>> peaks;
    for (size_t p = 0; p < results.size(); ++p)
    {
      for (size_t i = 0; i < results[p].size(); ++i)
      {
        peaks.emplace_back(p, results[p].getMZ(i), results[p].getRT(i), results[p].getPeak(i).size());
      }
    }
    // every blacklisted peak with its owners: the first and the last filtered peak (pattern, index in its result) that reserved
    // it, and the mass trace each one reserved it for
    std::vector<std::tuple<size_t, double, unsigned, size_t, int, unsigned, size_t, int>> owners;
    const MSExperiment blacklist = dense_filtering.getBlacklist();
    const std::vector<std::vector<MultiplexFiltering::BlacklistOwner>> blacklist_owners = dense_filtering.getBlacklistOwners();
    TEST_EQUAL(blacklist_owners.size(), blacklist.size())
    for (size_t s = 0; s < blacklist.size() && s < blacklist_owners.size(); ++s)
    {
      TEST_EQUAL(blacklist_owners[s].size(), blacklist[s].size())
      for (size_t k = 0; k < blacklist[s].size() && k < blacklist_owners[s].size(); ++k)
      {
        const MultiplexFiltering::BlacklistOwner& owner = blacklist_owners[s][k];
        owners.emplace_back(s, blacklist[s][k].getMZ(), owner.first.pattern_idx, owner.first.peak_idx, owner.first.mass_trace,
                            owner.last.pattern_idx, owner.last.peak_idx, owner.last.mass_trace);
        // the owners are filtered peaks that passed
        for (const MultiplexFiltering::BlacklistWriter& writer : {owner.first, owner.last})
        {
          TEST_EQUAL(writer.pattern_idx < results.size() && writer.peak_idx < results[writer.pattern_idx].size(), true)
        }
      }
    }
    return std::make_pair(peaks, owners);
  };

  const auto serial = run(1);
  TEST_EQUAL(serial.first.size() > 100, true)
  TEST_EQUAL(serial.second.size() > 100, true)
  // some peaks were reserved by two different filtered peaks
  TEST_EQUAL(std::count_if(serial.second.begin(), serial.second.end(), [](const auto& o)
                           { return std::make_pair(std::get<2>(o), std::get<3>(o)) != std::make_pair(std::get<5>(o), std::get<6>(o)); }) > 0, true)
  for (int threads : {2, 4, 8})
  {
    const auto parallel = run(threads);
    TEST_EQUAL(parallel.first == serial.first, true)
    TEST_EQUAL(parallel.second == serial.second, true)
  }
}
END_SECTION

START_SECTION([EXTRA] zeroth peak check: a satellite at the M-1 of a heavier peptide is found after an earlier pattern blacklisted peaks)
{
  // Three spectra with peptide A (3+, m/z 400, six peaks) and the doublet B (2+, light m/z 1500). Pattern 0 (3+ singlet)
  // finds A and blacklists its six peaks, so in the white spectra of pattern 1 (the 2+ doublet) every peak of B has an
  // index six lower than in the experiment. The zeroth peak check of B's heavy peptide finds the peak at its M-1 in the
  // white spectrum; whether it is a satellite has to be decided with its index in the experiment.
  const double c13 = Constants::C13C12_MASSDIFF_U;
  auto averagine = [](double mass)
  {
    std::vector<double> intensities;
    IsotopeDistribution distribution = CoarseIsotopePatternGenerator(6).estimateFromPeptideWeight(mass);
    for (Size i = 0; i < 6; ++i) { intensities.push_back(distribution[i].getIntensity()); }
    return intensities;
  };
  const std::vector<double> isotopes_b = averagine(2 * 1500.0);
  // heavy / light, the mass shift of the heavy peptide, whether a peak that is no satellite lies at the heavy M-1, the m/z of
  // A, the RT band, and the first spectrum that contains A
  auto filter = [&](double ratio, double shift, bool interferer, double mz_a = 400.0, double rt_band = 0.0, int first_a = 0)
  {
    const std::vector<double> isotopes_a = averagine(3 * mz_a);
    MSExperiment experiment;
    for (int s = 0; s < 3; ++s)
    {
      MSSpectrum spectrum;
      spectrum.setRT(10.0 + s);
      for (int i = 0; i < 6 && s >= first_a; ++i) { spectrum.push_back(Peak1D(mz_a + i * c13 / 3, float(1e6 * isotopes_a[i]))); }
      for (int i = 0; i < 6; ++i) { spectrum.push_back(Peak1D(1500.0 + i * c13 / 2, float(1e6 * isotopes_b[i]))); }
      for (int i = 0; i < 6; ++i) { spectrum.push_back(Peak1D(1500.0 + (shift + i * c13) / 2, float(1e6 * ratio * isotopes_b[i]))); }
      if (interferer) { spectrum.push_back(Peak1D(1500.0 + (shift - c13) / 2, float(0.8e6 * ratio * isotopes_b[0]))); }
      spectrum.sortByPosition();
      experiment.addSpectrum(spectrum);
    }
    experiment.updateRanges();
    MultiplexDeltaMasses singlet;
    singlet.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
    MultiplexDeltaMasses doublet;
    doublet.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
    doublet.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(shift, "heavy"));
    std::vector<MultiplexIsotopicPeakPattern> two_patterns{MultiplexIsotopicPeakPattern(3, 6, singlet, 0), MultiplexIsotopicPeakPattern(2, 6, doublet, 1)};
    MultiplexFilteringCentroided f(experiment, two_patterns, 3, 6, 1.0, rt_band, 10, true, 0.8, 0.75, 0.75, "peptide");
    f.setLogType(ProgressLogger::NONE);
    std::vector<MultiplexFilteredMSExperiment> results = f.filter();
    TEST_EQUAL(results[0].size(), 3 - first_a) // A in every spectrum that contains it
    return results[1].size();
  };
  // heavy shifted by six 13C: its M-1 is the light M+5, a satellite (mass trace 5). At heavy / light = 1 / 4 the light M+5
  // is 0.75 x the heavy M0, so above the 0.3 threshold. (Before: compared with index 5, B was missed in every spectrum.)
  TEST_EQUAL(filter(0.25, 6 * c13, false), 3)
  // At heavy / light = 1 the light M+5 is below the threshold either way.
  TEST_EQUAL(filter(1.0, 6 * c13, false), 3)
  // Lys8: the heavy M-1 is no light isotope position; a peak there of 0.8 x the heavy M0 rejects B. (Before: its white
  // index 6 was the index of B's light M0 in the experiment, a satellite, so B was taken.)
  TEST_EQUAL(filter(1.0, 8.0141988132, true), 0)
  // The same with A (m/z 1501.6) between B's light M0 and the heavy M-1: its six peaks lie between B's light M+3 and the
  // interferer, so the white index has to be mapped, not offset by the index of the primary peak (that would give B's light M+4).
  TEST_EQUAL(filter(1.0, 8.0141988132, true, 1501.6), 0)
  // Arg6 at 1 / 4 with an RT band of 3 s and A only in the second and third spectrum: in the second spectrum, the band begins
  // with the first, which has no blacklisted peak. The M-1 is mapped with the mapping of the spectrum it was found in.
  TEST_EQUAL(filter(0.25, 6 * c13, false, 400.0, 3.0, 1), 3)
}
END_SECTION

START_SECTION([EXTRA] zeroth peak check: a peak an earlier pattern has blacklisted does not reject a later pattern)
{
  // Three spectra with peptide A (3+, m/z 1000, six peaks) and peptide B (2+), whose M-1 is A's M+4. No peak of B is a peak
  // of A. A's M+4 is above 0.3 x B's M0, a significant peak at B's M-1 that is no satellite of B.
  const double c13 = Constants::C13C12_MASSDIFF_U;
  const double mz_a = 1000.0, mz_b = mz_a + 4 * c13 / 3 + c13 / 2;
  auto averagine = [](double mass)
  {
    std::vector<double> intensities;
    IsotopeDistribution distribution = CoarseIsotopePatternGenerator(6).estimateFromPeptideWeight(mass);
    for (Size i = 0; i < 6; ++i) { intensities.push_back(distribution[i].getIntensity()); }
    return intensities;
  };
  const std::vector<double> isotopes_a = averagine(3 * mz_a), isotopes_b = averagine(2 * mz_b);
  MSExperiment experiment;
  for (int s = 0; s < 3; ++s)
  {
    MSSpectrum spectrum;
    spectrum.setRT(10.0 + s);
    for (int i = 0; i < 6; ++i) { spectrum.push_back(Peak1D(mz_a + i * c13 / 3, float(1e6 * isotopes_a[i]))); }
    for (int i = 0; i < 6; ++i) { spectrum.push_back(Peak1D(mz_b + i * c13 / 2, float(1e5 * isotopes_b[i]))); }
    spectrum.sortByPosition();
    experiment.addSpectrum(spectrum);
  }
  experiment.updateRanges();
  MultiplexDeltaMasses singlet;
  singlet.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
  const MultiplexIsotopicPeakPattern pattern_a(3, 6, singlet, 0), pattern_b(2, 6, singlet, 1);
  // the number of filtered peaks of the last pattern at B's M0
  auto filter = [&](const std::vector<MultiplexIsotopicPeakPattern>& patterns)
  {
    MultiplexFilteringCentroided f(experiment, patterns, 3, 6, 1.0, 0.0, 10, true, 0.8, 0.75, 0.75, "peptide");
    f.setLogType(ProgressLogger::NONE);
    const std::vector<double> mzs = f.filter().back().getMZ();
    return std::count_if(mzs.begin(), mzs.end(), [&](double mz) { return std::fabs(mz - mz_b) < 1e-6; });
  };
  // B alone: A's M+4 rejects B in every spectrum.
  TEST_EQUAL(filter({pattern_b}), 0)
  // A first: its peaks are blacklisted, and the zeroth peak check of B looks in the white spectra, so B is found.
  TEST_EQUAL(filter({pattern_a, pattern_b}), 3)
}
END_SECTION

START_SECTION([EXTRA] getBlacklistOwners() keeps the first and the last filtered peak that reserved a peak)
{
  // two spectra of four peaks each, all above the intensity cutoff, so that the constructor sizes the blacklist and its owners to them
  MSExperiment two_spectra;
  for (int s = 0; s < 2; ++s)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(1);
    spectrum.setRT(100.0 + s);
    for (int k = 0; k < 4; ++k)
    {
      spectrum.push_back(Peak1D(500.0 + k, 1e5));
    }
    two_spectra.addSpectrum(spectrum);
  }
  two_spectra.updateRanges();
  BlacklistOwnerTester tester(two_spectra, patterns, isotopes_per_peptide_min, isotopes_per_peptide_max, intensity_cutoff, rt_band, mz_tolerance, mz_tolerance_unit, peptide_similarity, averagine_similarity, averagine_similarity_scaling, averagine_type);
  // a writer as "pattern index, index of the filtered peak in the result of the pattern, mass trace"
  auto writer = [](const MultiplexFiltering::BlacklistWriter& w)
  {
    return std::to_string(w.pattern_idx) + ", " + std::to_string(w.peak_idx) + ", " + std::to_string(w.mass_trace);
  };

  // A (filtered peak 0 of pattern 0) reserves peak 0 (mass trace 8) and peak 1 of spectrum 0, peak 1 twice (mass trace 1, then 3),
  // and peak 2 of spectrum 1 (mass trace 2)
  tester.reserve({{0, 0, 8}, {0, 1, 1}, {0, 1, 3}, {1, 2, 2}}, 0, 0);
  // A rewriting its own peak adds no record
  TEST_EQUAL(tester.recordCount(), 0)
  const std::vector<MultiplexFiltering::BlacklistOwner> owners_after_a = tester.getBlacklistOwners(0);
  // the entries of a list are checked only if the list has the right size, which the TEST_EQUAL before them checks; unlike
  // ABORT_IF, a wrong size does not skip the independent checks further down
  TEST_EQUAL(owners_after_a.size(), 2)
  if (owners_after_a.size() == 2)
  {
    // peak 1 of spectrum 0 has a single writer, A, which is its first and its last owner, with the mass trace it wrote last
    TEST_STRING_EQUAL(writer(owners_after_a[1].first), "0, 0, 3")
    TEST_STRING_EQUAL(writer(owners_after_a[1].last), "0, 0, 3")
  }

  // B (filtered peak 0 of pattern 1) reserves peaks 1 (mass trace 4) and 3 (mass trace 5) of spectrum 0; peak 2 stays white
  tester.reserve({{0, 1, 4}, {0, 3, 5}}, 1, 0);
  // B is the second writer of peak 1 of spectrum 0, which gets a record
  TEST_EQUAL(tester.recordCount(), 1)
  // C (filtered peak 1 of pattern 0) reserves peak 1 of spectrum 0 (mass trace 6) and peak 2 of spectrum 1 (mass trace 7)
  tester.reserve({{0, 1, 6}, {1, 2, 7}}, 0, 1);
  // C updates the record of peak 1 of spectrum 0, and is the second writer of peak 2 of spectrum 1, which gets a record
  TEST_EQUAL(tester.recordCount(), 2)

  const std::vector<std::vector<MultiplexFiltering::BlacklistOwner>> owners = tester.getBlacklistOwners();
  TEST_EQUAL(owners.size(), 2)
  if (owners.size() == 2)
  {
    // spectrum 0: peaks 0, 1 and 3 are reserved; the white peak 2 has no entry
    TEST_EQUAL(owners[0].size(), 3)
    if (owners[0].size() == 3)
    {
      // peak 0: A alone is the first and the last owner (mass trace 8, so that it differs from a value-initialized BlacklistWriter)
      TEST_STRING_EQUAL(writer(owners[0][0].first), "0, 0, 8")
      TEST_STRING_EQUAL(writer(owners[0][0].last), "0, 0, 8")
      // peak 1: the first owner is A with mass trace 3, which its own rewrite left (not 1, and not B's 4)
      TEST_STRING_EQUAL(writer(owners[0][1].first), "0, 0, 3")
      // peak 1: the last owner is C with mass trace 6; B, which reserved the peak in between, is not kept
      TEST_STRING_EQUAL(writer(owners[0][1].last), "0, 1, 6")
      // peak 3: B alone is the first and the last owner
      TEST_STRING_EQUAL(writer(owners[0][2].first), "1, 0, 5")
      TEST_STRING_EQUAL(writer(owners[0][2].last), "1, 0, 5")
    }
    // spectrum 1: only peak 2 is reserved, with a record of its own: A (mass trace 2) is the first, C (mass trace 7) the last owner
    TEST_EQUAL(owners[1].size(), 1)
    if (owners[1].size() == 1)
    {
      TEST_STRING_EQUAL(writer(owners[1][0].first), "0, 0, 2")
      TEST_STRING_EQUAL(writer(owners[1][0].last), "0, 1, 7")
    }
  }

  // the owners of spectrum 1 alone are the same
  const std::vector<MultiplexFiltering::BlacklistOwner> owners_1 = tester.getBlacklistOwners(1);
  TEST_EQUAL(owners_1.size(), 1)
  if (owners_1.size() == 1)
  {
    TEST_STRING_EQUAL(writer(owners_1[0].first), "0, 0, 2")
    TEST_STRING_EQUAL(writer(owners_1[0].last), "0, 1, 7")
  }
  // there are only two spectra
  TEST_EXCEPTION(Exception::IndexOverflow, tester.getBlacklistOwners(2))

  // entry k of the owners of a spectrum belongs to peak k of that spectrum of the blacklist, which skips the white peaks as well
  const MSExperiment blacklist = tester.getBlacklist();
  TEST_EQUAL(blacklist.size(), 2)
  if (blacklist.size() == 2)
  {
    TEST_EQUAL(blacklist[0].size(), 3)
    if (blacklist[0].size() == 3)
    {
      TEST_REAL_SIMILAR(blacklist[0][0].getMZ(), 500.0)
      TEST_REAL_SIMILAR(blacklist[0][1].getMZ(), 501.0)
      TEST_REAL_SIMILAR(blacklist[0][2].getMZ(), 503.0)
    }
    TEST_EQUAL(blacklist[1].size(), 1)
    if (blacklist[1].size() == 1)
    {
      TEST_REAL_SIMILAR(blacklist[1][0].getMZ(), 502.0)
    }
  }
}
END_SECTION

END_TEST
