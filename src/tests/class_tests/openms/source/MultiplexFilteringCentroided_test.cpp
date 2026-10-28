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

#include <OpenMS/FEATUREFINDER/MultiplexFilteringCentroided.h>
#include <OpenMS/FEATUREFINDER/MultiplexFilteredMSExperiment.h>

#include <OpenMS/FORMAT/MzMLFile.h>

#include <cmath>
#include <tuple>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;

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
    return peaks;
  };

  const auto serial = run(1);
  TEST_EQUAL(serial.size() > 100, true)
  for (int threads : {2, 4, 8})
  {
    TEST_EQUAL(run(threads) == serial, true)
  }
}
END_SECTION

END_TEST
