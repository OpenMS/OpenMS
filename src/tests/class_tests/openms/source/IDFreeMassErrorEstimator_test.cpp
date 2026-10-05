// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/QC/IDFreeMassErrorEstimator.h>

#include <cmath>

namespace
{
  using namespace OpenMS;

  MSSpectrum makeSpectrum(Size cycle, Size target, bool profile = false, double fragment_shift_amplitude = 0.0015)
  {
    const double phase = static_cast<double>(cycle) * 0.73 + static_cast<double>(target) * 0.41;
    const double precursor_center = 500.0 + static_cast<double>(target) * 20.0;
    const double precursor_mz = precursor_center + std::sin(phase) * 0.0007;
    const double fragment_shift = std::sin(phase * 1.31) * fragment_shift_amplitude;

    MSSpectrum spectrum;
    spectrum.setMSLevel(2);
    spectrum.setRT(static_cast<double>(cycle * 10 + target));
    spectrum.setType(profile ? SpectrumSettings::SpectrumType::PROFILE : SpectrumSettings::SpectrumType::CENTROID);

    Precursor precursor;
    precursor.setMZ(precursor_mz);
    precursor.setCharge(2);
    spectrum.setPrecursors({precursor});

    for (Size i = 0; i < 12; ++i)
    {
      const double center = 150.0 + static_cast<double>(i) * 35.0 + fragment_shift;
      if (profile)
      {
        spectrum.emplace_back(center - 0.01, 25.0f);
        spectrum.emplace_back(center, 100.0f + static_cast<float>(i));
        spectrum.emplace_back(center + 0.01, 25.0f);
      }
      else
      {
        spectrum.emplace_back(center, 100.0f + static_cast<float>(i));
      }
    }
    spectrum.sortByPosition();
    return spectrum;
  }
}

START_TEST(IDFreeMassErrorEstimator, "$Id$")

using namespace OpenMS;

START_SECTION((IDFreeMassErrorEstimator()))
{
  IDFreeMassErrorEstimator estimator;
  TEST_EQUAL(estimator.getParameters().top_peaks, 50)
}
END_SECTION

START_SECTION((Result compute(const MSExperiment& experiment)))
{
  IDFreeMassErrorEstimator::Parameters parameters;
  parameters.min_spectrum_pairs = 4;
  parameters.min_precursor_clusters = 2;
  parameters.min_tolerance_pairs = 4;
  parameters.min_tolerance_clusters = 2;
  parameters.min_fragment_pairs = 10;
  parameters.min_fragment_tolerance_pairs = 10;
  parameters.min_fragment_tolerance_spectra = 2;

  MSExperiment experiment;
  for (Size cycle = 0; cycle < 8; ++cycle)
  {
    for (Size target = 0; target < 4; ++target)
    {
      experiment.addSpectrum(makeSpectrum(cycle, target));
    }
  }

  IDFreeMassErrorEstimator estimator(parameters);
  const auto result = estimator.compute(experiment);

  TEST_EQUAL(result.precursor_ppm.has_value(), true)
  TEST_EQUAL(result.precursor_da.has_value(), true)
  TEST_EQUAL(result.precursor_tolerance_ppm.has_value(), true)
  TEST_EQUAL(result.fragment_ppm.has_value(), true)
  TEST_EQUAL(result.fragment_tolerance_ppm.has_value(), true)
  TEST_EQUAL(result.fragment_tolerance_da.has_value(), false)
  TEST_EQUAL(result.fragment_resolution_regime, IDFreeMassErrorEstimator::FragmentResolutionRegime::HIGH_RESOLUTION)
  TEST_EQUAL(result.diagnostics.precursor_clusters_used, 4)
  TEST_EQUAL(result.diagnostics.precursor_eligible_ms2, 32)
  TEST_EQUAL(result.diagnostics.fragment_centroid_spectra, 32)
  TEST_EQUAL(result.precursor_tolerance_ppm->unit, "ppm")
  TEST_EQUAL(result.fragment_tolerance_ppm->unit, "ppm")
  TEST_EQUAL(result.precursor_tolerance_ppm->tolerance > 0.0, true)
  TEST_EQUAL(result.fragment_tolerance_ppm->tolerance > 0.0, true)
}
END_SECTION


START_SECTION((low-resolution fragments select an uncensored Da tolerance))
{
  IDFreeMassErrorEstimator::Parameters parameters;
  parameters.min_spectrum_pairs = 4;
  parameters.min_precursor_clusters = 2;
  parameters.min_tolerance_pairs = 4;
  parameters.min_tolerance_clusters = 2;
  parameters.min_fragment_pairs = 10;
  parameters.min_fragment_tolerance_pairs = 10;
  parameters.min_fragment_tolerance_spectra = 2;

  IDFreeMassErrorEstimator estimator(parameters);
  for (Size cycle = 0; cycle < 8; ++cycle)
  {
    for (Size target = 0; target < 4; ++target)
    {
      estimator.consumeSpectrum(makeSpectrum(cycle, target, false, 0.08));
    }
  }

  const auto result = estimator.getResult();
  TEST_EQUAL(result.fragment_resolution_regime, IDFreeMassErrorEstimator::FragmentResolutionRegime::LOW_RESOLUTION)
  TEST_REAL_SIMILAR(result.fragment_match_window_da, 0.5)
  TEST_EQUAL(result.fragment_window_censored, false)
  TEST_EQUAL(result.fragment_tolerance_ppm.has_value(), false)
  TEST_EQUAL(result.fragment_tolerance_da.has_value(), true)
  TEST_REAL_SIMILAR(result.fragment_tolerance_da->tolerance, 0.4654778152135567)
}
END_SECTION

START_SECTION((void consumeSpectrum(const MSSpectrum& spectrum)))
{
  IDFreeMassErrorEstimator::Parameters parameters;
  parameters.min_spectrum_pairs = 2;
  parameters.min_precursor_clusters = 1;
  parameters.min_tolerance_pairs = 2;
  parameters.min_tolerance_clusters = 1;
  parameters.min_fragment_pairs = 4;
  parameters.min_fragment_tolerance_pairs = 4;
  parameters.min_fragment_tolerance_spectra = 1;

  IDFreeMassErrorEstimator estimator(parameters);
  for (Size cycle = 0; cycle < 5; ++cycle)
  {
    estimator.consumeSpectrum(makeSpectrum(cycle, 0, true));
  }
  const auto result = estimator.getResult();
  TEST_EQUAL(result.precursor_ppm.has_value(), true)
  TEST_EQUAL(result.fragment_ppm.has_value(), true)
  TEST_EQUAL(result.diagnostics.fragment_profile_spectra, 5)
  TEST_EQUAL(result.diagnostics.fragment_profile_centroids >= 60, true)
  TEST_EQUAL(result.diagnostics.excluded_fragment_unknown, 0)

  estimator.reset();
  const auto reset_result = estimator.getResult();
  TEST_EQUAL(reset_result.precursor_ppm.has_value(), false)
  TEST_EQUAL(reset_result.diagnostics.precursor_eligible_ms2, 0)
}
END_SECTION

START_SECTION((std::optional<RobustError> getPrecursorPrecisionPPM() const))
{
  IDFreeMassErrorEstimator::Parameters parameters;
  parameters.min_precursor_cluster_size = 3;
  IDFreeMassErrorEstimator estimator(parameters);
  for (Size cycle = 0; cycle < 5; ++cycle)
  {
    estimator.consumeSpectrum(makeSpectrum(cycle, 0));
  }
  const auto precision = estimator.getPrecursorPrecisionPPM();
  TEST_EQUAL(precision.has_value(), true)
  TEST_EQUAL(precision->single_measurement_sigma > 0.0, true)
}
END_SECTION

END_TEST
