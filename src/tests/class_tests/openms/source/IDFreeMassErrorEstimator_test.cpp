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

  MSSpectrum makeContaminatedSpectrum(Size cycle, Size target)
  {
    auto spectrum = makeSpectrum(cycle, target, false, 0.0015);
    for (Size i = 0; i < 24; ++i)
    {
      const Int pseudo = static_cast<Int>((cycle * 37 + target * 17 + i * 29) % 101);
      const double shift = (static_cast<double>(pseudo) / 100.0 - 0.5) * 0.30;
      spectrum.emplace_back(650.0 + static_cast<double>(i) * 2.0 + shift, 70.0f + static_cast<float>(i));
    }
    spectrum.sortByPosition();
    return spectrum;
  }

  MSSpectrum makeAdaptiveFallbackSpectrum(Size cycle, Size target)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(2);
    spectrum.setRT(static_cast<double>(cycle * 10 + target));
    spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);

    Precursor precursor;
    precursor.setMZ(500.0 + static_cast<double>(target) * 20.0);
    precursor.setCharge(0);
    spectrum.setPrecursors({precursor});

    const double narrow_shift = std::sin(static_cast<double>(cycle) * 0.71 + static_cast<double>(target) * 0.37) * 0.0015;
    for (Size i = 0; i < 12; ++i)
    {
      spectrum.emplace_back(150.25 + static_cast<double>(i) * 35.0 + narrow_shift, 1000.0f + static_cast<float>(i));
    }

    for (Size i = 0; i < 30; ++i)
    {
      const double phase = static_cast<double>(cycle) * 0.91 + static_cast<double>(target) * 0.43 + static_cast<double>(i) * 0.61;
      const double broad_shift = std::sin(phase) * 0.035;
      spectrum.emplace_back(650.25 + static_cast<double>(i) * 3.0 + broad_shift, 70.0f + static_cast<float>(i));
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


START_SECTION((out-of-order spectra do not pair outside the RT window))
{
  IDFreeMassErrorEstimator::Parameters parameters;
  parameters.rt_window_seconds = 5.0;

  auto later = makeSpectrum(0, 0);
  later.setRT(100.0);
  auto earlier = makeSpectrum(1, 0);
  earlier.setRT(0.0);
  auto nearby = makeSpectrum(2, 0);
  nearby.setRT(1.0);

  MSExperiment experiment;
  experiment.addSpectrum(later);
  experiment.addSpectrum(earlier);
  experiment.addSpectrum(nearby);

  IDFreeMassErrorEstimator estimator(parameters);
  const auto result = estimator.compute(experiment);

  // RT=100 must not join the precursor cluster formed at RT=0/1.
  // The two nearby observations do not reach the default three-spectrum
  // precursor-cluster support threshold, so no precursor pair is committed.
  TEST_EQUAL(result.diagnostics.precursor_paired_spectra, 0)

  // RT=0 and RT=1 are a legitimate repeated-spectrum pair; RT=100 is not.
  TEST_EQUAL(result.diagnostics.fragment_paired_spectra, 1)
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

START_SECTION((charge-less DIA spectra contribute fragment precision but not precursor precision))
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
      auto spectrum = makeSpectrum(cycle, target);
      auto precursors = spectrum.getPrecursors();
      precursors.front().setCharge(0);
      spectrum.setPrecursors(precursors);
      estimator.consumeSpectrum(spectrum);
    }
  }

  const auto result = estimator.getResult();
  TEST_EQUAL(result.precursor_ppm.has_value(), false)
  TEST_EQUAL(result.precursor_da.has_value(), false)
  TEST_EQUAL(result.precursor_tolerance_ppm.has_value(), false)
  TEST_EQUAL(result.fragment_ppm.has_value(), true)
  TEST_EQUAL(result.fragment_tolerance_ppm.has_value(), true)
  TEST_EQUAL(result.fragment_tolerance_da.has_value(), false)
  TEST_EQUAL(result.fragment_resolution_regime, IDFreeMassErrorEstimator::FragmentResolutionRegime::HIGH_RESOLUTION)
  TEST_EQUAL(result.diagnostics.precursor_eligible_ms2, 0)
  TEST_EQUAL(result.diagnostics.fragment_eligible_ms2, 32)
  TEST_EQUAL(result.diagnostics.fragment_centroid_spectra, 32)
  TEST_EQUAL(result.diagnostics.excluded_missing_precursor, 0)
}
END_SECTION

START_SECTION((background-aware fragment fit recovers a narrow component from broad false matches))
{
  IDFreeMassErrorEstimator::Parameters parameters;
  parameters.min_fragment_pairs = 20;
  parameters.min_fragment_tolerance_pairs = 20;
  parameters.min_fragment_tolerance_spectra = 2;

  IDFreeMassErrorEstimator estimator(parameters);
  for (Size cycle = 0; cycle < 12; ++cycle)
  {
    for (Size target = 0; target < 4; ++target)
    {
      auto spectrum = makeContaminatedSpectrum(cycle, target);
      auto precursors = spectrum.getPrecursors();
      precursors.front().setCharge(0);
      spectrum.setPrecursors(precursors);
      estimator.consumeSpectrum(spectrum);
    }
  }

  const auto result = estimator.getResult();
  TEST_EQUAL(result.fragment_resolution_regime, IDFreeMassErrorEstimator::FragmentResolutionRegime::HIGH_RESOLUTION)
  TEST_EQUAL(result.fragment_ppm.has_value(), true)
  TEST_EQUAL(result.fragment_tolerance_ppm.has_value(), true)
  TEST_EQUAL(result.diagnostics.fragment_mixture_converged, true)
  TEST_EQUAL(result.diagnostics.fragment_mixture_signal_fraction > 0.1, true)
  TEST_EQUAL(result.diagnostics.fragment_mixture_signal_fraction < 0.95, true)
  TEST_EQUAL(result.fragment_ppm->single_measurement_sigma < 10.0, true)
  TEST_EQUAL(result.diagnostics.fragment_mixture_signal_pairs >= parameters.min_fragment_tolerance_pairs, true)
}
END_SECTION


START_SECTION((high-intensity fallback rescues a narrow fragment component from weak broad evidence))
{
  IDFreeMassErrorEstimator::Parameters parameters;
  parameters.min_fragment_pairs = 20;
  parameters.min_fragment_tolerance_pairs = 20;
  parameters.min_fragment_tolerance_spectra = 2;
  // Keep this synthetic fixture focused on the fallback branch: its full fit is
  // slightly broader than its high-intensity fit, so place the configured gate
  // between the two instead of relying on the production 10 ppm threshold.
  parameters.fragment_high_res_max_sigma_ppm = 2.75;

  IDFreeMassErrorEstimator estimator(parameters);
  for (Size cycle = 0; cycle < 16; ++cycle)
  {
    for (Size target = 0; target < 4; ++target)
    {
      estimator.consumeSpectrum(makeAdaptiveFallbackSpectrum(cycle, target));
    }
  }

  const auto result = estimator.getResult();
  TEST_EQUAL(result.fragment_resolution_regime, IDFreeMassErrorEstimator::FragmentResolutionRegime::HIGH_RESOLUTION)
  TEST_EQUAL(result.fragment_ppm.has_value(), true)
  TEST_EQUAL(result.fragment_tolerance_ppm.has_value(), true)
  TEST_EQUAL(result.diagnostics.fragment_high_intensity_pairs > 0, true)
  TEST_EQUAL(result.diagnostics.fragment_high_intensity_signal_pairs >= parameters.min_fragment_pairs, true)
  TEST_EQUAL(result.diagnostics.fragment_high_intensity_mixture_converged, true)
  TEST_EQUAL(result.diagnostics.fragment_high_intensity_fallback_used, true)
  TEST_EQUAL(result.fragment_ppm->single_measurement_sigma < parameters.fragment_high_res_max_sigma_ppm, true)
}
END_SECTION

START_SECTION((quantized fragment differences fail closed))
{
  IDFreeMassErrorEstimator::Parameters parameters;
  parameters.min_fragment_pairs = 10;
  parameters.min_fragment_tolerance_pairs = 10;
  parameters.min_fragment_tolerance_spectra = 2;

  IDFreeMassErrorEstimator estimator(parameters);
  for (Size cycle = 0; cycle < 8; ++cycle)
  {
    for (Size target = 0; target < 4; ++target)
    {
      auto spectrum = makeSpectrum(cycle, target, false, 0.0);
      auto precursors = spectrum.getPrecursors();
      precursors.front().setCharge(0);
      spectrum.setPrecursors(precursors);
      estimator.consumeSpectrum(spectrum);
    }
  }

  const auto result = estimator.getResult();
  TEST_EQUAL(result.diagnostics.fragment_pairs > 0, true)
  TEST_EQUAL(result.diagnostics.fragment_zero_delta_fraction >= 0.5, true)
  TEST_EQUAL(result.diagnostics.fragment_mixture_rejected_zero_quantization, true)
  TEST_EQUAL(result.diagnostics.fragment_mixture_converged, false)
  TEST_EQUAL(result.fragment_resolution_regime, IDFreeMassErrorEstimator::FragmentResolutionRegime::UNAVAILABLE)
  TEST_EQUAL(result.fragment_tolerance_ppm.has_value(), false)
  TEST_EQUAL(result.fragment_tolerance_da.has_value(), false)
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
