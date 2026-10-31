/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "mzpeak/spectrum_metadata.h"

namespace MzPeak {

/**
 * Per-wavelength-spectrum (UV/Vis) descriptive metadata.
 *
 * On disk the primary struct column is named `spectrum`, exactly as for mass
 * spectra, but the observed-range and maximum fields describe a WAVELENGTH
 * axis in nanometres.  They are given their own names here rather than reusing
 * @ref SpectrumMetadata: a field called `lowest_observed_mz` holding 210 nm is
 * a lie that survives every code review.
 */
struct WavelengthSpectrumMetadata final {
  /// `spectrum.index` (uint64).
  uint64_t index = 0;

  /// `spectrum.id` — the native identifier.
  std::string id;

  /// Acquisition time in SECONDS.
  ///
  /// Sourced from `scan.MS_1000016_scan_start_time` when the scan facet is
  /// present, otherwise from `spectrum.time`; both are stored in MINUTES
  /// (UO:0000031) and converted here.  The two agree to float32 rounding
  /// (~7 µs on the bundled fixture).  Seconds is the OpenMS convention and
  /// matches @ref SpectrumMetadata::retention_time — mixing the two units is a
  /// silent 60x error.
  std::optional<double> time;

  /// `MS_1000559_spectrum_type` (CURIE).  Empty when null.
  std::string spectrum_type;

  /// `MS_1000525_spectrum_representation` (CURIE, e.g. "MS:1000128" profile).
  /// Empty when null.
  std::string representation;

  /// `MS_1000619_lowest_observed_wavelength`, in NANOMETRES (UO:0000018).
  ///
  /// @note The reference writer derives this from the input's first value
  /// AFTER sorting a separate copy for output, so on unsorted input it can
  /// disagree with the stored axis.  Trust the decoded array over this field
  /// when they conflict.
  std::optional<double> lowest_observed_wavelength;

  /// `MS_1000618_highest_observed_wavelength`, in NANOMETRES.  Carries the same
  /// caveat as @ref lowest_observed_wavelength.
  std::optional<double> highest_observed_wavelength;

  /// `MS_1003812_lambda_max` — wavelength of maximum intensity, in NANOMETRES.
  ///
  /// @note The reference writer seeds its running maximum at zero, so a
  /// spectrum whose intensities are all negative — ordinary for a
  /// baseline-subtracted absorbance trace — is written with lambda_max 0 and
  /// base_peak_intensity 0 rather than its true maximum.  A 0 here alongside a
  /// 0 base peak intensity should be read as "not computed".
  std::optional<double> lambda_max;

  /// `MS_1003060_number_of_data_points`.
  std::optional<uint64_t> number_of_data_points;

  /// `MS_1000505_base_peak_intensity`.  See the caveat on @ref lambda_max.
  std::optional<float> base_peak_intensity;

  /// `MS_1000285_total_ion_current`.
  std::optional<float> total_ion_current;

  /// `data_processing_ref`.  Empty string when null.
  std::string data_processing_ref;

  /// Per-spectrum CV parameters (`spectrum.parameters`).
  std::vector<CvParam> parameters;

  /// Scan-level CV parameters (`scan.parameters`).
  std::vector<CvParam> scan_parameters;
};

} // namespace MzPeak
