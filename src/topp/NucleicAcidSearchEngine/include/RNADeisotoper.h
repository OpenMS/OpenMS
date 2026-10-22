// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Sam Wein $
// $Authors: Sam Wein $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  /// FLASHDeconv isotope detection with RNA averagine; preserves observed peak intensities.
  /// Each worker must use its own instance. Input must be sorted MS2, in profile or centroided form.
  class RNADeisotoper
  {
  public:
    /**
      @brief Configure FLASH RNA isotope detection.
      @param[in] tolerance_ppm Fragment tolerance in ppm.
      @param[in] max_abs_charge Largest charge magnitude to consider, used unchanged.
      @param[in] positive Use positive ionization mode.
      @param[in] max_mass Largest neutral mass in Da (default: 1000). Explicit limits are used unchanged.
      @throws Exception::IllegalArgument If tolerance_ppm is not finite and positive, max_abs_charge is invalid,
      min_cos or min_snr is non-finite, or max_mass is non-finite or below the FLASH minimum mass (50 Da).
      @param[in] min_cos Minimum charge-specific isotope cosine.
      @param[in] min_snr Minimum charge-specific envelope SNR.
      @param[in] legacy_trim Reproduce stage 8's RNA averagine trimming.
    */
    RNADeisotoper(double tolerance_ppm, int max_abs_charge, bool positive,
                  double max_mass = 1000.0, double min_cos = 0.85, double min_snr = 0.25,
                  bool legacy_trim = false);
    /**
      @brief Remove confidently assigned isotope samples while preserving observed intensities.

      Keep unassigned peaks and envelopes without an observed monoisotope. Retain every
      monoisotope profile sample once. Ambiguous charges are not decharged.
      Higher-charge competitors within configured bounds may veto assignments but are not emitted
      above assignment bounds. Harmonic screening cannot see charges beyond max_abs_charge.
      @param[in,out] spectrum Sorted profile or centroided MS2 spectrum.
      @param[in] scan_number Scan index supplied to FLASH.
      @param[in] decharge Convert uniquely assigned monoisotopes to single charge.
    */
    void filter(MSSpectrum& spectrum, int scan_number, bool decharge);

    /// Preserve stage 8's normalized-input FLASH detection for RNOME reproduction.
    void filterLegacy(MSSpectrum& spectrum, int scan_number, bool decharge,
                      double intensity_scale);

  private:
    SpectralDeconvolution detector_;
    bool positive_;
    int max_abs_charge_;
    double min_cos_;
    double min_snr_;
  };
}
