// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Sam Wein $
// $Authors: Sam Wein $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/SpectralDeconvolution.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  /** @brief Correct MS2 precursor m/z from FLASH RNA envelopes in the parent MS1 scan.

    Construct precursor-related FLASH envelopes using the complete parent MS1.
    Match the selected ion to an observed isotope at the reported charge.
    Unknown charges can be inferred from a unique envelope. Missing monoisotopes
    are allowed. FLASH supplies the isotope index; subtract its mass shift from
    the instrument's selected-ion m/z to preserve the original mass precision.
    Ambiguous or unmatched precursors retain their original m/z and charge.
    Accepted corrections preserve the input charge sign convention.
  */
  class RNAPrecursorCorrector
  {
  public:
    struct Result
    {
      std::string spectrum, parent, parent_source, status = "no_envelope";
      double original_mz = 0, corrected_mz = 0, mono_mass = 0, cosine = 0, snr = 0;
      double selected_error_ppm = 0;
      int original_charge = 0, charge = 0, isotope = -1;
      Size parent_peak_count = 0, parent_query_count = 0;
      bool mono_observed = false;

      /// Whether a unique confident isotope or monoisotope was assigned.
      bool accepted() const { return status == "corrected" || status == "monoisotopic"; }
    };

    /**
      @brief Configure MS1-supported RNA precursor correction.
      @param[in] tolerance_ppm Selected-ion and envelope tolerance in ppm.
      @param[in] min_abs_charge Smallest charge magnitude to consider.
      @param[in] max_abs_charge Largest charge magnitude to consider.
      @param[in] positive Use positive ionization mode.
      @param[in] max_mass Largest neutral mass in Da (default: 1000). Explicit limits are used unchanged.
      @throws Exception::IllegalArgument If tolerance_ppm is not finite and positive, charge magnitudes
      are invalid, min_cos or min_snr is non-finite, or max_mass is non-finite or below the FLASH minimum mass (50 Da).
      @param[in] min_cos Minimum charge-specific isotope cosine.
      @param[in] min_snr Minimum charge-specific envelope SNR.
    */
    RNAPrecursorCorrector(double tolerance_ppm, int min_abs_charge, int max_abs_charge,
                          bool positive, double max_mass = 1000.0,
                          double min_cos = 0.85, double min_snr = 0.25);

    /**
      @brief Select a unique envelope independently of RNA sequence matches or scores.
      @param[in] groups RNA envelopes detected from complete parent MS1 evidence.
      @param[in] precursor Original selected ion and reported charge.
      @return Confident correction or an unresolved result preserving the input.
    */
    Result match(const DeconvolvedSpectrum& groups, const Precursor& precursor) const;

    /**
      @brief Correct selected-ion values without moving acquisition isolation windows.

      Batch all selected ions sharing a parent MS1 and deconvolve that parent once.
      Retain correction results only, without caching raw MS1 copies or envelopes.
      MS1 peaks may be unsorted; any sorting uses a copy and preserves the original scan.
      Prefer spectrum_ref; invalid references fall back to the preceding MS1 and are reported.
      Replace earlier FLASH annotations; unresolved results omit accepted-only fields.
      @param[in,out] experiment Original MS1/MS2 experiment in acquisition order.
      @return One audit result for each MS2, including unresolved scans.
    */
    std::vector<Result> correct(MSExperiment& experiment);

    /**
      @brief Write one tab-delimited audit row per MS2, including unresolved scans.
      @param[in] filename Audit output path.
      @param[in] results Correction results in acquisition order.
    */
    static void storeReport(const std::string& filename, const std::vector<Result>& results);

  private:
    SpectralDeconvolution detector_;
    double tolerance_ppm_, min_cos_, min_snr_;
    int min_abs_charge_, max_abs_charge_;
    bool positive_;
  };
}
