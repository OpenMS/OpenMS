// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <numeric>
#include <vector>

namespace OpenMS
{
  /**
    @brief Spectrum preprocessing and theoretical-vs-experimental peak alignment helpers used by the OpenPepXL cross-link search engines.

    Static utilities (the class carries no state) that the cross-link identification
    workflows in OpenPepXL build on:
      - Pre-filter an MS2 dataset before searching — drop spectra with too few peaks or
        with precursor charges outside the configured range, normalise intensities, and
        (optionally) deisotope.
      - Merge two annotated spectra into a single peak list while preserving paired
        FloatDataArray / StringDataArray / IntegerDataArray annotations.
      - Align a theoretical and an experimental fragment spectrum honouring per-peak
        charge annotations, with optional intensity-ratio filtering.

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI OPXLSpectrumProcessingAlgorithms
  {
    public:

      /**
       * @brief Merge two annotated spectra into one peak list, preserving paired DataArrays.
       *
       * Peaks of @p first_spectrum and @p second_spectrum are concatenated and the result
       * is sorted by m/z. For each kind of DataArray (Float / String / Integer) the i-th
       * array of @p first_spectrum is paired with the i-th array of @p second_spectrum and
       * their contents are concatenated; the output array inherits its name from
       * @p first_spectrum's i-th array. Extra arrays present only in @p second_spectrum are
       * dropped — pairing is positional, not by name.
       *
       * Despite the non-const references in the signature, neither input is modified.
       *
       * @param[in,out] first_spectrum  Spectrum whose DataArray names and ordering define the output.
       * @param[in,out] second_spectrum Spectrum whose peaks and (positionally paired) DataArrays are appended.
       * @return Merged spectrum, sorted by m/z.
       */
      static PeakSpectrum mergeAnnotatedSpectra(PeakSpectrum & first_spectrum, PeakSpectrum & second_spectrum);

      /**
       * @brief Preprocess an MSExperiment for cross-link search and return the surviving MS2 spectra.
       *
       * Two phases — first the input @p exp is modified in place: zero-intensity peaks are
       * removed (ThresholdMower), intensities are normalised, and the spectra are sorted by
       * retention time. Then the MS2 spectra are iterated (OpenMP-parallelised loop) and
       * those that pass the per-spectrum filters are copied into a freshly built PeakMap
       * that is returned. MS1 spectra are not copied to the output.
       *
       * For unlabeled data (@p labeled is @c false) a spectrum is retained only if it has a
       * single precursor whose charge lies in @c [min_precursor_charge, max_precursor_charge]
       * and at least @c 2 @c * @c peptide_min_size peaks. Such spectra are further reduced
       * by a WindowMower with hardcoded settings (window size @c 100, keep @c 20 peaks per
       * window, @c "jump" mode). The final peak count must again exceed
       * @c 2 @c * @c peptide_min_size.
       *
       * For labeled data (@p labeled is @c true) the precursor and peak-count filters are
       * bypassed, the WindowMower is not applied, and every MS2 spectrum of the input is
       * present in the output — keeping spectrum indices stable across the heavy/light
       * pairing performed downstream via consensusXML.
       *
       * When @p deisotope is @c true, each surviving spectrum is run through
       * Deisotoper::deisotopeAndSingleCharge with a simple averagine model, charge range
       * @c [1,7], isopeak counts in @c [3,10]; charge and isotopic-peak counts are
       * annotated and monoisotopic intensity is summed. The deisotoped result is kept only
       * if it still exceeds the post-filter peak count or @p labeled is @c true.
       *
       * @param[in,out] exp                              Input data (MS1 + MS2). Modified in place.
       * @param[in]     fragment_mass_tolerance          Peak mass tolerance used by deisotoping (ignored if @p deisotope is @c false).
       * @param[in]     fragment_mass_tolerance_unit_ppm Interpret @p fragment_mass_tolerance as ppm (true) or Da (false).
       * @param[in]     peptide_min_size                 Lower bound on peak count: spectra must have at least @c 2 @c * @c peptide_min_size peaks both before and after the WindowMower / Deisotoper step.
       * @param[in]     min_precursor_charge             Minimum allowed precursor charge for unlabeled data.
       * @param[in]     max_precursor_charge             Maximum allowed precursor charge for unlabeled data.
       * @param[in]     deisotope                        If @c true, deisotope each surviving spectrum.
       * @param[in]     labeled                          If @c true, bypass precursor/peak-count filters and the WindowMower, keeping every MS2 spectrum.
       *
       * @return PeakMap containing only the preprocessed MS2 spectra that passed all filters; for @p labeled inputs this contains every input MS2 spectrum.
       */
      static PeakMap preprocessSpectra(PeakMap& exp, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, Size peptide_min_size, Int min_precursor_charge, Int max_precursor_charge, bool deisotope, bool labeled);

      /**
       * @brief Align a theoretical and an experimental fragment spectrum using charge annotations and an intensity-ratio cut-off.
       *
       * For each theoretical peak, the closest experimental peak inside the mass-tolerance
       * window is picked, restricted to peaks whose charge and intensity also pass the
       * per-peak filters.
       *
       * Tolerance window half-width: @c theo_mz @c * @p fragment_mass_tolerance @c * @c 1e-6
       * when @p fragment_mass_tolerance_unit_ppm is @c true, else @p fragment_mass_tolerance Da.
       *
       * Charge filter: a pair (theoretical charge @c tz, experimental charge @c ez) matches
       * when @c tz @c == @c ez or either side is @c 0 (treated as "unknown"). If
       * @p theo_charges or @p exp_charges is empty, the charge filter degrades to permissive.
       *
       * Intensity filter: a pair (theoretical intensity @c ti, experimental intensity @c ei)
       * matches when @c min(ti,ei) @c / @c max(ti,ei) @c > @p intensity_cutoff. Pass
       * @c 0 to disable.
       *
       * Both spectra must be sorted by m/z; @p alignment and @p ppm_error_array must be
       * empty on entry (precondition). When either spectrum is empty, the function returns
       * with no output.
       *
       * @param[out] alignment                       Receives (theo-index, exp-index) match pairs. Must be empty on entry.
       * @param[in]  fragment_mass_tolerance         Tolerance window half-width.
       * @param[in]  fragment_mass_tolerance_unit_ppm Interpret @p fragment_mass_tolerance as ppm (true) or Da (false).
       * @param[in]  theo_spectrum                   Theoretical spectrum (sorted by m/z).
       * @param[in]  exp_spectrum                    Experimental spectrum (sorted by m/z).
       * @param[in]  theo_charges                    Per-peak charges for @p theo_spectrum; an empty array disables the charge filter.
       * @param[in]  exp_charges                     Per-peak charges for @p exp_spectrum; an empty array disables the charge filter.
       * @param[out] ppm_error_array                 Receives per-match ppm errors @c (exp_mz @c - @c theo_mz) @c / @c theo_mz @c * @c 1e6. Must be empty on entry.
       * @param[in]  intensity_cutoff                Minimum smaller-over-larger intensity ratio for a match; @c 0 disables the intensity filter.
       */
      static void getSpectrumAlignmentFastCharge(
            std::vector<std::pair<Size, Size> > & alignment, double fragment_mass_tolerance,
            bool fragment_mass_tolerance_unit_ppm,
            const PeakSpectrum& theo_spectrum,
            const PeakSpectrum& exp_spectrum,
            const DataArrays::IntegerDataArray& theo_charges,
            const DataArrays::IntegerDataArray& exp_charges,
            DataArrays::FloatDataArray& ppm_error_array,
            double intensity_cutoff = 0.0);

            /**
             * @brief Align a SimplePeak-based theoretical spectrum to an experimental spectrum using charge annotations only.
             *
             * Mirror of @ref getSpectrumAlignmentFastCharge but without the
             * intensity-ratio filter and without ppm-error output, and with the
             * theoretical side expressed as a @c vector<SimpleTSGXLMS::SimplePeak>
             * (charges carried per peak inside the @c SimplePeak struct, not as a
             * separate DataArray). @p alignment is cleared on entry — it does not need
             * to start empty.
             *
             * Charge filter rule is identical to the fast-charge variant: @c (theo_charge
             * @c == @c exp_charge or either side is @c 0) matches. An empty @p exp_charges
             * array makes the charge filter permissive.
             *
             * @param[out] alignment                       Receives (theo-index, exp-index) match pairs. Cleared on entry.
             * @param[in]  fragment_mass_tolerance         Tolerance window half-width.
             * @param[in]  fragment_mass_tolerance_unit_ppm Interpret @p fragment_mass_tolerance as ppm (true) or Da (false).
             * @param[in]  theo_spectrum                   Theoretical spectrum (vector of SimplePeak; per-peak m/z and charge).
             * @param[in]  exp_spectrum                    Experimental spectrum (sorted by m/z).
             * @param[in]  exp_charges                     Per-peak charges for @p exp_spectrum; an empty array disables the charge filter.
             */
            static void getSpectrumAlignmentSimple(
                  std::vector<std::pair<Size, Size> > & alignment,
                  double fragment_mass_tolerance,
                  bool fragment_mass_tolerance_unit_ppm,
                  const std::vector< SimpleTSGXLMS::SimplePeak >& theo_spectrum,
                  const PeakSpectrum& exp_spectrum,
                  const DataArrays::IntegerDataArray& exp_charges);
  };

}
