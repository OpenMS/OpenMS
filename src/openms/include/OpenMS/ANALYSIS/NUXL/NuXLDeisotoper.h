// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <string.h>

namespace OpenMS
{

class MSSpectrum;

/** @brief Deisotoping of nucleotide cross-link mass spectra.

  This class provides methods for deisotoping mass spectra. The algorithm
  considers each peak (starting from the right of a spectrum) and, for each
  peak
  **/
class OPENMS_DLLAPI NuXLDeisotoper
{
  public:

  /** @brief Detect isotopic clusters in a mass spectrum.

    Deisotoping is based on C13 abundance and will try to identify a simple
    model based on the C12-C13 distance and charge state. This is often a good
    approximation for peptide fragment ion spectra but may not work well for
    other spectra. The algorithm will consider each peak (starting from the
    right of a spectrum) and, for each peak, attempt to add isotopic peaks to
    its envelope until either no peak is found, the maximum number of isotopic
    peaks is reached or (only when using @p use_decreasing_model) the intensity
    of the peak is higher than the previous peak.

    Deisotoping is done in-place and if @p annotate_charge is true,
    an additional IntegerDataArray "charge" will be appended. If
    @p annotate_iso_peak_count is true, an additional IntegerDataArray
    "iso_peak_count" containing the number of isotopic peaks will be
    appended. If @p annotate_features is true, an addtional IntegerData Array
    "feature_number" containing the feature index will be appended.
    Existing DataArrays are kept and shrunken to the peaks which
    remain in the spectrum.

   * @param [spectrum] Input spectrum (sorted by m/z)
   * @param [fragment_tolerance] The tolerance used to match isotopic peaks
   * @param [fragment_unit_ppm] Whether ppm or m/z is used as tolerance
   * @param [min_charge] The minimum charge considered
   * @param [max_charge] The maximum charge considered
   * @param [keep_only_deisotoped] Only monoisotopic peaks of fragments with isotopic pattern are retained
   * @param [min_isopeaks] The minimum number of isotopic peaks (at least 2) required for an isotopic cluster
   * @param [max_isopeaks] The maximum number of isotopic peaks (at least 2) considered for an isotopic cluster
   * @param [make_single_charged] Convert deisotoped monoisotopic peak to single charge
   * @param [annotate_charge] Annotate the charge to the peaks in the IntegerDataArray: "charge" (0 for unknown charge)
   * @param [annotate_iso_peak_count] Annotate the number of isotopic peaks in a pattern for each monoisotopic peak in the IntegerDataArray: "iso_peak_count"
   * @param [use_decreasing_model] Use a simple averagine model that expects heavier isotopes to have less intensity. If false, no intensity checks are applied.
   * @param [start_intensity_check] Number of the isotopic peak from which the decreasing model should be applied. <= 1 will force the monoisotopic peak to be the most intense.
                                    2 will allow the monoisotopic peak to be less intense than the second peak.
                                    3 will allow the monoisotopic and the second peak to be less intense than the third, etc.
                                    A number higher than max_isopeaks will effectively disable use_decreasing_model completely.
   * @param [add_up_intensity] Sum up the total intensity of each isotopic pattern into the intensity of the reported monoisotopic peak
   * @param [annotate_features] Annotates the feature index in the IntegerDataArray: "feature_number".
   * @param [preserve_high_intensity_peaks] If true, the highest intensity peak of each isotopic pattern will never be filtered
   * @param [preserve_low_mz_peaks_threshold] If preserve_high_intensity_peaks is set, all peaks with smaller m/z will never be filtered
   *
   * Note: If @p make_single_charged is selected, the original charge (>=1) gets annotated.
   */
    static void deisotopeAndSingleCharge(MSSpectrum& spectrum,
                                         double fragment_tolerance,
                                         bool fragment_unit_ppm,
                                         int min_charge = 1,
                                         int max_charge = 3,
                                         bool keep_only_deisotoped = false,
                                         unsigned int min_isopeaks = 3,
                                         unsigned int max_isopeaks = 10,
                                         bool make_single_charged = true,
                                         bool annotate_charge = false,
                                         bool annotate_iso_peak_count = false,
                                         bool use_decreasing_model = true,
                                         unsigned int start_intensity_check = 2,
                                         bool add_up_intensity = false,
                                         bool annotate_features = false,
                                         bool preserve_high_intensity_peaks = false,
                                         double preserve_low_mz_peaks_threshold = -1e10);
};

}
