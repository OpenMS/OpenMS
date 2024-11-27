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

namespace OpenMS
{

class MSSpectrum;

class OPENMS_DLLAPI Deisotoper
{
  public:

  /** @brief Detect isotopic clusters in a mass spectrum.
    
    This algorithm is in parts taken from Guo Ci Teo et al, DOI: 10.1021/acs.jproteome.0c00544
    and is closely related to deisotopeAndSingleCharge by Timo Sachsenberg.

    Deisotoping is based on C13 abundance and will try to identify isotopic
    clusters fitting to an averagine model, taking in account the corresponding
    charge state. This only makes sense for peptide fragment ion spectra.
    The algorithm considers each peak in the spectrum and will try to form
    isotopic clusters for every charge state from @p min_charge to @p max_charge.
    Of the clusters that pass an averaginge check based on KL-divergence (for all 
    subclusters starting at the base peak as well), the cluster with most peaks (and
    in case of equality, also highest charge) is kept.

    Deisotoping is done in-place, and the algorithm removes peaks with 
    intensity 0. If @p rem_low_intensity is true,
    peaks not belonging to the highest 1000/5000 peaks are removed (see 
    @p used_for_open_search). If @p annotate_charge is true,
    an additional IntegerDataArray "charge" will be appended. If
    @p annotate_iso_peak_count is true, an additional IntegerDataArray
    "iso_peak_count" containing the number of isotopic peaks will be
    appended. Existing DataArrays are kept and shrunken to the peaks which
    remain in the spectrum.

   * @param [spectrum] Input spectrum (sorted by m/z)
   * @param [fragment_tolerance] The tolerance used to match isotopic peaks
   * @param [fragment_unit_ppm] Whether ppm or m/z is used as tolerance 
   * @param [number_of_final_peaks] Only the largest @p number_of_final_peaks peaks are kept in any spectrum. If 0, no filtering is performed. For open search, 1000 is recommended, else 5000.
   * @param [min_charge] The minimum charge considered
   * @param [max_charge] The maximum charge considered
   * @param [keep_only_deisotoped] Only monoisotopic peaks of fragments with isotopic pattern are retained
   * @param [min_isopeaks] The minimum number of isotopic peaks (at least 2) required for an isotopic cluster
   * @param [max_isopeaks] The maximum number of isotopic peaks (at least 2) considered for an isotopic cluster
   * @param [make_single_charged] Convert deisotoped monoisotopic peak to single charge, the original charge (>=1) gets annotated
   * @param [annotate_charge] Annotate the charge to the peaks in the IntegerDataArray: "charge" (0 for unknown charge)
   * @param [annotate_iso_peak_count] Annotate the number of isotopic peaks in a pattern for each monoisotopic peak in the IntegerDataArray: "iso_peak_count"
   * @param [add_up_intensity] Sum up the total intensity of each isotopic pattern into the intensity of the reported monoisotopic peak
   */
    static void deisotopeWithAveragineModel(MSSpectrum& spectrum,
                                            double fragment_tolerance,
                                            bool fragment_unit_ppm,
                                            int number_of_final_peaks = 5000,
                                            int min_charge = 1,
                                            int max_charge = 3,
                                            bool keep_only_deisotoped = false,
                                            unsigned int min_isopeaks = 2,
                                            unsigned int max_isopeaks = 10,
                                            bool make_single_charged = true,
                                            bool annotate_charge = false,
                                            bool annotate_iso_peak_count = false,
                                            bool add_up_intensity = false);

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
                                         bool annotate_features = false);
};

}
