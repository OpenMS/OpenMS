// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Macros.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <vector>

namespace OpenMS
{

/**
 *  @brief An implementation of the X!Tandem HyperScore PSM scoring function
 */               

struct OPENMS_DLLAPI HyperScore
{
  typedef std::pair<Size, double> IndexScorePair; 

  /** @brief compute the (ln transformed) X!Tandem HyperScore 
   *  1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
   *  2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
   * @note Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
   * @param[in] fragment_mass_tolerance mass tolerance applied left and right of the theoretical spectrum peak position
   * @param[in] fragment_mass_tolerance_unit_ppm Unit of the mass tolerance is: Thomson if false, ppm if true
   * @param[in] exp_spectrum measured spectrum
   * @param[in] theo_spectrum theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator.
   */
//  static double compute(double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, const PeakSpectrum& exp_spectrum, const RichPeakSpectrum& theo_spectrum);

  static double compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const PeakSpectrum& theo_spectrum);

  /** @brief compute the (ln transformed) X!Tandem HyperScore 
   *  overload that returns some additional information on the match
   */
  struct PSMDetail
  {
    size_t matched_prefix_ions = 0;  ///< N-terminal ions (a, b, c)
    size_t matched_suffix_ions = 0;  ///< C-terminal ions (x, y, z)
    double mean_error = 0.0;
  };

  static double computeWithDetail(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const PeakSpectrum& theo_spectrum,
                        PSMDetail& d
                       );

  /**
   * @brief Experimental intensity score with binomial fragment-match evidence.
   *
   * Replaces HyperScore's factorial rewards with negative log binomial tails,
   * accounting for the number of theoretical ions and experimental peak density.
   * Fragment matches are dependent, so this score is not a calibrated PSM p-value.
   * Spectra must be sorted by m/z and theoretical peaks need IonNames annotations.
   *
   * @param[in] fragment_mass_tolerance Fragment matching tolerance (Da or ppm).
   * @param[in] fragment_mass_tolerance_unit_ppm Whether tolerance is in ppm.
   * @param[in] exp_spectrum Experimental spectrum.
   * @param[in] theo_spectrum Annotated theoretical spectrum.
   * @param[out] detail Counts and mean absolute error for the selected matches.
   * @return Natural-log intensity plus prefix/suffix match evidence; zero without matches.
   * @throws Exception::InvalidParameter if the tolerance is not finite and positive.
   * @throws Exception::InvalidValue if the theoretical ion annotations are missing or incomplete.
   */
  static double computeCalibrated(double fragment_mass_tolerance,
                                 bool fragment_mass_tolerance_unit_ppm,
                                 const PeakSpectrum& exp_spectrum,
                                 const PeakSpectrum& theo_spectrum,
                                 PSMDetail& detail);

  /**
   * @brief Experimental HyperScore with mass-accuracy-weighted fragment evidence.
   *
   * A matched ion with signed mass error e (ppm) contributes exp(-0.5*((e-shift)/sigma)^2)
   * to its terminal ion count and intensity product. Fractional factorial rewards
   * use max(0, lgamma(count+1)), agreeing with HyperScore for matches at the kernel center.
   * The kernel center defaults to zero; a search can center it on the systematic fragment
   * error of a run. This ranking statistic is not a PSM p-value.
   *
   * @param[in] fragment_mass_tolerance Positive finite match tolerance (Da or ppm).
   * @param[in] fragment_mass_tolerance_unit_ppm Whether matching tolerance is in ppm.
   * @param[in] exp_spectrum Experimental spectrum, sorted by m/z.
   * @param[in] theo_spectrum Theoretical spectrum, sorted by m/z, with ion names.
   * @param[in] mass_error_sd_ppm Positive finite Gaussian standard deviation in ppm.
   * @param[out] detail Unweighted match counts and mean absolute error in matching units.
   * @param[in] mass_error_shift_ppm Finite kernel center in ppm (signed, observed minus theoretical).
   * @return Nonnegative weighted log HyperScore; zero without matches.
   * @throws Exception::InvalidParameter for invalid tolerance, standard deviation or shift.
   * @throws Exception::InvalidValue for incomplete theoretical ion annotations.
   */
  static double computeMassAccuracy(double fragment_mass_tolerance,
                                    bool fragment_mass_tolerance_unit_ppm,
                                    const PeakSpectrum& exp_spectrum,
                                    const PeakSpectrum& theo_spectrum,
                                    double mass_error_sd_ppm,
                                    PSMDetail& detail,
                                    double mass_error_shift_ppm = 0.0);

  /**
   * @brief Signed ppm errors (observed minus theoretical) of the theoretical ions matched within the tolerance.
   *
   * Uses the same closest-peak matching as computeWithDetail() and computeMassAccuracy(), so the
   * errors describe exactly the matches those scorers count. One value per matched theoretical
   * ion is appended to @p errors_ppm; nothing is appended for unmatched ions or empty spectra.
   * Both spectra must be sorted by m/z.
   *
   * @param[in] fragment_mass_tolerance Positive finite match tolerance (Da or ppm).
   * @param[in] fragment_mass_tolerance_unit_ppm Whether the tolerance is in ppm.
   * @param[in] exp_spectrum Experimental spectrum, sorted by m/z.
   * @param[in] theo_spectrum Theoretical spectrum, sorted by m/z.
   * @param[out] errors_ppm Receives the signed errors of the matched ions (appended).
   * @throws Exception::InvalidParameter if the tolerance is not finite and positive.
   */
  static void matchedFragmentErrorsPpm(double fragment_mass_tolerance,
                                       bool fragment_mass_tolerance_unit_ppm,
                                       const PeakSpectrum& exp_spectrum,
                                       const PeakSpectrum& theo_spectrum,
                                       std::vector<double>& errors_ppm);

  /* @brief compute the (ln transformed) X!Tandem HyperScore only matching peaks that match in charge
   *  1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
   *  2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
   * @note Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
   * @param[in] fragment_mass_tolerance mass tolerance applied left and right of the theoretical spectrum peak position
   * @param[in] fragment_mass_tolerance_unit_ppm Unit of the mass tolerance is: Thomson if false, ppm if true
   * @param[in] exp_spectrum measured spectrum
   * @param[in] exp_charges charges of measured peaks
   * @param[in] theo_spectrum theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator.
   * @param[in] theo_charges charges of theoretical peaks
  */
  static double compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const DataArrays::IntegerDataArray& exp_charges,
                        const PeakSpectrum& theo_spectrum,
                        const DataArrays::IntegerDataArray& theo_charges);

  /* @brief compute the (ln transformed) X!Tandem HyperScore only matching peaks that match in charge
   *  1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
   *  2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
   * @note Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
   * @param[in] fragment_mass_tolerance mass tolerance applied left and right of the theoretical spectrum peak position
   * @param[in] fragment_mass_tolerance_unit_ppm Unit of the mass tolerance is: Thomson if false, ppm if true
   * @param[in] exp_spectrum measured spectrum
   * @param[in] exp_charges charges of measured peaks
   * @param[in] theo_spectrum theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator.
   * @param[in] theo_charges charges of theoretical peaks
   * @param[in] intensity_sum summed intensity for observed bond indices (e.g., b3=123 -> intensity_sum[2]=123)
   * Note: intensity_sum must be zeroed and of size #AA in peptide
  */
  static double compute(double fragment_mass_tolerance, 
                        bool fragment_mass_tolerance_unit_ppm, 
                        const PeakSpectrum& exp_spectrum, 
                        const DataArrays::IntegerDataArray& exp_charges,
                        const PeakSpectrum& theo_spectrum,
                        const DataArrays::IntegerDataArray& theo_charges,
                        std::vector<double>& intensity_sum);

  private:
    /// helper to compute the log factorial
    static double logfactorial_(const int x, int base = 2);
};

}


