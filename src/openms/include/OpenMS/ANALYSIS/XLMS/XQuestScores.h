// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <vector>

namespace OpenMS
{

/**
 *  @brief An implementation of the scores for cross-link identification from the xQuest algorithm (O. Rinner et al., 2008, "Identification of cross-linked peptides from large sequence databases")
 */

  class OPENMS_DLLAPI XQuestScores
  {

  public:
   /**
     @brief compute a simple and fast to compute pre-score for a cross-link spectrum match
     @param matched_alpha number of experimental peaks matched to theoretical linear ions from the alpha peptide
     @param ions_alpha number of theoretical ions from the alpha peptide
     @param matched_beta number of experimental peaks matched to theoretical linear ions from the beta peptide
     @param ions_beta number of theoretical ions from the beta peptide
    */
    static float preScore(Size matched_alpha, Size ions_alpha, Size matched_beta, Size ions_beta);

   /**
     @brief compute a simple and fast to compute pre-score for a mono-link spectrum match
     @param matched_alpha number of experimental peaks matched to theoretical linear ions from the alpha peptide
     @param ions_alpha number of theoretical ions from the alpha peptide
    */
    static float preScore(Size matched_alpha, Size ions_alpha);

   /**
     @brief compute the match-odds score, a score based on the probability of getting the given number of matched peaks by chance
     @param theoretical_spec theoretical spectrum, sorted by position
     @param matched_size alignment between the theoretical and the experimental spectra
     @param fragment_mass_tolerance fragment mass tolerance of the alignment
     @param fragment_mass_tolerance_unit_ppm fragment mass tolerance unit of the alignment, true = ppm, false = Da
     @param is_xlink_spectrum type of cross-link, true = cross-link, false = mono-link
     @param n_charges number of considered charges in the theoretical spectrum
    */
    static double matchOddsScore(const PeakSpectrum& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum = false, Size n_charges = 1);

    static double matchOddsScoreSimpleSpec(const std::vector< SimpleTSGXLMS::SimplePeak >& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum = false, Size n_charges = 1);


    /**
      @brief compute the logOccupancyProb score, similar to the match_odds, a score based on the probability of getting the given number of matched peaks by chance
      @param theoretical_spec theoretical spectrum, sorted by position
      @param matched_size number of matched peaks between experimental and theoretical spectra
      @param fragment_mass_tolerance the tolerance of the alignment
      @param fragment_mass_tolerance_unit_ppm the tolerance unit of the alignment, true = ppm, false = Da
     */
    static double logOccupancyProb(const PeakSpectrum& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm);

   /**
    * @brief compute the weighted total ion current score for a cross-link. Reimplementation from xQuest.
    * @param alpha_size sequence length of alpha peptide
    * @param beta_size  sequence length of beta peptide
    * @param intsum_alpha intensity sum of matched peaks from alpha peptide
    * @param intsum_beta intensity sum of matched peaks from beta peptide
    * @param total_current sum of peak intensities of the experimental spectrum
    * @param type_is_cross_link type of cross-link, true = cross-link, false = mono-link
    * @return true = cross-link, false = mono-link. in case of a mono-link, beta_size and intsum_beta should be 0
    */
    static double weightedTICScoreXQuest(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link);

   /**
    * @brief compute the weighted total ion current score for a cross-link. Scaling changed from original xQuest.
    * @param alpha_size sequence length of alpha peptide
    * @param beta_size sequence length of beta peptide
    * @param intsum_alpha intensity sum of matched peaks from alpha peptide
    * @param intsum_beta intensity sum of matched peaks from beta peptide
    * @param total_current Sum of peak intensities of the experimental spectrum
    * @param type_is_cross_link Type of cross-link, true = cross-link, false = mono-link
    * @return true = cross-link, false = mono-link. in case of a mono-link, beta_size and intsum_beta should be 0
    */
    static double weightedTICScore(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link);

   /**
    * @brief computes sum of peak intensities of matched peaks for either the alpha or the beta peptide
    * @param matched_spec_linear alignment between linear alpha or beta ions and linear experimental peaks
    * @param matched_spec_xlinks alignment between xlink alpha or beta ions and xlink experimental peaks
    * @param spectrum_linear_peaks experimental linear ion spectrum
    * @param spectrum_xlink_peaks experimental xlink spectrum
    */
    static double matchedCurrentChain(const std::vector< std::pair< Size, Size > >& matched_spec_linear, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_linear_peaks, const PeakSpectrum& spectrum_xlink_peaks);

   /**
    * @brief computes sum of peak intensities of all matched peaks
    * @param matched_spec_linear_alpha alignment between linear alpha ions and linear experimental peaks
    * @param matched_spec_linear_beta alignment between linear beta ions and linear experimental peaks
    * @param matched_spec_xlinks_alpha alignment between xlink alpha ions and xlink experimental peaks
    * @param matched_spec_xlinks_beta alignment between xlink beta ions and xlink experimental peaks
    * @param spectrum_linear_peaks experimental linear ion spectrum
    * @param spectrum_xlink_peaks experimental xlink spectrum
    */
    static double totalMatchedCurrent(const std::vector< std::pair< Size, Size > >& matched_spec_linear_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_linear_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_linear_peaks, const PeakSpectrum& spectrum_xlink_peaks);

   /**
    * @brief computes a crude cross-correlation between two spectra. Crude, because it uses a static binsize based on a tolerance in Da and it uses equal intensities for all peaks
    * @param spec1 first spectrum
    * @param spec2 second spectrum
    * @param maxshift Number of bins, that should be considered for shifting the second spectrum. the second spectrum is shifted from -maxshift to +maxshift of tolerance bins and a correlation is computed for each position.
    * @param tolerance or binsize in Da
    */
    static std::vector< double > xCorrelation(const PeakSpectrum & spec1, const PeakSpectrum & spec2, Int maxshift, double tolerance);

    /**
     * @brief computes a crude dot product between two spectra. Crude, because it uses a static binsize based on a tolerance in Da and it uses equal intensities for all peaks
     * @param spec1 first spectrum
     * @param spec2 second spectrum
     * @param tolerance tolerance or binsize in Da
     */
    static double xCorrelationPrescore(const PeakSpectrum & spec1, const PeakSpectrum & spec2, double tolerance);
  };

}
