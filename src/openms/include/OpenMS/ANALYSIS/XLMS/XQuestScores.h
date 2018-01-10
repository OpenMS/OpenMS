// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_XLMS_XQUESTSCORES_H
#define OPENMS_ANALYSIS_XLMS_XQUESTSCORES_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Types.h>
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
    * @brief compute a simple and fast to compute pre-score for a cross-link spectrum match
    * @param number of experimental peaks matched to theoretical common ions from the alpha peptide
    * @param number of theoretical ions from the alpha peptide
    * @param number of experimental peaks matched to theoretical common ions from the beta peptide
    * @param number of theoretical ions from the beta peptide
    */
    static float preScore(Size matched_alpha, Size ions_alpha, Size matched_beta, Size ions_beta);

   /**
    * @brief compute a simple and fast to compute pre-score for a mono-link spectrum match
    * @param number of experimental peaks matched to theoretical common ions from the alpha peptide
    * @param number of theoretical ions from the alpha peptide
    */
    static float preScore(Size matched_alpha, Size ions_alpha);

   /**
    * @brief compute the match-odds score, a score based on the probability of getting the given number of matched peaks by chance
    * @param theoretical spectrum, sorted by position
    * @param alignment between the theoretical and the experimental spectra
    * @param fragment mass tolerance of the alignment
    * @param fragment mass tolerance unit of the alignment, true = ppm, false = Da
    * @param type of cross-link, true = cross-link, false = mono-link
    * @param number of considered charges in the theoretical spectrum
    */
    static double matchOddsScore(const PeakSpectrum& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm, bool is_xlink_spectrum = false, Size n_charges = 1);

    static double logOccupancyProb(const PeakSpectrum& theoretical_spec,  const Size matched_size, double fragment_mass_tolerance, bool fragment_mass_tolerance_unit_ppm);


   /**
    * @brief compute the weighted total ion current score for a cross-link. Reimplementation from xQuest.
    * @param sequence length of alpha peptide
    * @param sequence length of beta peptide
    * @param intensity sum of matched peaks from alpha peptide
    * @param intensity sum of matched peaks from beta peptide
    * @param type of cross-link, true = cross-link, false = mono-link
    * @param sum of peak intensities of the experimental spectrum
    * @param true = cross-link, false = mono-link. in case of a mono-link, beta_size and intsum_beta should be 0
    */
    static double weightedTICScoreXQuest(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link);

   /**
    * @brief compute the weighted total ion current score for a cross-link. Scaling changed from original xQuest.
    * @param sequence length of alpha peptide
    * @param sequence length of beta peptide
    * @param intensity sum of matched peaks from alpha peptide
    * @param intensity sum of matched peaks from beta peptide
    * @param type of cross-link, true = cross-link, false = mono-link
    * @param sum of peak intensities of the experimental spectrum
    * @param true = cross-link, false = mono-link. in case of a mono-link, beta_size and intsum_beta should be 0
    */
    static double weightedTICScore(Size alpha_size, Size beta_size, double intsum_alpha, double intsum_beta, double total_current, bool type_is_cross_link);

   /**
    * @brief computes sum of peak intensities of matched peaks for either the alpha or the beta peptide
    * @param alignment between common alpha or beta ions and common experimental peaks
    * @param alignment between xlink alpha or beta ions and xlink experimental peaks
    * @param experimental common ion spectrum
    * @param experimental xlink spectrum
    */
    static double matchedCurrentChain(const std::vector< std::pair< Size, Size > >& matched_spec_common, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks);

   /**
    * @brief computes sum of peak intensities of all matched peaks
    * @param alignment between common alpha ions and common experimental peaks
    * @param alignment between common beta ions and common experimental peaks
    * @param alignment between xlink alpha ions and xlink experimental peaks
    * @param alignment between xlink beta ions and xlink experimental peaks
    * @param experimental common ion spectrum
    * @param experimental xlink spectrum
    */
    static double totalMatchedCurrent(const std::vector< std::pair< Size, Size > >& matched_spec_common_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_common_beta, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_alpha, const std::vector< std::pair< Size, Size > >& matched_spec_xlinks_beta, const PeakSpectrum& spectrum_common_peaks, const PeakSpectrum& spectrum_xlink_peaks);

   /**
    * @brief computes a crude cross-correlation between two spectra. Crude, because it uses a static binsize based on a tolerance in Da and it uses equal intensities for all peaks
    * @param first spectrum
    * @param second spectrum
    * @param number of bins, that should be considered for shifting the second spectrum. the second spectrum is shifted from -maxshift to +maxshift of tolerance bins and a correlation is computed for each position.
    * @param tolerance or binsize in Da
    */
    static std::vector< double > xCorrelation(const PeakSpectrum & spec1, const PeakSpectrum & spec2, Int maxshift, double tolerance);

  };

}

#endif // XQUESTSCORES_H
