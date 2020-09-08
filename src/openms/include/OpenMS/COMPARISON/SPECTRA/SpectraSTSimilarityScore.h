// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>

namespace OpenMS
{

  /**
      @brief Similarity score of SpectraST.

      Unlike the other similarity scores this score is used for matching a
      spectrum against a whole library, although the dot product seems to be an
      effective method for scoring on its own. For calculating the SpectraST
      score, first preprocess the spectra if not already done. Transform them
      and calculate the dot product and the dot bias.  Afterwards get the best
      two hits and calculate delta_D. Now for every spectrum from the library
      you can calculate the final score.

      The details of the score can be found in:
      H. Lam et al., Development and validation of a spectral library searching
      method for peptide identification from MS/MS,
      Proteomics, 7 , 655-667, 2007

      @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI SpectraSTSimilarityScore :
    public PeakSpectrumCompareFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    SpectraSTSimilarityScore();

    /// copy constructor
    SpectraSTSimilarityScore(const SpectraSTSimilarityScore & source);

    /// destructor
    ~SpectraSTSimilarityScore() override;
    // @}

    /// assignment operator
    SpectraSTSimilarityScore & operator=(const SpectraSTSimilarityScore & source);

    /**
        @brief: calculates the dot product of the two spectra
    */
    double operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const override;
    /**
        @brief: calculates the dot product of the two spectra
    */
    double operator()(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2)   const;
    /**
        @brief: calculates the dot product of itself
    */
    double operator()(const PeakSpectrum & spec) const override;

    /**
        @brief Preprocesses the spectrum

        The preprocessing removes peak below a intensity threshold, reject spectra that does
        not have enough peaks, and cuts peaks exceeding the max_peak_number most intense peaks.

        @return true if spectrum passes filtering
    */
    bool preprocess(PeakSpectrum & spec, float remove_peak_intensity_threshold = 2.01, UInt cut_peaks_below = 1000, Size min_peak_number = 5, Size max_peak_number = 150);


    ///spectrum is transformed into a binned spectrum with bin size 1 and spread 1 and the intensities are normalized.
    BinnedSpectrum transform(const PeakSpectrum & spec);

    /**
        @brief Calculates how much of the dot product is dominated by a few peaks

        @param dot_product if -1 this value will be calculated as well.
        @param bin1 first spectrum in binned representation
        @param bin2 second spectrum in binned representation
    */
    double dot_bias(const BinnedSpectrum & bin1, const BinnedSpectrum & bin2, double dot_product = -1) const;

    /**
        @brief calculates the normalized distance between top_hit and runner_up.
        @param top_hit is the best score for a given match.
        @param runner_up a match with a worse score than top_hit, e.g. the second best score.

        @return normalized distance
        @throw DividedByZero exception if top_hit is 0.

        @note Range of the dot products is between 0 and 1.
    */
    double delta_D(double top_hit, double runner_up);

    /**
        @brief: computes the overall all score
        @param dot_product of a match
        @param delta_D should be calculated after all dot products for a unidentified spectrum are computed
        @param dot_bias

        @return the SpectraST similarity score
    */
    double compute_F(double dot_product, double delta_D, double dot_bias);



    ///
    static PeakSpectrumCompareFunctor * create() { return new SpectraSTSimilarityScore(); }

    ///Reimplemented from PeakSpectrumCompareFunctor.
    static const String getProductName()
    {
      return "SpectraSTSimilarityScore";
    }

protected:


  };

}
