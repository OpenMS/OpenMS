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
// $Maintainer: Timo Sachsenberg $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_STEINSCOTTIMPROVESCORE_H
#define OPENMS_COMPARISON_SPECTRA_STEINSCOTTIMPROVESCORE_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <cmath>

namespace OpenMS
{
  /**
      @brief Similarity score based of Stein & Scott

      This is a pairwise based score function. The spectrum contains peaks, and
      each peak can be defined by two values (mz and the intensity).  The score
      function takes the sum of the product of the peak intensities from
      Spectrum 1 and Spectrum 2, only if the mz-ratio distance between the two
      spectra is smaller than a given window size. In the default status, the
      window size is (accuracy of the mass spectrometer). This sum is
      normalised by dividing it with a distance function. sqrt(sum of the
      Intensity of square Spectrum1 sum of the Intensity of square Spectrum2).
      This is all based on SteinScott score.  To distinguish the close from the
      distant spectra an additional term is calculated. It denotes the expected
      value of both Spectra under the random placement of all peaks, within
      the given mass-to-charge range. The probability that two peaks with
      randomized intensity values lie within two epsilon of each other is a
      constant. This constant is proportional to epsilon. So the additional
      term is the sum over all peaks of Spectrum 1 and Spectrum 2 of the
      products of their intensities multiplied with the constant.


      The details of the score can be found in:
      Signal Maps for Mass Spectrometry-based
      Comparative Proteomics
      Amol Prakash, Parag Mallick , Jeffrey Whiteaker, Heidi Zhang,
      Amanda Paulovich, Mark Flory, Hookeun Lee, Ruedi Aebersold,
      and Benno Schwikowski

      @htmlinclude OpenMS_SteinScottImproveScore.parameters

      @ingroup SpectraComparison
  */
  class OPENMS_DLLAPI SteinScottImproveScore :
    public PeakSpectrumCompareFunctor
  {

public:
    /// default constructor
    SteinScottImproveScore();
    /// copy constructor
    SteinScottImproveScore(const SteinScottImproveScore & source);
    /// destructor
    ~SteinScottImproveScore() override;
    /// assignment operator
    SteinScottImproveScore & operator=(const SteinScottImproveScore & source);
    /**
        @brief Similarity pairwise score

        This function return the similarity score of two Spectra based on SteinScott.
    */
    double operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const override;
    /**
        @brief Similarity pairwise score itself

        This function return the similarity score of itself based on SteinScott.
    */
    double operator()(const PeakSpectrum & spec) const override;
    static PeakSpectrumCompareFunctor * create()
    {
      return new SteinScottImproveScore();
    }

    static const String getProductName()
    {
      return "SteinScottImproveScore";
    }

  };
}


#endif /*OPENMS_COMPARISON_SPECTRA_STEINSCOTTIMPROVESCORE_H*/
