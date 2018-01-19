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
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_ZHANGSIMILARITYSCORE_H
#define OPENMS_COMPARISON_SPECTRA_ZHANGSIMILARITYSCORE_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

namespace OpenMS
{

  /**
      @brief Similarity score of Zhang

        The details of the score can be found in:
        Z. Zhang, Prediction of Low-Energy Collision-Induced Dissociation Spectra of Peptides,
        Anal. Chem., 76 (14), 3908 - 3922, 2004

        @htmlinclude OpenMS_ZhangSimilarityScore.parameters

        @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI ZhangSimilarityScore :
    public PeakSpectrumCompareFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    ZhangSimilarityScore();

    /// copy constructor
    ZhangSimilarityScore(const ZhangSimilarityScore & source);

    /// destructor
    ~ZhangSimilarityScore() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    ZhangSimilarityScore & operator=(const ZhangSimilarityScore & source);

    ///
    double operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const override;

    double operator()(const PeakSpectrum & spec) const override;
    // @}

    // @name Accessors
    // @{
    ///
    static PeakSpectrumCompareFunctor * create() { return new ZhangSimilarityScore(); }

    ///
    static const String getProductName()
    {
      return "ZhangSimilarityScore";
    }

    // @}

protected:

    /// returns the factor associated with the m/z tolerance and m/z difference of the peaks
    double getFactor_(double mz_tolerance, double mz_difference, bool is_gaussian = false) const;


  };

}
#endif //OPENMS_COMPARISON_SPECTRA_ZHANGSIMILARTIYSCORE_H
