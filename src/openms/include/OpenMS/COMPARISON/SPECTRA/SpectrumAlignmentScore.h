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
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENTSCORE_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENTSCORE_H

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

#include <boost/math/special_functions/erf.hpp>
#include <cmath>

namespace OpenMS
{

  /**
      @brief Similarity score via spectra alignment

        This class implements a simple scoring based on the alignment of spectra. This alignment
        is implemented in the SpectrumAlignment class and performs a dynamic programming alignment
        of the peaks, minimizing the distances between the aligned peaks and maximizing the number
        of peak pairs.

        The scoring is done via the simple formula score = sum / (sqrt(sum1 * sum2)). sum is the
        product of the intensities of the aligned peaks, with the given exponent (default is 2).
        sum1 and sum2 are the sum of the intensities squared for each peak of both spectra respectively.

        A binned version of this scoring is implemented in the ZhangSimilarityScoring class.

        @htmlinclude OpenMS_SpectrumAlignmentScore.parameters

        @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI SpectrumAlignmentScore :
    public PeakSpectrumCompareFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    SpectrumAlignmentScore();

    /// copy constructor
    SpectrumAlignmentScore(const SpectrumAlignmentScore & source);

    /// destructor
    ~SpectrumAlignmentScore() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    SpectrumAlignmentScore & operator=(const SpectrumAlignmentScore & source);

    ///
    double operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const override;

    double operator()(const PeakSpectrum & spec) const override;
    // @}

    // @name Accessors
    // @{
    ///
    static PeakSpectrumCompareFunctor * create() { return new SpectrumAlignmentScore(); }

    ///
    static const String getProductName()
    {
      return "SpectrumAlignmentScore";
    }

    // @}

protected:

    /// returns the factor associated with the m/z tolerance and m/z difference of the peaks
    double getFactor_(double mz_tolerance, double mz_difference, bool is_gaussian = false) const;
  };

}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENTSCORE_H
