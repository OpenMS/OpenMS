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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_PEAKALIGNMENT_H
#define OPENMS_COMPARISON_SPECTRA_PEAKALIGNMENT_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>
#include <vector>

namespace OpenMS
{

  /**
      @brief make a PeakAlignment of two PeakSpectra

      The alignment is done according to the Needleman-Wunsch Algorithm (local alignment considering gaps).

      @htmlinclude OpenMS_PeakAlignment.parameters

      @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI PeakAlignment :
    public PeakSpectrumCompareFunctor
  {
public:

    /// default constructor
    PeakAlignment();

    /// copy constructor
    PeakAlignment(const PeakAlignment & source);

    /// destructor
    ~PeakAlignment() override;

    /// assignment operator
    PeakAlignment & operator=(const PeakAlignment & source);

    /** function call operator, calculates the similarity of the given arguments

        @param spec1 First spectrum given in a binned representation
        @param spec2 Second spectrum given in a binned representation
    */
    double operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const override;

    /// function call operator, calculates self similarity
    double operator()(const PeakSpectrum & spec) const override;

    ///
    static PeakSpectrumCompareFunctor * create() { return new PeakAlignment(); }

    /// get the identifier for this DefaultParamHandler
    static const String getProductName()
    {
      return "PeakAlignment";
    }

    /// make alignment and get the traceback
    std::vector<std::pair<Size, Size> > getAlignmentTraceback(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const;

private:

    /// calculates the score for aligning two peaks
    double peakPairScore_(double & pos1, double & intens1, double & pos2, double & intens2, const double & sigma) const;


  };

}
#endif //OPENMS_COMPARISON_SPECTRA_PEAKALIGNMENT_H
