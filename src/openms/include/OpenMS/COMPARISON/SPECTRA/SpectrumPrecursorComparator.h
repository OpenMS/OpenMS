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
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMPRECURSORCOMPARATOR_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMPRECURSORCOMPARATOR_H

#include <OpenMS/COMPARISON/SPECTRA/PeakSpectrumCompareFunctor.h>

namespace OpenMS
{

  /**
    @brief SpectrumPrecursorComparator compares just the parent mass of two spectra

        @htmlinclude OpenMS_SpectrumPrecursorComparator.parameters

        @ingroup SpectraComparison
  */
  class OPENMS_DLLAPI SpectrumPrecursorComparator :
    public PeakSpectrumCompareFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    SpectrumPrecursorComparator();

    /// copy constructor
    SpectrumPrecursorComparator(const SpectrumPrecursorComparator & source);

    /// destructor
    ~SpectrumPrecursorComparator() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    SpectrumPrecursorComparator & operator=(const SpectrumPrecursorComparator & source);

    double operator()(const PeakSpectrum & a, const PeakSpectrum & b) const override;

    double operator()(const PeakSpectrum & a) const override;
    // @}

    // @name Accessors
    // @{
    ///
    static PeakSpectrumCompareFunctor * create() { return new SpectrumPrecursorComparator(); }

    ///
    static const String getProductName()
    {
      return "SpectrumPrecursorComparator";
    }

    // @}

  };

}

#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMPRECURSORCOMPARATOR_H
