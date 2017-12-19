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
#ifndef OPENMS_COMPARISON_SPECTRA_PEAKSPECTRUMCOMPAREFUNCTOR_H
#define OPENMS_COMPARISON_SPECTRA_PEAKSPECTRUMCOMPAREFUNCTOR_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{

  /**

      @brief Base class for compare functors of spectra, that return a similarity value for two spectra.

  PeakSpectrumCompareFunctor classes return a similarity value for a pair of PeakSpectrum objects.
  The value should be greater equal 0.

      @ingroup SpectraComparison
*/
  class OPENMS_DLLAPI PeakSpectrumCompareFunctor :
    public DefaultParamHandler
  {

public:

    /// default constructor
    PeakSpectrumCompareFunctor();

    /// copy constructor
    PeakSpectrumCompareFunctor(const PeakSpectrumCompareFunctor & source);

    /// destructor
    ~PeakSpectrumCompareFunctor() override;

    /// assignment operator
    PeakSpectrumCompareFunctor & operator=(const PeakSpectrumCompareFunctor & source);

    /// function call operator, calculates the similarity
    virtual double operator()(const PeakSpectrum & a, const PeakSpectrum & b) const = 0;

    /// calculates self similarity
    virtual double operator()(const PeakSpectrum & a) const = 0;

    /// registers all derived products
    static void registerChildren();

    ///
    static const String getProductName()
    {
      return "PeakSpectrumCompareFunctor";
    }

  };

}
#endif // OPENMS_COMPARISON_SPECTRA_COMPAREFUNCTOR_H
