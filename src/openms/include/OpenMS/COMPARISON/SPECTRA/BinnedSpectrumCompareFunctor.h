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
#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUMCOMPAREFUNCTOR_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUMCOMPAREFUNCTOR_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>

#include <cmath>

namespace OpenMS
{

  /**
    @brief Base class for compare functors of BinnedSpectra

    BinnedSpectrumCompareFunctor classes return a value for a pair of BinnedSpectrum objects (or a single one with itself).
    Ideally the value should reflect the similarity of the pair. For methods of computing the similarity see the
    documentation of the concrete functors.
    Functors normalized in the range [0,1] are identifiable at the set "normalized" parameter of the ParameterHandler

    @ingroup SpectraComparison
  */
  class OPENMS_DLLAPI BinnedSpectrumCompareFunctor :
    public DefaultParamHandler
  {

private:

public:

    /**
      @brief Exception thrown if compared spectra are incompatible

      the compared spectra have different settings in binsize and/or binspread
      due to which comparison would fail
    */
    class OPENMS_DLLAPI IncompatibleBinning :
      public Exception::BaseException
    {
public:
      IncompatibleBinning(const char* file, int line, const char* function, const char* message
                            = "compared spectra have different settings in binsize and/or binspread")  throw();
      ~IncompatibleBinning() throw() override;
    };

    /// default constructor
    BinnedSpectrumCompareFunctor();

    /// copy constructor
    BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor& source);

    /// destructor
    ~BinnedSpectrumCompareFunctor() override;

    /// assignment operator
    BinnedSpectrumCompareFunctor& operator=(const BinnedSpectrumCompareFunctor& source);

    /// function call operator, calculates the similarity of the given arguments
    virtual double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const = 0;

    /// function call operator, calculates self similarity
    virtual double operator()(const BinnedSpectrum& spec) const = 0;

    /// registers all derived products
    static void registerChildren();

    /// get the identifier for a DefaultParamHandler
    static const String getProductName()
    {
      return "BinnedSpectrumCompareFunctor";
    }

  };

}

#endif // OPENMS_COMPARISON_SPECTRA_BINNEDSPECTRUMCOMPAREFUNCTOR_H
