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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>

#include <cmath>
#include <cfloat>

namespace OpenMS
{

  /**
    @brief Sum of agreeing intensities for similarity measurement

    Transformation and other factors of the peptide mass spectrometry pairwise peak-list comparison process
    Witold E Wolski , Maciej Lalowski* , Peter Martus* , Ralf Herwig* , Patrick Giavalisco , Johan Gobom , Albert Sickmann , Hans Lehrach and Knut Reinert*
    BMC Bioinformatics 2005, 6:285 doi:10.1186/1471-2105-6-285

    Bins whose intensity differences are larger than their average intensity receive a weight of zero.

    Prefect agreement results in a similarity score of 1.0

    @htmlinclude OpenMS_BinnedSumAgreeingIntensities.parameters

    @see BinnedSpectrumCompareFunctor
    @see BinnedSpectrum

    @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI BinnedSumAgreeingIntensities :
    public BinnedSpectrumCompareFunctor
  {
public:

    /// default constructor
    BinnedSumAgreeingIntensities();

    /// copy constructor
    BinnedSumAgreeingIntensities(const BinnedSumAgreeingIntensities& source);

    /// destructor
    ~BinnedSumAgreeingIntensities() override;

    /// assignment operator
    BinnedSumAgreeingIntensities& operator=(const BinnedSumAgreeingIntensities& source);

    /** function call operator, calculates the similarity of the given arguments

      @param spec1 First spectrum given as a binned representation
      @param spec2 Second spectrum given as a binned representation
      @throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same
    */
    double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const override;

    /// function call operator, calculates self similarity
    double operator()(const BinnedSpectrum& spec) const override;

    ///
    static BinnedSpectrumCompareFunctor* create() { return new BinnedSumAgreeingIntensities(); }

    /// get the identifier for this DefaultParamHandler
    static const String getProductName()
    {
      return "BinnedSumAgreeingIntensities";
    }

protected:
    void updateMembers_() override;
    double precursor_mass_tolerance_;
  };

}
