// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/COMPARISON/BinnedSpectrumCompareFunctor.h>

#include <cfloat>
#include <cmath>

namespace OpenMS
{

  /**
    @brief Compare functor scoring the shared peaks for similarity measurement

    The details of the score can be found in:
    K. Wan, I. Vidavsky, and M. Gross. Comparing similar spectra: from
    similarity index to spectral contrast angle. Journal of the American Society
    for Mass Spectrometry, 13(1):85{88, January 2002.

    @htmlinclude OpenMS_BinnedSharedPeakCount.parameters

    @see BinnedSpectrumCompareFunctor @see BinnedSpectrum

    @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI BinnedSharedPeakCount :
    public BinnedSpectrumCompareFunctor
  {
public:

    /// default constructor
    BinnedSharedPeakCount();

    /// copy constructor
    BinnedSharedPeakCount(const BinnedSharedPeakCount& source);

    /// destructor
    ~BinnedSharedPeakCount() override;

    /// assignment operator
    BinnedSharedPeakCount& operator=(const BinnedSharedPeakCount& source);

    /** function call operator, calculates the similarity of the given arguments

      @param spec1 First spectrum given as a binned representation
      @param spec2 Second spectrum given as a binned representation
      @throw IncompatibleBinning is thrown if the binning of the two input spectra are not the same
    */
    double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const override;

    /// function call operator, calculates self similarity
    double operator()(const BinnedSpectrum& spec) const override;

protected:
    void updateMembers_() override;
    double precursor_mass_tolerance_;
  };

}
