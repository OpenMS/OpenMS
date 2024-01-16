// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>

namespace OpenMS
{

  /**
    @brief Compare functor scoring the spectral contrast angle for similarity measurement

    The details of the score can be found in:
    K. Wan, I. Vidavsky, and M. Gross. Comparing similar spectra: from
    similarity index to spectral contrast angle. Journal of the American Society
    for Mass Spectrometry, 13(1):85-88, January 2002.

    @htmlinclude OpenMS_BinnedSpectralContrastAngle.parameters

    @see BinnedSpectrumCompareFunctor
    @see BinnedSpectrum

    @ingroup SpectraComparison
  */
  class OPENMS_DLLAPI BinnedSpectralContrastAngle :
    public BinnedSpectrumCompareFunctor
  {

public:

    /// default constructor
    BinnedSpectralContrastAngle();

    /// copy constructor
    BinnedSpectralContrastAngle(const BinnedSpectralContrastAngle& source);

    /// destructor
    ~BinnedSpectralContrastAngle() override;

    /// assignment operator
    BinnedSpectralContrastAngle& operator=(const BinnedSpectralContrastAngle& source);

    /** function call operator, calculates the similarity of the given arguments

      @param spec1 First spectrum given in a binned representation
      @param spec2 Second spectrum given in a binned representation
      @throw IncompatibleBinning is thrown if the bins of the spectra are not the same
    */
    double operator()(const BinnedSpectrum& spec1, const BinnedSpectrum& spec2) const override;

    /// function call operator, calculates self similarity
    double operator()(const BinnedSpectrum& spec) const override;

    ///
    static BinnedSpectrumCompareFunctor* create() { return new BinnedSpectralContrastAngle(); }

    /// get the identifier for this DefaultParamHandler
    static const String getProductName()
    {
      return "BinnedSpectralContrastAngle";
    }

protected:
    void updateMembers_() override;
    double precursor_mass_tolerance_;
  };

}
