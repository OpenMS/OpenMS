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

protected:
    void updateMembers_() override;
    double precursor_mass_tolerance_;
  };

}
