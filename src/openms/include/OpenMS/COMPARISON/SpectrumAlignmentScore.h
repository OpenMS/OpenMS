// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h>

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

    // @}

  };

}
