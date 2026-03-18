// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h>

namespace OpenMS
{

  /**
      @brief Similarity score of Zhang

        The details of the score can be found in:
        Z. Zhang, Prediction of Low-Energy Collision-Induced Dissociation Spectra of Peptides,
        Anal. Chem., 76 (14), 3908 - 3922, 2004

        @htmlinclude OpenMS_ZhangSimilarityScore.parameters

        @ingroup SpectraComparison
  */

  class OPENMS_DLLAPI ZhangSimilarityScore :
    public PeakSpectrumCompareFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    ZhangSimilarityScore();

    /// copy constructor
    ZhangSimilarityScore(const ZhangSimilarityScore & source);

    /// destructor
    ~ZhangSimilarityScore() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    ZhangSimilarityScore & operator=(const ZhangSimilarityScore & source);

    ///
    double operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const override;

    double operator()(const PeakSpectrum & spec) const override;
    // @}

protected:

    /// returns the factor associated with the m/z tolerance and m/z difference of the peaks
    double getFactor_(double mz_tolerance, double mz_difference, bool is_gaussian = false) const;


  };

}
