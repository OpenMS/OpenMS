// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h>
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

    /// make alignment and get the traceback
    std::vector<std::pair<Size, Size> > getAlignmentTraceback(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const;

private:

    /// calculates the score for aligning two peaks
    double peakPairScore_(double & pos1, double & intens1, double & pos2, double & intens2, const double & sigma) const;


  };

}
