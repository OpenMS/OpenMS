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

  };

}

