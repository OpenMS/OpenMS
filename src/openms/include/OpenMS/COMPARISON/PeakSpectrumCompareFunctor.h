// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

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

  };

}
