// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

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

