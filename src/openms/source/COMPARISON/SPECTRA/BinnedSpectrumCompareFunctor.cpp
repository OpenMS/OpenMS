// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectralContrastAngle.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSumAgreeingIntensities.h>
#include <OpenMS/CONCEPT/Factory.h>



using namespace std;

namespace OpenMS
{
  BinnedSpectrumCompareFunctor::BinnedSpectrumCompareFunctor() :
    DefaultParamHandler(BinnedSpectrumCompareFunctor::getProductName())
  {
  }

  BinnedSpectrumCompareFunctor::BinnedSpectrumCompareFunctor(const BinnedSpectrumCompareFunctor & source) = default;

  BinnedSpectrumCompareFunctor::~BinnedSpectrumCompareFunctor() = default;

  BinnedSpectrumCompareFunctor & BinnedSpectrumCompareFunctor::operator=(const BinnedSpectrumCompareFunctor & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

}
