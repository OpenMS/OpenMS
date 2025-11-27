// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/BinnedSpectrumCompareFunctor.h>
#include <OpenMS/COMPARISON/BinnedSharedPeakCount.h>
#include <OpenMS/COMPARISON/BinnedSpectralContrastAngle.h>
#include <OpenMS/COMPARISON/BinnedSumAgreeingIntensities.h>

using namespace std;

namespace OpenMS
{
  BinnedSpectrumCompareFunctor::BinnedSpectrumCompareFunctor() :
    DefaultParamHandler("BinnedSpectrumCompareFunctor")
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
