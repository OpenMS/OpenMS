// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>

using namespace std;

namespace OpenMS
{
  SqrtMower::SqrtMower() :
    DefaultParamHandler("SqrtMower")
  {
  }

  SqrtMower::~SqrtMower() = default;

  SqrtMower::SqrtMower(const SqrtMower & source) = default;

  SqrtMower & SqrtMower::operator=(const SqrtMower & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void SqrtMower::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void SqrtMower::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
