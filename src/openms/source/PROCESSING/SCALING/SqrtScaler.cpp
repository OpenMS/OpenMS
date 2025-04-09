// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/PROCESSING/SCALING/SqrtScaler.h>

using namespace std;

namespace OpenMS
{
  SqrtScaler::SqrtScaler() :
    DefaultParamHandler("SqrtScaler")
  {
  }

  SqrtScaler::~SqrtScaler() = default;

  SqrtScaler::SqrtScaler(const SqrtScaler & source) = default;

  SqrtScaler & SqrtScaler::operator=(const SqrtScaler & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void SqrtScaler::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void SqrtScaler::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
