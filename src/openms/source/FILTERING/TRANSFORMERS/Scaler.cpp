// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>

using namespace std;
namespace OpenMS
{
  Scaler::Scaler() :
    DefaultParamHandler("Scaler")
  {
  }

  Scaler::~Scaler() = default;

  Scaler::Scaler(const Scaler & source) = default;

  Scaler & Scaler::operator=(const Scaler & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void Scaler::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void Scaler::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
