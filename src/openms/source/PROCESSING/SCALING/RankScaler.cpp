// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/PROCESSING/SCALING/RankScaler.h>

using namespace std;
namespace OpenMS
{
  RankScaler::RankScaler() :
    DefaultParamHandler("RankScaler")
  {
  }

  RankScaler::~RankScaler() = default;

  RankScaler::RankScaler(const RankScaler & source) = default;

  RankScaler & RankScaler::operator=(const RankScaler & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void RankScaler::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void RankScaler::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
