// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/MarkerMower.h>

using namespace std;

namespace OpenMS
{
  MarkerMower::MarkerMower() :
    DefaultParamHandler("MarkerMower")
  {
  }

  MarkerMower::~MarkerMower() = default;

  MarkerMower::MarkerMower(const MarkerMower & source) :
    DefaultParamHandler(source)
  {
  }

  MarkerMower & MarkerMower::operator=(const MarkerMower & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void MarkerMower::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void MarkerMower::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

  ///@todo violates DefaultParamHandler interface (Andreas)
  void MarkerMower::insertmarker(PeakMarker * pm)
  {
    markers_.push_back(pm);
  }

}
