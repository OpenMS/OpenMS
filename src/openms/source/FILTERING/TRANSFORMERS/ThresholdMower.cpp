// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>

using namespace std;
namespace OpenMS
{

  ThresholdMower::ThresholdMower() :
    DefaultParamHandler("ThresholdMower")
  {
    defaults_.setValue("threshold", 0.05, "Intensity threshold, peaks below this threshold are discarded");
    defaultsToParam_();
  }

  ThresholdMower::~ThresholdMower() = default;

  ThresholdMower::ThresholdMower(const ThresholdMower & source) :
    DefaultParamHandler(source)
  {
  }

  ThresholdMower & ThresholdMower::operator=(const ThresholdMower & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void ThresholdMower::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void ThresholdMower::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
