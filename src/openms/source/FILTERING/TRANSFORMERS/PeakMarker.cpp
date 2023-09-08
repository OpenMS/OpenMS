// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>
namespace OpenMS
{
  PeakMarker::PeakMarker() :
    DefaultParamHandler(PeakMarker::getProductName())
  {
  }

  PeakMarker::PeakMarker(const PeakMarker & source) = default;

  PeakMarker::~PeakMarker() = default;

  PeakMarker & PeakMarker::operator=(const PeakMarker & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

}
