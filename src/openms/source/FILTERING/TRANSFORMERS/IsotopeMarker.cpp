// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>

using namespace std;

namespace OpenMS
{

  IsotopeMarker::IsotopeMarker() :
    PeakMarker()
  {
    setName(IsotopeMarker::getProductName());
    defaults_.setValue("marks", 1, "How often a peak must be marked to be reported");
    defaults_.setValue("mz_variation", 0.1, "variation in m/z direction");
    defaults_.setValue("in_variation", 0.5, "variation in intensity");
    defaultsToParam_();
  }

  IsotopeMarker::IsotopeMarker(const IsotopeMarker & source) = default;

  IsotopeMarker::~IsotopeMarker() = default;

  IsotopeMarker & IsotopeMarker::operator=(const IsotopeMarker & source)
  {
    if (this != &source)
    {
      PeakMarker::operator=(source);
    }
    return *this;
  }

}
