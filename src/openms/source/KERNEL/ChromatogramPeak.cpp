// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ChromatogramPeak.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  std::ostream & operator<<(std::ostream & os, const ChromatogramPeak & point)
  {
    os << "POS: " << point.getRT() << " INT: " << point.getIntensity();
    return os;
  }
}
