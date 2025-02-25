// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MobilityPeak1D.h>

namespace OpenMS
{
  std::ostream & operator<<(std::ostream & os, const MobilityPeak1D & point)
  {
    os << "POS: " << point.getMobility() << " INT: " << point.getIntensity();
    return os;
  }
}
