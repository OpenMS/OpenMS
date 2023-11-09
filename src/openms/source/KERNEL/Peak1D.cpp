// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/Peak1D.h>

namespace OpenMS
{
  std::ostream & operator<<(std::ostream & os, const Peak1D & point)
  {
    os << "POS: " << point.getMZ() << " INT: " << point.getIntensity();
    return os;
  }

}
