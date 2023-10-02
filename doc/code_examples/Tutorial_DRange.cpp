// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <iostream>

using namespace OpenMS;

Int main()
{
  DRange<2> range;
  range.setMin(DPosition<2>(2.0, 3.0));
  range.setMax(DPosition<2>(1.0, 5.0));

  for (UInt i = 0; i < DRange<2>::DIMENSION; ++i)
  {
    std::cout << "min " << i << ": " << range.minPosition()[i] << std::endl;
    std::cout << "max " << i << ": " << range.maxPosition()[i] << std::endl;
  }

  return 0;
} //end of main
