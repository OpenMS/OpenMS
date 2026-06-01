// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <iostream>

using namespace OpenMS;

Int main()
{
  DPosition<2> pos {-8.15, 47.11};
  static_assert(pos.size() == 2);

  std::cout << "largest possible value: " << DPosition<2>::maxPositive() << '\n';
  // make values in all dimensions positive and print
  std::cout << "abs: " << pos.abs() << '\n';

  // manipulate individual dimensions
  pos[0] = -3.15;
  pos[1] = 7.11;

  for (Size i = 0; i < pos.DIMENSION; ++i)
  {
    std::cout << "Dimension " << i << ": " << pos[i] << std::endl;
  }
  // same thing
  int i = 0;
  for (const auto e : pos)
  {
    std::cout << "Dimension " << i++ << ": " << e << std::endl;
  }

  return 0;
} //end of main
