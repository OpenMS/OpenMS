// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <iostream>

using namespace OpenMS;

Int main()
{
  // A D-dimensional range, without units;
  // Note: if you want something more modern with dimensions for RT, m/z, intensity and mobility, then use RangeManager.
  // 
  // You can use any dimension you like; for D=2 and D=3 there are some convenience overloads though, especially for C'tors
  // a 2-dimensional, i.e. [x_min..x_max, y_min..y_max],  range
  DRange<2> range;
  range.setMin(DPosition<2>(2.0, 3.0)); // for (x_min, y_min)
  range.setMax(DPosition<2>(4.0, 5.0)); // for (x_max, y_max)
  std::cout << "values:\n" << range; // prints [2..4, 3..5]

  // Note: the class maintains the invariant min<=max for each dimension
  //       Thus, setting a 'min' which is larger than the current 'max', also adjusts 'max' to the same value
  range.setMin(DPosition<2>(10.0, 2.0)); // for (x_max, y_max)
  std::cout << "\nadjusted max:\n" << range; // prints [10..10, 2..5]

  // you can also set each dimension's min/max: 0 = X, 1 = Y
  range.setDimMinMax(0, {0.6, 6.6});
  std::cout << "\nnew X range:\n" << range;

  // print values using a custom format
  for (UInt i = 0; i < range.DIMENSION; ++i)
  {
    std::cout << "DIM " << i << ": " << range.minPosition()[i] << " ... " << range.maxPosition()[i] << '\n';
  }


} //end of main
