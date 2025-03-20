// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


// STL
#include <iostream>
#include <cmath>
#include <numeric>

// OpenMS
#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/MATH/MathFunctions.h>

using namespace std;

namespace OpenMS
{
  using namespace Math;

  void AxisTickCalculator::calcGridLines(double x1, double x2, GridVector & grid)
  {
    grid.clear();

    if (std::isnan(x1) || std::isnan(x2))
    {
      return;
    }
    if (x1 > -0.0001 && x1 < 0.0001)
    {
      x1 = 0.0001;
    }

    double dx = x2 - x1;

    if (dx < 0.0000001)
    {
      return;
    }
    double epsilon = dx / 200;

    double sDecPow = floor(log10(dx));
    double sDec = pow(10.0, sDecPow);

    UInt n_max_big_gridlines = (UInt)floor(dx / sDec);

    std::vector<double> big;
    double currGL = ceilDecimal(x1, (UInt)sDecPow);

    // big grid lines
    while (currGL < (x2 + epsilon))
    {
      big.push_back(currGL);
      currGL += sDec;
    }
    grid.push_back(big);

    // only one major grid line
    if (n_max_big_gridlines <= 3)
    {
      std::vector<double> small;
      currGL = grid[0][0] - sDec * 9 / 10;
      while (currGL < (x2 + epsilon))
      {
        // check if in visible area
        if (currGL > x1)
        {
          // check for overlap with big gridlines
          bool overlap = false;
          for (Size i = 0; i != big.size(); ++i)
          {
            if (fabs(big[i] - currGL) < epsilon)
            {
              overlap = true;
            }
          }
          if (!overlap)
          {
            small.push_back(currGL);
          }
        }
        currGL += sDec / 10;
      }
      grid.push_back(small);
      return;
    }

    // four or more major grid lines
    std::vector<double> small;
    currGL = grid[0][0] - sDec / 2;
    while (currGL < (x2 + epsilon))
    {
      if (currGL > x1)
      {
        small.push_back(currGL);
      }
      currGL += sDec;
    }

    grid.push_back(small);
  }

  void AxisTickCalculator::calcLogGridLines(double x1, double x2, GridVector& grid)
  {
    if (std::isnan(x1)) {
      x1 = 0; // may happen for negative values
    }
    if (std::isnan(x2))
    {
      x2 = 0; // may happen for negative values
    }
    double dx = x2 - x1;
    if (dx < 0.00000001)
    {
      return;
    }

    // small ticks covering one log order
    static const double TICK_VALUES[] = {log10(2.0), log10(3.0), log10(4.0), log10(5.0), log10(6.0), log10(7.0), log10(8.0), log10(9.0)};
    static constexpr int TICK_COUNT = sizeof(TICK_VALUES) / sizeof(decltype(TICK_VALUES[0]));

    grid.clear();
    grid.resize(2);
    Int x1floor = (Int)floor(x1);
    Int x2ceil = (Int)ceil(x2);
    std::vector<double>& big = grid[0];
    std::vector<double>& small = grid[1];
    big.resize(x2ceil - x1floor);
    std::iota(big.begin(), big.end(), x1floor);

    for (const auto b : big)
    {
      for (Int j = 0; j != TICK_COUNT; ++j)
      {
        const auto mini_tick = b + TICK_VALUES[j];
        if (mini_tick > x2)
        {
          break;
        }
        small.push_back(mini_tick);
      }
    }
  }

} //namespace
