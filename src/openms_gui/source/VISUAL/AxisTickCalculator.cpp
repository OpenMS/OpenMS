// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
#include <OpenMS/MATH/MISC/MathFunctions.h>

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
