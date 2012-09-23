// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

// OpenMS
#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

namespace OpenMS
{
  using namespace Math;

  void AxisTickCalculator::calcGridLines(DoubleReal x1, DoubleReal x2, GridVector & grid)
  {
    grid.clear();

    if (boost::math::isnan(x1) || boost::math::isnan(x2))
      return;

    if (x1 > -0.0001 && x1 < 0.0001)
    {
      x1 = 0.0001;
    }

    DoubleReal dx = x2 - x1;

    if (dx < 0.0000001)
    {
      return;
    }
    DoubleReal epsilon = dx / 200;

    DoubleReal sDecPow = floor(log10(dx));
    DoubleReal sDec = pow(10.0, sDecPow);

    UInt n_max_big_gridlines = (UInt)floor(dx / sDec);

    std::vector<DoubleReal> big;
    DoubleReal currGL = ceilDecimal(x1, (UInt)sDecPow);

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
      std::vector<DoubleReal> small;
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
    std::vector<DoubleReal> small;
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

  void AxisTickCalculator::calcLogGridLines(DoubleReal x1, DoubleReal x2, GridVector & grid)
  {
    grid.clear();
    DoubleReal scalValues[8];
    scalValues[0] = log10(2.0);
    scalValues[1] = log10(3.0);
    scalValues[2] = log10(4.0);
    scalValues[3] = log10(5.0);
    scalValues[4] = log10(6.0);
    scalValues[5] = log10(7.0);
    scalValues[6] = log10(8.0);
    scalValues[7] = log10(9.0);
    DoubleReal dx = x2 - x1;

    if (dx < 0.00000001)
    {
      return;
    }

    Int x1ceil = (Int)floor(x1);
    Int x2floor = (Int)ceil(x2);
    std::vector<DoubleReal> big;
    for (Int i = x1ceil; i != x2floor; ++i)
    {
      big.push_back(i);
    }
    grid.push_back(big);
    std::vector<DoubleReal> small;
    for (Size i = 0; i != grid[0].size(); ++i)
    {
      DoubleReal currGL = grid[0][i];
      for (Int j = 0; j != 8; ++j)
      {
        if (currGL + scalValues[j] > x2)
        {
          break;
        }
        small.push_back(currGL + scalValues[j]);
      }
    }
    grid.push_back(small);
  }

} //namespace
