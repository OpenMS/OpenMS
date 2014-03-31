// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#include <OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h>

#include <algorithm>
#include <cstdlib>

#include <OpenMS/MATH/STATISTICS/QuadraticRegression.h>


namespace OpenMS
{
  LowessSmoothing::LowessSmoothing() :
    DefaultParamHandler("LowessSmoothing")
  {
    defaults_.setValue("window_size", 10, "The number of peaks to be included for local fitting in one window.");
    defaultsToParam_();
  }

  LowessSmoothing::~LowessSmoothing()
  {
  }

  void LowessSmoothing::smoothData(const DoubleVector & input_x, const DoubleVector & input_y, DoubleVector & smoothed_output)
  {
    if (input_x.size() != input_y.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Sizes of x and y values not equal! Aborting... ", String(input_x.size()));
    }

    // unable to smooth over 2 or less data points (we need at least 3)
    if (input_x.size() <= 2)
    {
      smoothed_output = input_y;
      return;
    }

    Size input_size = input_y.size();

    // const Size q = floor( input_size * alpha );
    const Size q = (window_size_ < input_size) ? window_size_ : input_size - 1;

    DoubleVector distances(input_size, 0.0);
    DoubleVector sortedDistances(input_size, 0.0);

    for (Size outer_idx = 0; outer_idx < input_size; ++outer_idx)
    {
      // Compute distances.
      // Size inner_idx = 0;
      for (Size inner_idx = 0; inner_idx < input_size; ++inner_idx)
      {
        distances[inner_idx] = std::fabs(input_x[outer_idx] - input_x[inner_idx]);
        sortedDistances[inner_idx] = distances[inner_idx];
      }

      // Sort distances in order from smallest to largest.
      std::sort(sortedDistances.begin(), sortedDistances.end());

      // Compute weigths.
      std::vector<double> weigths (input_size, 0);
      for (Size inner_idx = 0; inner_idx < input_size; ++inner_idx)
      {
        weigths.at(inner_idx) = tricube_(distances[inner_idx], sortedDistances[q]);
      }

      //calculate regression
      Math::QuadraticRegression qr;
      std::vector<double>::const_iterator w_begin = weigths.begin();
      qr.computeRegressionWeighted(input_x.begin(), input_x.end(), input_y.begin(), w_begin);

      //smooth y-values
      DoubleReal rt = input_x[outer_idx];
      smoothed_output.push_back(qr.eval(rt));
    }

    return;
  }

  DoubleReal LowessSmoothing::tricube_(DoubleReal u, DoubleReal t)
  {
    // In our case, u represents a distance and hence should be strictly positive.
    if (u < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Value of u must be strictly positive! Aborting...", String(u));
    }

    // 0 <= u < t; u is regarded as 0.0 if fabs(u) falls below epsilon
    if ((fabs(u) < std::numeric_limits<double>::epsilon() || (0.0 < u)) && (u < t))
    {
      // (1 - (u/t)^3)^3
      // return pow( ( 1.0 - pow(u/t, 3.0)), 3.0 );
      DoubleReal quot(u / t);
      DoubleReal inner_term(1.0 - quot * quot * quot);

      return inner_term * inner_term * inner_term;
    }
    // u >= t
    else
    {
      return 0.0;
    }
  }

  void LowessSmoothing::updateMembers_()
  {
    window_size_ = (Size)param_.getValue("window_size");
  }

} //namespace OpenMS
