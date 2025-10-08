// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#include <OpenMS/PROCESSING/SMOOTHING/LowessSmoothing.h>
#include <OpenMS/ML/REGRESSION/QuadraticRegression.h>

#include <algorithm>
#include <cmath>

namespace OpenMS
{
  LowessSmoothing::LowessSmoothing() :
    DefaultParamHandler("LowessSmoothing")
  {
    defaults_.setValue("window_size", 10, "The number of peaks to be included for local fitting in one window.");
    defaultsToParam_();
  }

  LowessSmoothing::~LowessSmoothing() = default;

  void LowessSmoothing::smoothData(const DoubleVector& input_x, const DoubleVector& input_y, DoubleVector& smoothed_output)
  {
    if (input_x.size() != input_y.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Sizes of x and y values not equal! Aborting... ", String(input_x.size()));
    }

    // unable to smooth over 2 or less data points (we need at least 3)
    if (input_x.size() <= 2)
    {
      smoothed_output = input_y;
      return;
    }

    Size input_size = input_y.size();

    // const Size q = floor( input_size * alpha );
    const Size q = (window_size_ < input_size) ? static_cast<Size>(window_size_) : input_size - 1;

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
      std::vector<double> weigths(input_size, 0);
      for (Size inner_idx = 0; inner_idx < input_size; ++inner_idx)
      {
        weigths.at(inner_idx) = tricube_(distances[inner_idx], sortedDistances[q]);
      }

      //calculate regression
      Math::QuadraticRegression qr;
      std::vector<double>::const_iterator w_begin = weigths.begin();
      qr.computeRegressionWeighted(input_x.begin(), input_x.end(), input_y.begin(), w_begin);

      //smooth y-values
      double rt = input_x[outer_idx];
      smoothed_output.push_back(qr.eval(rt));
    }

    return;
  }

  double LowessSmoothing::tricube_(double u, double t)
  {
    // In our case, u represents a distance and hence should be strictly positive.
    if (u < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Value of u must be strictly positive! Aborting...", String(u));
    }

    // 0 <= u < t; u is regarded as 0.0 if fabs(u) falls below epsilon
    if ((std::fabs(u) < std::numeric_limits<double>::epsilon() || (0.0 < u)) && (u < t))
    {
      // (1 - (u/t)^3)^3
      // return pow( ( 1.0 - pow(u/t, 3.0)), 3.0 );
      double quot(u / t);
      double inner_term(1.0 - quot * quot * quot);

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
