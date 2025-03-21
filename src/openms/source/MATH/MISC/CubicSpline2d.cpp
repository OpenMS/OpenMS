// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <functional>

using namespace std;

namespace OpenMS
{
  CubicSpline2d::CubicSpline2d(const std::vector<double>& x, const std::vector<double>& y)
  {
    if (x.size() != y.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "x and y vectors are not of the same size.");
    }

    if (x.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "x and y vectors need to contain two or more elements.");
    }

    // assert spectrum is sorted
    if (std::adjacent_find(x.begin(), x.end(), std::greater<double>()) != x.end())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "x vector is not sorted.");
    }

    init_(x, y);
  }

  CubicSpline2d::CubicSpline2d(const std::map<double, double>& m)
  {
    if (m.size() < 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Map needs to contain two or more elements.");
    }

    std::vector<double> x;
    std::vector<double> y;

    x.reserve(m.size());
    y.reserve(m.size());

    std::map<double, double>::const_iterator map_it;
    for (map_it = m.begin(); map_it != m.end(); ++map_it)
    {
      x.push_back(map_it->first);
      y.push_back(map_it->second);
    }

    init_(x, y);
  }

  double CubicSpline2d::eval(double x) const
  {
    if (x < x_.front() || x > x_.back())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Argument out of range of spline interpolation.");
    }

    // determine index of closest node left of (or exactly at) x
    unsigned i = static_cast<unsigned>(std::lower_bound(x_.begin(), x_.end(), x) - x_.begin());
    if (x_[i] > x || x_.back() == x)
    {
        --i;
    }
    
    const double xx = x - x_[i];
    return ((d_[i] * xx + c_[i]) * xx + b_[i]) * xx + a_[i];
  }

  double CubicSpline2d::derivative(const double x) const
  {
    return derivatives(x, 1);
  }

  double CubicSpline2d::derivatives(double x, unsigned order) const
  {
    if (x < x_.front() || x > x_.back())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Argument out of range of spline interpolation.");
    }

    if (order < 1 || order > 3)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Only first, second and third derivative defined on cubic spline");
    }

    // determine index of closest node left of (or exactly at) x
    unsigned i = static_cast<unsigned>(std::lower_bound(x_.begin(), x_.end(), x) - x_.begin());
    if (x_[i] > x || x_.back() == x) // also, i must not point to last index in 'x_', since all other vectors are one element shorter
    {
      --i;
    }
    
    const double xx = x - x_[i];
    if (order == 1)
    {
      return b_[i] + 2 * c_[i] * xx + 3 * d_[i] * xx * xx;
    }
    else if (order == 2)
    {
      return 2 * c_[i] + 6 * d_[i] * xx;
    }
    else
    {
      return 6 * d_[i];
    }
  }

  void CubicSpline2d::init_(const std::vector<double>& x, const std::vector<double>& y)
  {
    const size_t n = x.size() - 1;

    std::vector<double> h;
    h.reserve(n);
    a_.reserve(n);
    x_.reserve(n + 1);
    // do the 0'th element manually, since the loop below only starts at 1
    h.push_back(x[1] - x[0]);
    x_.push_back(x[0]);
    a_.push_back(y[0]);

    std::vector<double> mu(n, 0.0);
    std::vector<double> z(n, 0.0);
    for (unsigned i = 1; i < n; ++i)
    {
      h.push_back(x[i + 1] - x[i]);
      const double l = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
      mu[i] = h[i] / l;
      z[i] = (3 * (y[i + 1] * h[i - 1] - y[i] * (x[i + 1] - x[i - 1]) + y[i - 1] * h[i]) / (h[i - 1] * h[i]) - h[i - 1] * z[i - 1]) / l;
      // store x,y -- required for evaluation later on
      x_.push_back(x[i]);
      a_.push_back(y[i]);
    }
    // 'x_' needs to be full length (all other member vectors (except c_) are one element shorter)
    x_.push_back(x[n]);

    b_.resize(n); 
    d_.resize(n);
    c_.resize(n+1);
    c_.back() = 0;
    for (int j = static_cast<int>(n) - 1; j >= 0; --j)
    {
      c_[j] = z[j] - mu[j] * c_[j + 1];
      b_[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c_[j + 1] + 2 * c_[j]) / 3;
      d_[j] = (c_[j + 1] - c_[j]) / (3 * h[j]);
    }

  }

}
