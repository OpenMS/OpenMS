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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>

#include <vector>
#include <map>

using namespace std;

namespace OpenMS
{
  CubicSpline2d::CubicSpline2d(const std::vector<double>& x, const std::vector<double>& y)
  {
    if (x.empty() || x.size() != y.size())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "x and y vectors either not of the same size or empty.");
    }

    init_(x, y);
  }

  CubicSpline2d::CubicSpline2d(const std::map<double, double>& m)
  {
    if (m.empty())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Map is empty.");
    }

    std::vector<double> x;
    std::vector<double> y;

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
    if (x < x_[0] || x > x_[x_.size() - 1])
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Argument out of range of spline interpolation.");
    }

    int i = std::lower_bound(x_.begin(), x_.end(), x) - x_.begin() - 1;
    double xx = x - x_[i];

    return ((d_[i] * xx + c_[i]) * xx + b_[i]) * xx + a_[i];
  }

  double CubicSpline2d::derivatives(double x, unsigned order) const
  {
    if (x < x_[0] || x > x_[x_.size() - 1])
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Argument out of range of spline interpolation.");
    }

    if (order != 1 && order != 2)
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Only first and second derivative defined on cubic spline");
    }

    int i = std::lower_bound(x_.begin(), x_.end(), x) - x_.begin() - 1;
    double xx = x - x_[i];

    if (order == 1)
    {
      return (d_[i] * xx + c_[i]) * xx + b_[i];
    }
    else
    {
      return d_[i] * xx + c_[i];
    }
  }

  void CubicSpline2d::init_(const std::vector<double>& x, const std::vector<double>& y)
  {
    const size_t n = x.size() - 1;

    std::vector<double> h(n, 0.0);
    for (unsigned i = 0; i < n; ++i)
    {
      h[i] = x[i + 1] - x[i];
    }

    std::vector<double> mu(n, 0.0);
    std::vector<double> z(n + 1, 0.0);
    for (unsigned i = 1; i < n; ++i)
    {
      double l = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
      mu[i] = h[i] / l;
      z[i] = (3 * (y[i + 1] * h[i - 1] - y[i] * (x[i + 1] - x[i - 1]) + y[i - 1] * h[i]) / (h[i - 1] * h[i]) - h[i - 1] * z[i - 1]) / l;
    }

    std::vector<double> b(n, 0.0);
    std::vector<double> c(n + 1, 0.0);
    std::vector<double> d(n, 0.0);
    for (int j = n - 1; j >= 0; --j)
    {
      c[j] = z[j] - mu[j] * c[j + 1];
      b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3;
      d[j] = (c[j + 1] - c[j]) / (3 * h[j]);
    }

    for (unsigned i = 0; i < n; ++i)
    {
      a_.push_back(y[i]);
      b_.push_back(b[i]);
      c_.push_back(c[i]);
      d_.push_back(d[i]);
      x_.push_back(x[i]);
    }
    x_.push_back(x[n]);
  }

}
