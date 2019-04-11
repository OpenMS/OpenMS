// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#pragma once

#include <OpenMS/OpenMSConfig.h>

#include <vector>
#include <map>

namespace OpenMS
{
  /**
    @brief cubic spline interpolation
    as described in R.L. Burden, J.D. Faires, Numerical Analysis, 4th ed.
    PWS-Kent, 1989, ISBN 0-53491-585-X, pp. 126-131.
    
    Construction of the spline takes by far the most time. Evaluating it is rather fast 
    (one evaluation is about 50x faster than construction -- depending on number of points etc.).
   
   */
  class OPENMS_DLLAPI CubicSpline2d
  {

    std::vector<double> a_; ///< constant spline coefficients
    std::vector<double> b_; ///< linear spline coefficients
    std::vector<double> c_; ///< quadratic spline coefficients
    std::vector<double> d_; ///< cubic spline coefficients
    std::vector<double> x_; ///< knots

public:

    /**
     * @brief constructor of spline interpolation
     *
     * The coordinates must match by index. Both vectors must be
     * the same size and sorted in x. Sortedness in x is required
     * for @see SplinePackage.
     *
     * @param x x-coordinates of input data points (knots)
     * @param y y-coordinates of input data points
     */
    CubicSpline2d(const std::vector<double>& x, const std::vector<double>& y);

    /**
     * @brief constructor of spline interpolation
     *
     * @param m (x,y) coordinates of input data points
     */
    CubicSpline2d(const std::map<double, double>& m);

    /**
     * @brief evaluates the spline at position x
     *
     * @param x x-position
     */
    double eval(double x) const;

    /**
     * @brief evaluates derivative of spline at position x
     *
     * @param x x-position
     * @param order order of the derivative
     * Only order 1 or 2 make sense for cubic splines.
     */
    double derivatives(double x, unsigned order) const;

private:

    /**
     * @brief initialize the spline
     *
     * @param x x-coordinates of input data points (knots)
     * @param y y-coordinates of input data points
     */
    void init_(const std::vector<double>& x, const std::vector<double>& y);

  };

}

