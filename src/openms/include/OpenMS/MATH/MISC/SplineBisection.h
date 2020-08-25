// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <limits> 
#include <cmath> 

namespace OpenMS
{
  /**
   * @brief Uses bisection to find the maximum point of a spline
   *
   * Should work with BSpline2d and CubicSpline2d
   *
   */
  namespace Math
  {

    template <class T>
    void spline_bisection(const T & peak_spline, 
        double const left_neighbor_mz,
        double const right_neighbor_mz,
        double & max_peak_mz,
        double & max_peak_int,
        double const threshold = 1e-6)
    {
      // calculate maximum by evaluating the spline's 1st derivative
      // (bisection method)
      double lefthand = left_neighbor_mz;
      double righthand = right_neighbor_mz;

      bool lefthand_sign = 1;
      double eps = std::numeric_limits<double>::epsilon();

      // bisection
      do
      {
        double mid = (lefthand + righthand) / 2.0;
        double midpoint_deriv_val = peak_spline.derivatives(mid, 1);

        // if deriv nearly zero then maximum already found
        if (!(std::fabs(midpoint_deriv_val) > eps))
        {
          break;
        }

        bool midpoint_sign = (midpoint_deriv_val < 0.0) ? 0 : 1;

        if (lefthand_sign ^ midpoint_sign)
        {
          righthand = mid;
        }
        else
        {
          lefthand = mid;
        }
      }
      while (righthand - lefthand > threshold);

      max_peak_mz = (lefthand + righthand) / 2;
      max_peak_int = peak_spline.eval(max_peak_mz);
    }

  }
}
