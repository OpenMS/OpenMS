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

#include <OpenMS/CONCEPT/Macros.h>
#include <vector>
#include <algorithm>    // std::min, std::max
#include <functional>

namespace OpenMS
{

  /**
    @brief LOWESS (locally weighted scatterplot smoothing).

    A non-parametric smoothing technique that fits a simple linear regression
    model to localized subsets of the data, point by point. This is often used
    for retention time alignments.

    The implementation here is optimized for speed and many datapoints.  Note
    that it performs a linear fit, it does not implement quadratic fits. It is
    based on the initial FORTRAN code by W. S. Cleveland published at NETLIB. 

    Note that this should work best for for large datasets with mostly linear
    behavior. For small datasets with non-linear behavior, use the
    LowessSmoothing class.

    @ingroup SignalProcessing
  */
  namespace FastLowessSmoothing
  {

    /**
      @brief Computes a lowess smoothing fit on the input vectors

      This is a fast implementation of a lowess fit that is based on the original
      Fortran code by W. S. Cleveland and it uses some optimizations. 

      @param x The input vector in the first dimension
      @param y The input vector in the second dimension
      @param f Fraction of datapoints to use for each local regression (the span, recommended value: 2/3) 
      @param nsteps The number of robustifying iterations (recommended value: 3)
      @param delta nonnegative parameter which may be used to save computations (recommended value: 0.01 * range of x)
      @param result Result of fit

      \pre The size of the vectors x and y needs to be equal
      \pre The vector needs to have at least 2 elements
      \pre The vector x needs to be sorted
      \pre The f value needs to be between 0 and 1
      \pre The nsteps parameter needs to be zero or larger
      \pre The delta parameter needs to be zero or larger

      The delta parameter allows the algorithm to not perform the regression at
      every data point, as it assumes that points that are close to each other
      will have the same regression parameters. A linear interpolation is used
      to fill in the skipped points, larger values lead to increased speed up.

      The f parameter allows the caller to influence the smoothness. A larger
      values will increase smoothness (recommended value: 2/3) It is the
      fraction of points used to compute each fitted value. Choosing F in the
      range .2 to .8 usually results in a good fit

      The nsteps parameter controls how many iterations are performed in the
      robust fit (setting it to zero turns off the robust fit and the nonrobust
      fit is returned). A value of 2 or 3 should be sufficient for most purposes.

    */
    int OPENMS_DLLAPI lowess(const std::vector<double>& x, const std::vector<double>& y,
               double f, int nsteps, double delta, std::vector<double>& result);

    /**
      @brief Computes a lowess smoothing fit on the input vectors with the recommended values

      @param x The input vector in the first dimension
      @param y The input vector in the second dimension
      @param result Result of fit

      \pre The size of the vectors x and y needs to be equal
      \pre The vector needs to have at least 2 elements
      \pre The vector x needs to be sorted
      
    */
    inline int OPENMS_DLLAPI lowess(const std::vector<double>& x, const std::vector<double>& y,
               std::vector<double>& result)
    {
      OPENMS_PRECONDITION(x.size() == y.size(), "Vectors x and y must have the same length")
      OPENMS_PRECONDITION(x.size() >= 2, "Need at least two points for smoothing")
      OPENMS_PRECONDITION(std::adjacent_find(x.begin(), x.end(), std::greater<double>()) == x.end(),
          "The vector x needs to be sorted")

      double delta = 0.01 * (x[ x.size()-1 ] - x[0]); // x is sorted
      return lowess(x, y, 2.0/3, 3, delta, result);
    }
  }

} // namespace OpenMS
