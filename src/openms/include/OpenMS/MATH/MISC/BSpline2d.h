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
// $Maintainer: .. $
// $Authors: .. $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_MISC_BSPLINE2D_H
#define OPENMS_MATH_MISC_BSPLINE2D_H

#include <vector>

#include <OpenMS/config.h>

// forward declaration of impl class BSpline
template <class T>
struct BSpline;

namespace OpenMS
{
  /**
   * @brief b spline interpolation
   */
  class OPENMS_DLLAPI BSpline2d
  {
public:
    /**
    */
    BSpline2d(const std::vector<double>& x, const std::vector<double>& y);

    /**
     * Solve the spline curve for a new set of y values.  Returns false
     * if the solution fails.
     *
     * @param y The array of y values corresponding to each of the nX()
     *		x values in the domain.
     */
    bool solve (const std::vector<double> y);

    /**
     * Return the evaluation of the smoothed curve
     * at a particular @p x value.  If current state is not ok(), returns 0.
     */
    double evaluate (double x);

    /**
     * Return the first derivative of the spline curve at the given @p x.
     * Returns zero if the current state is not ok().
     */
    double slope (double x);

    /**
     * Return the @p n-th basis coefficient, from 0 to M.  If the current
     * state is not ok(), or @p n is out of range, the method returns zero.
     */
    double coefficient (int n);
private:

    BSpline<double>* spline_;

  };

}

#endif
