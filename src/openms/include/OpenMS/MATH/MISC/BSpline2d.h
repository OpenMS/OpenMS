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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_MISC_BSPLINE2D_H
#define OPENMS_MATH_MISC_BSPLINE2D_H

#include <vector>

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/StandardTypes.h>

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
    
    // Note: Don't change boundary coniditon constants ase these are passed through to the eol-bspline implementation.
    enum BoundaryCondition
    {
      /// Set the endpoints of the spline to zero.
      BC_ZERO_ENDPOINTS = 0,
      /// Set the first derivative of the spline to zero at the endpoints.
      BC_ZERO_FIRST = 1,
      /// Set the second derivative to zero.
      BC_ZERO_SECOND = 2
    };

    /**
     * Create a single spline with the parameters required to set up
     * the domain and subsequently smooth the given set of y values.
     * The y values must correspond to each of the values in the x array.
     * If either the domain setup fails or the spline cannot be solved,
     * the state will be set to not ok.
     *
     * @see ok().
     *
     * @param x		The array of x values in the domain.
     * @param y		The array of y values corresponding to each of the
     *			nX() x values in the domain.
     * @param wavelength	The cutoff wavelength, in the same units as the
     *				@p x values.  A wavelength of zero disables
     *				the derivative constraint.
     * @param bc_type	The enumerated boundary condition type.  If
     *			omitted it defaults to BC_ZERO_SECOND.
     * @param num_nodes The number of nodes to use for the cubic b-spline.
     *			If less than 2 a "reasonable" number will be
     *			calculated automatically, taking into account
     *			the given cutoff wavelength.                                                                                                               
     **/
    BSpline2d(const std::vector<double>& x, const std::vector<double>& y, double wave_length = 0, BoundaryCondition boundary_condition = BC_ZERO_ENDPOINTS, Size num_nodes = 0);

    /**
     * Destructor
     */
    virtual ~BSpline2d();

    /**
     * Solve the spline curve for a new set of y values.  Returns false
     * if the solution fails.
     *
     * @param y The array of y values corresponding to each of the nX()
     *		x values in the domain.
     */
    bool solve (const std::vector<double>& y);

    /**
     * Return the evaluation of the smoothed curve
     * at a particular @p x value.  If current state is not ok(), returns 0.
     */
    double eval (double x);

    /**
     * Return the first derivative of the spline curve at the given @p x.
     * Returns zero if the current state is not ok().
     */
    double derivative (double x);

    /**
     * Return the @p n-th basis coefficient, from 0 to M.  If the current
     * state is not ok(), or @p n is out of range, the method returns zero.
     */
    double coefficient (int n);
private:
    BSpline2d();

    BSpline<double>* spline_;

  };

}

#endif
