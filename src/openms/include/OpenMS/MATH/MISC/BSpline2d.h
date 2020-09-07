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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <vector>

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/StandardDeclarations.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

// forward declaration of impl class BSpline
namespace eol_bspline
{
template <class T>
class BSpline;
}

namespace OpenMS
{
  /**
   * @brief b spline interpolation
   */
  class OPENMS_DLLAPI BSpline2d
  {
public:

    // Note: Don't change boundary condition constants as these are passed through to the eol-bspline implementation.
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
     * The y values must correspond to each of the x values.
     * If either the domain setup fails or the spline cannot be solved,
     * the state will be set to "not OK".
     *
     * @see ok().
     *
     * @param x		The array of x values in the domain.
     * @param y		The array of y values corresponding to each of the
     *			x values in the domain.
     * @param wavelength	The cutoff wavelength, in the same units as the
     *				@p x values.  A wavelength of zero disables
     *				the derivative constraint.
     * @param boundary_condition	The boundary condition type. If
     *			omitted it defaults to BC_ZERO_SECOND.
     * @param num_nodes The number of nodes to use for the cubic b-spline.
     *			If less than 2, a "reasonable" number will be
     *			calculated automatically, taking into account
     *			the given cutoff wavelength.
     * @pre x and y must be of the same dimensions.
     **/
    BSpline2d(const std::vector<double>& x, const std::vector<double>& y,
              double wavelength = 0, BoundaryCondition boundary_condition = BC_ZERO_SECOND, 
              Size num_nodes = 0);

    /**
     * Destructor
     */
    virtual ~BSpline2d();

    /**
     * Solve the spline curve for a new set of y values.  Returns false
     * if the solution fails.
     *
     * @param y The array of y values corresponding to each of the
     *		x values in the domain.
     */
    bool solve(const std::vector<double>& y);

    /**
     * Return the evaluation of the smoothed curve
     * at a particular @p x value. If current state is not ok(), returns zero.
     */
    double eval(const double x) const;

    /**
     * Return the first derivative of the spline curve at the given @p x.
     * Returns zero if the current state is not ok().
     */
    double derivative(const double x) const;

    /**
     * Return the first derivative of the spline curve at the given @p x.
     * Returns zero if the current state is not ok().
     */
    double derivatives(double x, unsigned order = 1) const;

    /**
     * Return whether the spline fit was successful.
     */
    bool ok() const;

    /**
     * Enable or disable debug messages from the B-spline library.
     */
    static void debug(bool enable);

private:

    // Pointer to actual implementation. Note: This class follows the PIMPL idiom hiding the actual 
    // B-spline implementation behind this pointer to avoid any dependency of the interface to the 
    // implementation. Thus, the eol splines are only required during compilation of OpenMS and 
    // not when linking against OpenMS.
    eol_bspline::BSpline<double>* spline_;
  };

}

