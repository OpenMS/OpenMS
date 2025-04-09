// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <vector>

#include <OpenMS/config.h>
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

