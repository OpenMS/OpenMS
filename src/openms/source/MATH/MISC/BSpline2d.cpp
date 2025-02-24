// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/MISC/BSpline2d.h>

#include <BSpline/BSplineBase.cpp>
#include <BSpline/BSpline.cpp>

namespace OpenMS
{

  BSpline2d::BSpline2d(const std::vector<double>& x, const std::vector<double>& y, 
                       double wavelength, BoundaryCondition boundary_condition, 
                       Size num_nodes)
  {
    OPENMS_PRECONDITION(x.size() == y.size(), "x and y vectors passed to BSpline2d constructor must have the same size.")
    spline_ = new eol_bspline::BSpline<double>(&x[0], static_cast<int>(x.size()), &y[0], wavelength, boundary_condition, num_nodes);
  }

  BSpline2d::~BSpline2d()
  {
    delete spline_;
  }

  bool BSpline2d::solve(const std::vector<double>& y)
  {
    OPENMS_PRECONDITION(static_cast<Size>(spline_->nX()) == y.size(), "y vector passed to 'BSpline2d::solve' must match size of x.")
    // pass vector as array
    return spline_->solve(&y[0]);
  }

  double BSpline2d::eval(const double x) const
  {
    OPENMS_PRECONDITION(ok(), "Spline was not initialized properly.")
    return spline_->evaluate(x);
  }

  double BSpline2d::derivative(const double x) const
  {
    OPENMS_PRECONDITION(ok(), "Spline was not initialized properly.")
    return spline_->slope(x);
  }

  bool BSpline2d::ok() const
  {
    return spline_->ok();
  }

  void BSpline2d::debug(bool enable)
  {
    eol_bspline::BSplineBase<double>::Debug(int(enable));
  }

}
