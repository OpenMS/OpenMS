// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // is this needed?

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/MATH/MISC/BSpline2d.h>

namespace OpenMS
{

  /**
    @brief B-spline (non-linear) model for transformations

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModelBSpline :
    public TransformationModel
  {
public:
    /**
      @brief Constructor

      @exception Exception::IllegalArgument is thrown if a parameter is invalid.
      @exception Exception::UnableToFit is thrown if the B-spline fit fails.
    */
    TransformationModelBSpline(const DataPoints& data, const Param& params);

    /// Destructor
    ~TransformationModelBSpline() override;

    /// Evaluates the model at the given value
    double evaluate(double value) const override;

    using TransformationModel::getParameters;

    /// Gets the default parameters
    static void getDefaultParameters(Param& params);

protected:
    /// Pointer to the actual B-spline
    BSpline2d* spline_;

    /// Min./max. x value (endpoints of the data range)
    double xmin_, xmax_;

    /// Method to use for extrapolation (beyond 'xmin_'/'xmax_')
    enum { EX_LINEAR, EX_BSPLINE, EX_CONSTANT, EX_GLOBAL_LINEAR } extrapolate_;

    /// Parameters for constant or linear extrapolation 
    double offset_min_, offset_max_, slope_min_, slope_max_;
  };
} // namespace

