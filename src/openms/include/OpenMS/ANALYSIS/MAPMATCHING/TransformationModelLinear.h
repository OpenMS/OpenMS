// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{

  /**
    @brief Linear model for transformations

    The model can be inferred from data or specified using explicit parameters. 
    If data is given, a least squares fit is used to find the model parameters (slope and intercept). 
    Depending on parameter @p symmetric_regression, a normal regression (@e y on @e x) or
    symmetric regression (@f$ y - x @f$ on @f$ y + x @f$) is performed.

    Without data, the model can be specified by giving the parameters @p slope, @p intercept, 
    @p x_weight, @p y_weight explicitly.

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModelLinear :
    public TransformationModel
  {
public:
    /**
      @brief Constructor

      @exception IllegalArgument is thrown if neither data points nor explicit parameters (slope/intercept) are given.
    */
    TransformationModelLinear(const DataPoints& data, const Param& params);

    /// Destructor
    ~TransformationModelLinear() override = default;

    /// Evaluates the model at the given value
    double evaluate(double value) const override;

    using TransformationModel::getParameters;

    /// Gets the "real" parameters
    void getParameters(double& slope, double& intercept, String& x_weight, String& y_weight, double& x_datum_min, double& x_datum_max, double& y_datum_min, double& y_datum_max) const;

    /// Gets the default parameters
    static void getDefaultParameters(Param& params);

    /**
     @brief Computes the inverse

     @exception DivisionByZero is thrown if the slope is zero.
    */
    void invert();

protected:
    /// Parameters of the linear model
    double slope_, intercept_;
    /// Was the model estimated from data?
    bool data_given_;
    /// Use symmetric regression?
    bool symmetric_;
  };
} // namespace

