// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // is this needed?
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>

namespace OpenMS
{

  /**
    @brief Lowess (non-linear) model for transformations

    @ingroup MapAlignment
  */
  class OPENMS_DLLAPI TransformationModelLowess :
    public TransformationModel
  {
public:
    /**
      @brief Constructor

      @exception Exception::IllegalArgument is thrown if too few data points are provided.
    */
    TransformationModelLowess(const DataPoints& data, const Param& params);

    /// Destructor
    ~TransformationModelLowess() override;

    /// Evaluates the model at the given value
    double evaluate(double value) const override
    {
      return model_->evaluate(value);
    }

    using TransformationModel::getParameters;

    /// Gets the default parameters
    static void getDefaultParameters(Param& params);

protected:
    /// Pointer to the underlying interpolation
    TransformationModelInterpolated* model_;

  };
} // namespace

