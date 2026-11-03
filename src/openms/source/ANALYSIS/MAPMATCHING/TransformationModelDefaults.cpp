// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Clemens Groepl, Hendrik Weisser, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelDefaults.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelBSpline.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLowess.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelInterpolated.h>

namespace OpenMS
{
  Param TransformationModelDefaults::getDefaults(const std::string& default_model)
  {
    Param params;
    params.setValue("type", default_model, "Type of model");
    // TODO: avoid referring to each TransformationModel subclass explicitly
    std::vector<std::string> model_types = {"linear","b_spline","lowess","interpolated"};
    if (!ListUtils::contains(model_types, default_model))
    {
      model_types.insert(model_types.begin(), default_model);
    }
    params.setValidStrings("type", model_types);

    Param model_params;
    TransformationModelLinear::getDefaultParameters(model_params);
    params.insert("linear:", model_params);
    params.setSectionDescription("linear", "Parameters for 'linear' model");

    TransformationModelBSpline::getDefaultParameters(model_params);
    params.insert("b_spline:", model_params);
    params.setSectionDescription("b_spline", "Parameters for 'b_spline' model");

    TransformationModelLowess::getDefaultParameters(model_params);
    params.insert("lowess:", model_params);
    params.setSectionDescription("lowess", "Parameters for 'lowess' model");

    TransformationModelInterpolated::getDefaultParameters(model_params);
    params.insert("interpolated:", model_params);
    params.setSectionDescription("interpolated",
                                "Parameters for 'interpolated' model");
    return params;
  }

} // namespace OpenMS
