// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
    @brief _MRMFeaturePicker_ defines the structures containing parameters to be used in
    [MRMTransitionGroupPicker](@ref MRMTransitionGroupPicker) for components and components groups.

    This data can be loaded from a file with [MRMFeaturePickerFile](@ref MRMFeaturePickerFile).

    Examples of parameters are:
    "TransitionGroupPicker:compute_peak_quality"
    "TransitionGroupPicker:stop_after_feature"
    "TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise"
    "TransitionGroupPicker:PeakPickerChromatogram:sn_bin_count"
  */
  class OPENMS_DLLAPI MRMFeaturePicker
  {
public:
    /// Constructor
    MRMFeaturePicker() = default;

    /// Destructor
    ~MRMFeaturePicker() = default;

    /// Structure to contain information about a single component with its parameters
    struct ComponentParams
    {
      String component_name; ///< The component_name can't be an empty string
      String component_group_name; ///< The component_group_name can't be an empty string
      Param params; ///< The parameters pertaining a single component
    };

    /// Structure to contain information about a component group with its parameters
    struct ComponentGroupParams
    {
      String component_group_name; ///< The component_group_name can't be an empty string
      Param params; ///< The parameters pertaining a component group
    };
  };
}

