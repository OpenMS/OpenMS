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
    "TransitionGroupPicker:PeakPickerMRM:signal_to_noise"
    "TransitionGroupPicker:PeakPickerMRM:sn_bin_count"
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

