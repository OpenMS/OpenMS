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

#include <OpenMS/config.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{

  /**
    @brief AbsoluteQuantitationStandards is a class to handle the relationship between
    runs, components, and their actual concentrations.

    A mapping between a run, the components in the run, and the actual concentration
    of the components in the run are required to build a calibration curve that is
    required for absolute quantitation.
  */
  class OPENMS_DLLAPI AbsoluteQuantitationStandards
  {

public:
    /// Constructor
    AbsoluteQuantitationStandards() = default;

    /// Destructor
    ~AbsoluteQuantitationStandards() = default;

    /// Structure to map runs to components to known concentrations.
    struct runConcentration
    {
      String sample_name;
      String component_name;
      String IS_component_name;
      double actual_concentration;
      double IS_actual_concentration;
      String concentration_units;
      double dilution_factor;
    };

    /// Structure to hold a single component and its corresponding known concentration.
    struct featureConcentration
    {
      Feature feature;
      Feature IS_feature;
      double actual_concentration;
      double IS_actual_concentration;
      String concentration_units;
      double dilution_factor;
    };

    /**
      @brief Method to map runs to components to known concentrations.

      @warning The method checks for the FeatureMaps' sample names with FeatureMap::getPrimaryMSRunPath()

      @param[in] run_concentrations A list of runConcentration structs (e.g., from file upload).
      @param[in] feature_maps The method maps to these features.
      @param[out] components_to_concentrations A map that links run data to feature data.
    */
    void mapComponentsToConcentrations(
      const std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations,
      const std::vector<FeatureMap>& feature_maps,
      std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>>& components_to_concentrations
    ) const;

    /**
      @brief Get the feature concentrations from a single component.

      This method internally calls `mapComponentsToConcentrations()`, but takes in consideration only those
      elements of `run_concentrations` that have the passed `component_name`.

      @warning The method checks for the FeatureMaps' sample names with FeatureMap::getPrimaryMSRunPath()

      @param[in] run_concentrations A list of runConcentration structs (e.g., from file upload).
      @param[in] feature_maps The method maps to these features.
      @param[in] component_name Only runConcentration with this name will be considered.
      @param[out] feature_concentrations The list of feature concentrations found.
    */
    void getComponentFeatureConcentrations(
      const std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations,
      const std::vector<FeatureMap>& feature_maps,
      const String& component_name,
      std::vector<AbsoluteQuantitationStandards::featureConcentration>& feature_concentrations
    ) const;

private:
    /**
      @brief Finds a feature for a given component name.

      @param[in] feature_map The container of features.
      @param[in] component_name The feature must have this name as its "native_id".
      @param[out] feature_found If found, the feature is saved in this parameter.
    */
    bool findComponentFeature_(
      const FeatureMap& feature_map,
      const String& component_name,
      Feature& feature_found
    ) const;
  };
}

