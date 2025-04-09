// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>

namespace OpenMS
{
  bool AbsoluteQuantitationStandards::findComponentFeature_(
    const FeatureMap& feature_map,
    const String& component_name,
    Feature& feature_found
  ) const
  {
    for (const Feature& feature : feature_map)
    {
      for (const Feature& subordinate : feature.getSubordinates())
      {
        if (subordinate.metaValueExists("native_id") && subordinate.getMetaValue("native_id") == component_name)
        {
          feature_found = subordinate;
          return true;
        }
      }
    }
    return false;
  }

  void AbsoluteQuantitationStandards::mapComponentsToConcentrations(
    const std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations,
    const std::vector<FeatureMap>& feature_maps,
    std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>>& components_to_concentrations
  ) const
  {
    components_to_concentrations.clear();
    for (const AbsoluteQuantitationStandards::runConcentration& run : run_concentrations)
    {
      if (run.sample_name.empty() || run.component_name.empty())
      {
        continue;
      }
      for (const FeatureMap& fmap : feature_maps) // not all elements are necessarily processed (break; is present inside the loop)
      {
        StringList filename;
        fmap.getPrimaryMSRunPath(filename);
        if (!filename.empty()) // if the FeatureMap doesn't have a sample_name, or if it is not the one we're looking for: skip.
        {
          if (filename[0].hasSuffix(".mzML"))
          {
            filename[0].resize(filename[0].size() - 5);
          }
          else if (filename[0].hasSuffix(".txt"))
          {
            filename[0].resize(filename[0].size() - 4);
          }
          if (filename[0] != run.sample_name)
          {
            continue;
          }
        }
        AbsoluteQuantitationStandards::featureConcentration fc;
        if (!findComponentFeature_(fmap, run.component_name, fc.feature)) // if there was no match: skip.
        {
          continue;
        }
        if (!run.IS_component_name.empty())
        {
          findComponentFeature_(fmap, run.IS_component_name, fc.IS_feature);
        }
        // fill the rest of the information from the current runConcentration
        fc.actual_concentration = run.actual_concentration;
        fc.IS_actual_concentration = run.IS_actual_concentration;
        fc.concentration_units = run.concentration_units;
        fc.dilution_factor = run.dilution_factor;
        // add to the map
        std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>>::iterator p;
        p = components_to_concentrations.find(run.component_name);
        if (p == components_to_concentrations.end()) // if the key doesn't exist, insert it and create a new vector with fc as its only element
        {
          components_to_concentrations.insert({run.component_name, {fc}});
        }
        else // otherwise, add the element to the existing vector
        {
          (p->second).push_back(fc);
        }
        break; // because there won't be another FeatureMap with the same sample_name
      }
    }
  }

  void AbsoluteQuantitationStandards::getComponentFeatureConcentrations(
    const std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations,
    const std::vector<FeatureMap>& feature_maps,
    const String& component_name,
    std::vector<AbsoluteQuantitationStandards::featureConcentration>& feature_concentrations
  ) const
  {
    std::vector<AbsoluteQuantitationStandards::runConcentration> filtered_rc;
    for (const AbsoluteQuantitationStandards::runConcentration& run : run_concentrations)
    {
      if (run.component_name == component_name)
      {
        filtered_rc.push_back(run);
      }
    }
    std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>> components_to_concentrations;
    mapComponentsToConcentrations(filtered_rc, feature_maps, components_to_concentrations);
    if (components_to_concentrations.count(component_name))
    {
      feature_concentrations = components_to_concentrations.at(component_name);
    }
  }
} // namespace

