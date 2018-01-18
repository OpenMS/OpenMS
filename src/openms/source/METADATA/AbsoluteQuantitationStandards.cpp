// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
    std::map<String, AbsoluteQuantitationStandards::featureConcentration>& components_to_concentrations
  ) const
  {
    components_to_concentrations.clear();
    for (const AbsoluteQuantitationStandards::runConcentration& run : run_concentrations)
    {
      if (run.sample_name == "" || run.component_name == "")
      {
        continue;
      }
      for (const FeatureMap& fmap : feature_maps) // not all elements are necessarily processed (break; is present inside the loop)
      {
        StringList sample_name;
        fmap.getPrimaryMSRunPath(sample_name);
        if (!sample_name.size() || sample_name[0] != run.sample_name) // if the FeatureMap doesn't have a sample_name, or if it is not the one we're looking for: skip.
        {
          continue;
        }
        AbsoluteQuantitationStandards::featureConcentration fc;
        if (!findComponentFeature_(fmap, run.component_name, fc.feature)) // if there was no match: skip.
        {
          continue;
        }
        if (run.IS_component_name != "")
        {
          findComponentFeature_(fmap, run.IS_component_name, fc.IS_feature);
        }
        // fill the rest of the information from the current runConcentration
        fc.actual_concentration = run.actual_concentration;
        fc.IS_actual_concentration = run.IS_actual_concentration;
        fc.concentration_units = run.concentration_units;
        fc.dilution_factor = run.dilution_factor;
        // add to the map
        std::pair<std::map<String, AbsoluteQuantitationStandards::featureConcentration>::const_iterator, bool> p;
        p = components_to_concentrations.insert({run.component_name, fc});
        if (p.second == false) // check that the key was not already present
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Key '" + run.component_name + "' was already present.");
        }
        break;
      }
    }
  }
} // namespace

