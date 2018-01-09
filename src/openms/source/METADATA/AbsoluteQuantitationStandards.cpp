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
  void AbsoluteQuantitationStandards::mapComponentsToConcentrations(
    const std::vector<AbsoluteQuantitationStandards::runConcentration>& run_concentrations,
    const std::vector<FeatureMap>& feature_maps,
    std::map<String, AbsoluteQuantitationStandards::featureConcentration>& components_to_concentrations
  ) const
  {
    components_to_concentrations.clear();
    for (const FeatureMap& fmap : feature_maps)
    {
      if (!fmap.metaValueExists("sample_name"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The FeatureMap misses the 'sample_name' MetaValue information.");
      }
      const String sample_name = fmap.getMetaValue("sample_name");
      for (const Feature& feature : fmap)
      {
        const std::vector<Feature>& subordinates = feature.getSubordinates();
        if (subordinates.size() != 2)
        {
          continue;
        }
        const Feature& f1 = subordinates[0];
        const Feature& f2 = subordinates[1];
        std::vector<AbsoluteQuantitationStandards::runConcentration>::const_iterator it;
        it = std::find_if(
          run_concentrations.begin(), run_concentrations.end(), [&sample_name, &f1, &f2] (AbsoluteQuantitationStandards::runConcentration run)
          {
            if (!(f1.metaValueExists("native_id") && f2.metaValueExists("native_id")))
            {
              //throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The feature misses the 'native_id' MetaValue information.");
              LOG_DEBUG << "The feature misses the 'native_id' MetaValue information.\n";
              return false;
            }
            return sample_name == run.sample_name &&
              (
                (f1.getMetaValue("native_id") == run.component_name && f2.getMetaValue("native_id") == run.IS_component_name) ||
                (f2.getMetaValue("native_id") == run.component_name && f1.getMetaValue("native_id") == run.IS_component_name)
              );
          }
        );
        if (it != run_concentrations.end())
        {
          AbsoluteQuantitationStandards::featureConcentration fc;
          if (f1.getMetaValue("native_id") == it->component_name)
          {
            fc.feature = subordinates[0];
            fc.IS_feature = subordinates[1];
          }
          else
          {
            fc.feature = subordinates[1];
            fc.IS_feature = subordinates[0];
          }
          fc.actual_concentration = it->actual_concentration;
          fc.IS_actual_concentration = it->IS_actual_concentration;
          fc.concentration_units = it->concentration_units;
          fc.dilution_factor = it->dilution_factor;
          std::pair<std::map<String, AbsoluteQuantitationStandards::featureConcentration>::const_iterator, bool> p;
          p = components_to_concentrations.insert({it->component_name, fc});
          if (p.second == false)
          {
            throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Key '" + it->component_name + "' was already present.");
          }
        }
      }
    }
  }
} // namespace

