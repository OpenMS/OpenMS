// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLFeatureAugmentation.h>

#include <stdexcept>
#include <vector>

namespace OpenMS
{
  class OPENMS_DLLAPI NuXLFeatureAugmentation
  {
    public:
      static void augment(std::vector<PeptideIdentification>& pep_ids, 
      std::vector<std::string> positive_weights,
      std::vector<std::string> negative_weights)
      {
        // only for XLs? because they are fewer?
        if (pep_ids.empty()) return;

        if (pep_ids[0].getHits().empty()) return;

        // feature names may not intersect between postive/negative constrained
        std::sort(positive_weights.begin(), positive_weights.end());
        std::sort(negative_weights.begin(), negative_weights.end());
        std::vector<std::string> v_intersection;
 
        std::set_intersection(positive_weights.begin(), positive_weights.end(),
                          negative_weights.begin(), negative_weights.end(),
                          std::back_inserter(v_intersection));

        if (!v_intersection.empty()) throw std::runtime_error("Positive and negative weights may not overlap.");

        // use first PSM as template
        auto p_template = pep_ids[0].getHits()[0];
        p_template.setScore(0);
        std::vector<String> keys;
        p_template.getKeys(keys);
        p_template.setMetaValue("NuXL:augmented", "true");

        // clear scores of template
        for (const auto& k : keys)
        { 
          if (p_template.getMetaValue(k).valueType() == DataValue::INT_VALUE) p_template.setMetaValue(k, 0);
          if (p_template.getMetaValue(k).valueType() == DataValue::DOUBLE_VALUE) p_template.setMetaValue(k, 0.0);
        }

        // determine minium and maximum for each of the positive constrained features
        std::map<String, double> minima;
        std::map<String, double> maxima;

        for (const auto& k : positive_weights) 
        {
          minima[k] = 1e32;
          maxima[k] = -1e32;
        }

        for (const auto& pid : pep_ids)
        {
          for (const auto& ph : pid.getHits())
          {
            for (const auto& k : positive_weights) 
            {
              auto dv = ph.getMetaValue(k);
              if (dv.valueType() == DataValue::INT_VALUE) 
              { 
                if (minima[k] > (int)dv) minima[k] = (int)dv; 
                if (maxima[k] < (int)dv) maxima[k] = (int)dv; 
              };
              if (dv.valueType() == DataValue::DOUBLE_VALUE) 
              { 
                if (minima[k] > (double)dv) minima[k] = (double)dv; 
                if (maxima[k] < (double)dv) maxima[k] = (double)dv; 
              };
            }    
          }
        }

        size_t c = 0;
        // for each positive_weight feature, create one example with that feature set to max value
        for (const auto& s : positive_weights) 
        {
          auto p = p_template;
          p.setMetaValue(s, maxima[s] + 1000.0 * (maxima[s] - minima[s])); // set feature value to >> maximum of observed ones
          std::vector<PeptideHit> phs;
          phs.push_back(p);
          PeptideIdentification pid = pep_ids[0];
          pid.setRT(1e6 + c); // RT of augmented example
          pid.setHits(phs);
          pep_ids.push_back(pid);
          ++c;
        } 

        // for each negative_weight feature, create one example with that feature set to min value
        for (const auto& s : negative_weights) 
        {
          auto p = p_template;
          p.setMetaValue(s, minima[s] - 1000.0 * (maxima[s] - minima[s])); // set feature value to << min of observed ones
          std::vector<PeptideHit> phs;
          phs.push_back(p);
          PeptideIdentification pid = pep_ids[0];
          pid.setRT(1e6 + c); // RT of augmented example
          pid.setHits(phs);
          pep_ids.push_back(pid);
          ++c;
        } 
      }
 
      static void removeAugmented(std::vector<PeptideIdentification>& pep_ids)
      {
        // remove augmented features again
        for (auto& pid : pep_ids)
        {
          auto& phs = pid.getHits();
          phs.erase(remove_if(phs.begin(), phs.end(), [](const PeptideHit& ph){ return ph.metaValueExists("NuXL:augmented"); }), phs.end());
        }
      }
  };
}

