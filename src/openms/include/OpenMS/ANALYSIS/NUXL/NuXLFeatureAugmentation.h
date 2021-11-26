// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
          p.setMetaValue(s, maxima[s]); // set feature value to maximum of observed ones
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
          p.setMetaValue(s, minima[s]); // set feature value to min of observed ones
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

