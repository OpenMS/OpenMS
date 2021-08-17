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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/TARGETED/MetaboTargetedTargetDecoy.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <iostream>
#include <regex>

namespace OpenMS
{
  std::map<String, std::vector<OpenMS::ReactionMonitoringTransition> > MetaboTargetedTargetDecoy::constructTransitionsMap_(const TargetedExperiment& t_exp)
  {
    // mapping of the transitions to a specific compound reference
    std::map<String, std::vector<OpenMS::ReactionMonitoringTransition> > TransitionsMap;
    for (const auto& tr_it : t_exp.getTransitions())
    {
      auto pair_it_success = TransitionsMap.emplace(tr_it.getCompoundRef(), std::vector<OpenMS::ReactionMonitoringTransition>());
      pair_it_success.first->second.push_back(tr_it);
    }
    return TransitionsMap;
  }

  std::vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping> MetaboTargetedTargetDecoy::constructTargetDecoyMassMapping(const TargetedExperiment& t_exp)
  {
    std::vector<String> identifier;
    for (const auto &it : t_exp.getCompounds())
    {
      // only need to extract identifier from the targets, since targets and decoys have the same
      if (it.getMetaValue("decoy") == DataValue(0))
      {
        identifier.emplace_back(it.getMetaValue("m_ids_id"));
      }
    }

    std::vector<ReactionMonitoringTransition> rmts = t_exp.getTransitions();
    std::vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping> mappings;
    for (const auto& it : identifier)
    {
      MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping mapping;
      mapping.identifier = it;
      auto it_target = rmts.begin();
      while ((it_target = find_if(it_target,
                                  rmts.end(),
                                  [&it](ReactionMonitoringTransition &rts)
                                  {
                                    return rts.getMetaValue("m_ids_id") == it &&
                                           rts.getDecoyTransitionType() == ReactionMonitoringTransition::TARGET;
                                  })) != rmts.end())
      {
        mapping.target_product_masses.emplace_back(it_target->getProductMZ());
        mapping.target_compound_ref = it_target->getCompoundRef();
        ++it_target;
      }
      auto it_decoy = rmts.begin();
      while ((it_decoy = find_if(it_decoy,
                                 rmts.end(),
                                 [&it](ReactionMonitoringTransition &rts)
                                 {
                                   return rts.getMetaValue("m_ids_id") == it &&
                                          rts.getDecoyTransitionType() == ReactionMonitoringTransition::DECOY;
                                 })) != rmts.end())
      {
        mapping.decoy_product_masses.emplace_back(it_decoy->getProductMZ());
        mapping.decoy_compound_ref = it_decoy->getCompoundRef();
        ++it_decoy;
      }
      mappings.emplace_back(mapping);
    }
    return mappings;
  }

  void MetaboTargetedTargetDecoy::resolveOverlappingTargetDecoyMassesByIndividualMassShift(TargetedExperiment& t_exp, std::vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping>& mappings, const double& mass_to_add)
  {
    // compare targets and decoy masses
    for (auto &it : mappings)
    {
      std::vector<double> intersection;
      std::vector<double> replace_intersection;
      std::vector<double> target_product_masses = it.target_product_masses;
      std::vector<double> decoy_product_masses = it.decoy_product_masses;
      std::sort(target_product_masses.begin(), target_product_masses.end());
      std::sort(decoy_product_masses.begin(), decoy_product_masses.end());

      std::set_intersection(target_product_masses.begin(),
                            target_product_masses.end(),
                            decoy_product_masses.begin(),
                            decoy_product_masses.end(),
                            std::back_inserter(intersection));

      if (!intersection.empty())
      {
        std::transform(intersection.begin(),
                       intersection.end(),
                       std::back_inserter(replace_intersection),
                       [mass_to_add](double d) -> double { return d + mass_to_add; });

        for (Size i = 0; i < replace_intersection.size(); ++i)
        {
          std::replace(it.decoy_product_masses.begin(),
                       it.decoy_product_masses.end(),
                       intersection[i],
                       replace_intersection[i]);
        }
      }
    }

    std::map<String, std::vector<OpenMS::ReactionMonitoringTransition> > TransitionsMap = MetaboTargetedTargetDecoy::constructTransitionsMap_(t_exp);

    // resolve mappings and add to current TargetedExperiment
    std::vector<OpenMS::ReactionMonitoringTransition> transitions;

    for (const auto& it: mappings)
    {
      if (!it.target_compound_ref.empty())
      {
        std::vector<OpenMS::ReactionMonitoringTransition> target_transitions = TransitionsMap[it.target_compound_ref];
        transitions.insert(transitions.end(), target_transitions.begin(), target_transitions.end());
      }
      if (!it.decoy_compound_ref.empty())
      {
        std::vector<OpenMS::ReactionMonitoringTransition> current_decoy_transitions = TransitionsMap[it.decoy_compound_ref];
        if (it.decoy_product_masses.size() == current_decoy_transitions.size())
        {
          for (size_t i = 0; i < it.decoy_product_masses.size(); ++i)
          {
            ReactionMonitoringTransition tr = current_decoy_transitions[i]; // old
            tr.setProductMZ(it.decoy_product_masses[i]); // new
            transitions.push_back(tr);
          }
        }
      }
    }
    t_exp.setTransitions(transitions);
  }

  void MetaboTargetedTargetDecoy::generateMissingDecoysByMassShift(TargetedExperiment& t_exp, std::vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping>& mappings, const double& mass_to_add)
  {
    // Add a decoy based on the target masses + mass_to_add (e.g. CH2) if fragmentation tree re-rooting was not possible
    for (auto &it : mappings)
    {
      if (it.decoy_product_masses.size() != it.target_product_masses.size())
      {
        // decoy was not generated by passatutto
        if (it.decoy_compound_ref.empty())
        {
          // add a potential decoy with the new decoy masses to the mapping
          it.decoy_compound_ref = std::regex_replace(it.target_compound_ref, std::regex(R"(_\[)"), "_decoy_[");
          std::transform(it.target_product_masses.begin(),
                         it.target_product_masses.end(),
                         std::back_inserter(it.decoy_product_masses),
                         [mass_to_add](double d) -> double { return d + mass_to_add; });
        }
      }
    }

    std::map<String, std::vector<OpenMS::ReactionMonitoringTransition> > TransitionsMap = MetaboTargetedTargetDecoy::constructTransitionsMap_(t_exp);

    std::vector<TargetedExperiment::Compound> compounds;
    std::vector<ReactionMonitoringTransition> transitions;
    // look if compounds exists as target and decoy
    // add it to the current TargetedExperiment
    for (const auto &it : mappings)
    {
      const auto it_target = std::find_if(t_exp.getCompounds().begin(),
                                          t_exp.getCompounds().end(),
                                          [&it](const TargetedExperiment::Compound &comp)
                                          {
                                            return comp.id == it.target_compound_ref;
                                          });
      const auto it_decoy = std::find_if(t_exp.getCompounds().begin(),
                                         t_exp.getCompounds().end(),
                                         [&it](const TargetedExperiment::Compound &comp)
                                         {
                                           return comp.id == it.decoy_compound_ref;
                                         });

      // if targets and decoy exists add them to the new datastructure
      if (it_target != t_exp.getCompounds().end())
      {
        compounds.emplace_back(*it_target);
        if (TransitionsMap.find(it_target->id) != TransitionsMap.end())
        {
          transitions.insert(transitions.end(),
                               TransitionsMap[it_target->id].begin(),
                               TransitionsMap[it_target->id].end());
        }
      }
      if (it_decoy != t_exp.getCompounds().end())
      {
        compounds.emplace_back(*it_decoy);
        if (TransitionsMap.find(it_decoy->id) != TransitionsMap.end())
        {
          transitions.insert(transitions.end(),
                               TransitionsMap[it_decoy->id].begin(),
                               TransitionsMap[it_decoy->id].end());
        }
      }
      else // decoy does not exists in TargetedExperimentCompound
      {
        // use the corresponding target compound to generate a new decoy
        // and add the decoy transitions.
        TargetedExperiment::Compound potential_decoy_compound = *it_target;
        std::vector<ReactionMonitoringTransition> potential_decoy_transitions;

        String current_compound_name = potential_decoy_compound.getMetaValue("CompoundName");
        potential_decoy_compound.setMetaValue("CompoundName", String(current_compound_name + "_decoy"));
        potential_decoy_compound.id = it.decoy_compound_ref;
        potential_decoy_compound.setMetaValue("decoy", DataValue(1));

        if (TransitionsMap.find(it_target->id) != TransitionsMap.end())
        {
          potential_decoy_transitions = TransitionsMap[it_target->id];
          for (size_t i = 0; i < potential_decoy_transitions.size(); ++i)
          {
            potential_decoy_transitions[i]
                .setNativeID(std::regex_replace(potential_decoy_transitions[i].getNativeID(),
                                                std::regex(R"(_\[)"),
                                                "_decoy_["));
            potential_decoy_transitions[i]
                .setDecoyTransitionType(ReactionMonitoringTransition::DecoyTransitionType::DECOY);
            potential_decoy_transitions[i].setMetaValue("annotation", "NA");
            potential_decoy_transitions[i].setProductMZ(it.decoy_product_masses[i]);
            potential_decoy_transitions[i].setCompoundRef(it.decoy_compound_ref);
          }
        }
        else
        {
          OPENMS_LOG_WARN << "Add_shift method failed: " << current_compound_name << "_decoy could not be generated." << std::endl;
        }
        compounds.emplace_back(potential_decoy_compound);
        transitions.insert(transitions.end(), potential_decoy_transitions.begin(), potential_decoy_transitions.end());
      }
    }
    t_exp.setCompounds(compounds);
    t_exp.setTransitions(transitions);
  }

} // namespace OpenMS
