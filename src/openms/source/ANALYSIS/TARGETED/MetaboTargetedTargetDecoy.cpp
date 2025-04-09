// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

  void MetaboTargetedTargetDecoy::resolveOverlappingTargetDecoyMassesByDecoyMassShift(TargetedExperiment& t_exp, std::vector<MetaboTargetedTargetDecoy::MetaboTargetDecoyMassMapping>& mappings, const double& mass_to_add, const double& mz_tol, const String& mz_tol_unit)
  {
    // Define a map to hold compound references and their corresponding sets of decoy m/z values.
    std::map<String, std::set<double>> match_compound_refs_decoy_mz;

    // Iterate over each mapping in the provided mappings list.
    for (const auto& map : mappings) {
        // Create a set to store m/z values that match the criterion.
        std::set<double> matched;

        // Iterate over each decoy m/z value in the current mapping.
        for (double decoy_mz : map.decoy_product_masses) {
            // Compare each decoy m/z value with each target m/z value.
            for (double target_mz : map.target_product_masses) {
                // Calculate the difference between decoy and target m/z values.
                // The calculation differs based on whether the tolerance is in ppm or Da.
                double difference = (mz_tol_unit == "ppm") ?
                    std::abs(decoy_mz - target_mz) / target_mz * 1e6 :
                    std::abs(decoy_mz - target_mz);

                // If the difference is small than mz_tol, the masses are too similar.
                if (difference <= mz_tol) {
                    // Add the decoy m/z to the matched set.
                    matched.insert(decoy_mz);
                    break; // Move to the next decoy m/z value after finding a match.
                }
            }
        }

        // Associate the set of matched decoy m/z values with the compound reference in the map.
        match_compound_refs_decoy_mz[map.decoy_compound_ref] = std::move(matched);
    }

    // Prepare a new vector to store updated ReactionMonitoringTransition objects.
    std::vector<ReactionMonitoringTransition> v_rmt_new;
    v_rmt_new.reserve(t_exp.getTransitions().size()); // Reserve space to optimize memory allocation.

    // Iterate over each transition in the experiment.
    for (const auto& tr : t_exp.getTransitions()) {
        // Look for the current transition's compound reference in the map.
        auto found = match_compound_refs_decoy_mz.find(tr.getCompoundRef());

        // Check if the compound reference is found and if the product m/z matches any in the set.
        if (found != match_compound_refs_decoy_mz.end() && found->second.count(tr.getProductMZ()) > 0) {
            // Create a new transition object based on the current transition.
            ReactionMonitoringTransition new_tr = tr;

            // Modify the product m/z of the new transition.
            new_tr.setProductMZ(tr.getProductMZ() + mass_to_add);

            // Add the updated transition to the new vector.
            v_rmt_new.push_back(std::move(new_tr));
        } else {
            // If no match is found, add the original transition to the new vector.
            v_rmt_new.push_back(tr);
        }
    }

    // Update the experiment's transitions with the new vector of updated transitions.
    t_exp.setTransitions(std::move(v_rmt_new));
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
