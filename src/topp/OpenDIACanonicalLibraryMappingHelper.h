// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <algorithm>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS::Internal
{
  struct OpenSwathCanonicalLibraryMapping
  {
    std::unordered_map<std::string, int64_t> compound_to_precursor;
    std::unordered_map<int64_t, double> precursor_mz_by_id;
    std::unordered_map<int64_t, bool> precursor_decoy_by_id;
    std::unordered_map<std::string, int64_t> transition_to_id;
  };

  inline OpenSwathCanonicalLibraryMapping buildOpenSwathCanonicalLibraryMapping(
    const OpenSwath::LightTargetedExperiment& targeted_exp)
  {
    OpenSwathCanonicalLibraryMapping mapping;

    std::vector<std::string> group_ids;
    group_ids.reserve(targeted_exp.compounds.size());
    for (const auto& compound : targeted_exp.compounds)
    {
      group_ids.push_back(compound.id);
    }
    std::sort(group_ids.begin(), group_ids.end());
    group_ids.erase(std::unique(group_ids.begin(), group_ids.end()), group_ids.end());

    mapping.compound_to_precursor.reserve(group_ids.size());
    for (Size i = 0; i < group_ids.size(); ++i)
    {
      mapping.compound_to_precursor.emplace(group_ids[i], static_cast<int64_t>(i));
    }

    mapping.precursor_mz_by_id.reserve(group_ids.size());
    mapping.precursor_decoy_by_id.reserve(group_ids.size());
    mapping.transition_to_id.reserve(targeted_exp.transitions.size());
    for (Size i = 0; i < targeted_exp.transitions.size(); ++i)
    {
      const auto& transition = targeted_exp.transitions[i];
      const auto precursor_it = mapping.compound_to_precursor.find(transition.peptide_ref);
      if (precursor_it == mapping.compound_to_precursor.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Transition references unknown peptide_ref '" +
                                            std::string(transition.peptide_ref) + "'");
      }

      const int64_t precursor_id = precursor_it->second;
      if (!mapping.precursor_mz_by_id.contains(precursor_id))
      {
        mapping.precursor_mz_by_id.emplace(precursor_id, transition.precursor_mz);
      }
      if (!mapping.precursor_decoy_by_id.contains(precursor_id) && transition.isDetectingTransition())
      {
        mapping.precursor_decoy_by_id.emplace(precursor_id, transition.getDecoy());
      }

      // Retain the first matching transition name to mirror the sqlite writer's
      // stable input-order assignment as closely as possible.
      mapping.transition_to_id.try_emplace(transition.transition_name, static_cast<int64_t>(i));
    }

    return mapping;
  }
} // namespace OpenMS::Internal
