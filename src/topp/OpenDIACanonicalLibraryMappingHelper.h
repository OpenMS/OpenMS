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
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <cstdint>
#include <string>
#include <unordered_map>

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

    // IDs are already canonical dense integer strings. Preserve that exact ID
    // domain because Parquet writers persist compound.id directly.
    mapping.compound_to_precursor.reserve(targeted_exp.compounds.size());
    for (const auto& compound : targeted_exp.compounds)
    {
      mapping.compound_to_precursor.emplace(compound.id, StringUtils::toInt64(compound.id));
    }

    mapping.precursor_mz_by_id.reserve(targeted_exp.compounds.size());
    mapping.precursor_decoy_by_id.reserve(targeted_exp.compounds.size());
    mapping.transition_to_id.reserve(targeted_exp.transitions.size());
    for (const auto& transition : targeted_exp.transitions)
    {
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

      // transition_name is canonicalized to the integer transition ID.
      mapping.transition_to_id.try_emplace(transition.transition_name,
                                           StringUtils::toInt64(transition.transition_name));
    }

    return mapping;
  }
} // namespace OpenMS::Internal
