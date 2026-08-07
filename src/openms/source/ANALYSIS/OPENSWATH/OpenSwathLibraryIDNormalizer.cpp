// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryIDNormalizer.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <algorithm>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace OpenMS
{
  namespace
  {
    [[noreturn]] void throwInvalidID_(const std::string& message, const std::string& value)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message, value);
    }

    int64_t parseCanonicalID_(const std::string& text, const std::string& entity)
    {
      if (text.empty())
      {
        throwInvalidID_(entity + " ID must not be empty", text);
      }

      int64_t value = -1;
      try
      {
        value = StringUtils::toInt64(text);
      }
      catch (Exception::ConversionError&)
      {
        throwInvalidID_(entity + " ID must be a non-negative Int64 in canonical decimal form", text);
      }

      if (value < 0 || StringUtils::toStr(value) != text)
      {
        throwInvalidID_(entity + " ID must be a non-negative Int64 in canonical decimal form", text);
      }
      return value;
    }
  } // namespace

  void OpenSwathLibraryIDNormalizer::normalizeSourceIDs(OpenSwath::LightTargetedExperiment& exp)
  {
    // Validate source compound identifiers before constructing the operational subset.
    std::unordered_set<std::string> known_compound_ids;
    known_compound_ids.reserve(exp.compounds.size());
    for (const auto& compound : exp.compounds)
    {
      const std::string source_id = compound.id;
      if (source_id.empty())
      {
        throwInvalidID_("Source compound ID must not be empty", source_id);
      }
      if (!known_compound_ids.insert(source_id).second)
      {
        throwInvalidID_("Source compound IDs must be unique", source_id);
      }
    }

    // Only compounds referenced by transitions participate in the operational OpenSWATH
    // precursor ID space. Unreferenced source records cannot provide an extraction coordinate
    // and must not shift the canonical IDs of extractable precursor groups.
    std::unordered_set<std::string> referenced_compound_ids;
    referenced_compound_ids.reserve(exp.transitions.size());
    for (const auto& transition : exp.transitions)
    {
      const std::string source_ref = transition.peptide_ref;
      if (source_ref.empty())
      {
        throwInvalidID_("Transition precursor reference must not be empty", source_ref);
      }
      if (!known_compound_ids.contains(source_ref))
      {
        throwInvalidID_("Transition references an unknown source compound ID", source_ref);
      }
      referenced_compound_ids.insert(source_ref);
    }

    std::vector<std::string> source_compound_ids;
    source_compound_ids.reserve(referenced_compound_ids.size());
    for (const auto& source_id : referenced_compound_ids)
    {
      source_compound_ids.push_back(source_id);
    }
    std::sort(source_compound_ids.begin(), source_compound_ids.end());

    std::unordered_map<std::string, std::string> precursor_id_map;
    precursor_id_map.reserve(source_compound_ids.size());
    for (Size i = 0; i < source_compound_ids.size(); ++i)
    {
      precursor_id_map.emplace(source_compound_ids[i], StringUtils::toStr(static_cast<int64_t>(i)));
    }

    // Build a fresh experiment instead of mutating compound IDs in place. LightTargetedExperiment
    // maintains an internal compound-reference lookup cache; rebuilding the object guarantees that
    // no cache entries keyed by pre-normalization source IDs survive canonicalization.
    OpenSwath::LightTargetedExperiment normalized;
    normalized.proteins = exp.proteins;
    normalized.compounds.reserve(referenced_compound_ids.size());
    normalized.transitions = exp.transitions;

    for (const auto& source_compound : exp.compounds)
    {
      const std::string source_id = source_compound.id;
      const auto precursor_it = precursor_id_map.find(source_id);
      if (precursor_it == precursor_id_map.end())
      {
        continue;
      }

      OpenSwath::LightCompound compound = source_compound;
      compound.id = precursor_it->second;
      normalized.compounds.push_back(std::move(compound));
    }

    for (Size i = 0; i < normalized.transitions.size(); ++i)
    {
      auto& transition = normalized.transitions[i];
      const std::string source_ref = transition.peptide_ref;
      const auto precursor_it = precursor_id_map.find(source_ref);
      if (precursor_it == precursor_id_map.end())
      {
        // The source-reference validation above makes this unreachable unless the experiment
        // is modified concurrently while normalization is running.
        throwInvalidID_("Transition references an unknown source compound ID", source_ref);
      }

      transition.peptide_ref = precursor_it->second;
      transition.transition_name = StringUtils::toStr(static_cast<int64_t>(i));
    }

    validateCanonicalIDs(normalized);
    exp = std::move(normalized);
  }

  void OpenSwathLibraryIDNormalizer::validateCanonicalIDs(const OpenSwath::LightTargetedExperiment& exp)
  {
    std::unordered_set<std::string> precursor_ids;
    precursor_ids.reserve(exp.compounds.size());

    for (const auto& compound : exp.compounds)
    {
      const std::string id = compound.id;
      parseCanonicalID_(id, "Precursor");
      if (!precursor_ids.insert(id).second)
      {
        throwInvalidID_("Canonical precursor IDs must be unique", id);
      }
    }

    std::unordered_set<std::string> transition_ids;
    transition_ids.reserve(exp.transitions.size());

    for (const auto& transition : exp.transitions)
    {
      const std::string transition_id = transition.transition_name;
      parseCanonicalID_(transition_id, "Transition");
      if (!transition_ids.insert(transition_id).second)
      {
        throwInvalidID_("Canonical transition IDs must be unique", transition_id);
      }

      const std::string precursor_ref = transition.peptide_ref;
      parseCanonicalID_(precursor_ref, "Transition precursor reference");
      if (!precursor_ids.contains(precursor_ref))
      {
        throwInvalidID_("Transition references an unknown canonical precursor ID", precursor_ref);
      }
    }
  }
} // namespace OpenMS
