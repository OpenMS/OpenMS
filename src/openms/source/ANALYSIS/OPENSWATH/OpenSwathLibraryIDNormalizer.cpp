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
#include <charconv>
#include <cstdint>
#include <string>
#include <system_error>
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

    bool hasConventionalDecoyPrefix_(const std::string& id)
    {
      return id.rfind("DECOY", 0) == 0 || id.rfind("Decoy", 0) == 0 || id.rfind("decoy", 0) == 0;
    }

    bool hasCanonicalIDFormat_(const std::string& text) noexcept
    {
      if (text.empty() || (text.size() > 1 && text.front() == '0'))
      {
        return false;
      }

      int64_t value = -1;
      const char* begin = text.data();
      const char* end = begin + text.size();
      const auto [ptr, ec] = std::from_chars(begin, end, value);
      return ec == std::errc{} && ptr == end && value >= 0;
    }

    void validateCanonicalIDFormat_(const std::string& text, const std::string& entity)
    {
      if (text.empty())
      {
        throwInvalidID_(entity + " ID must not be empty", text);
      }
      if (!hasCanonicalIDFormat_(text))
      {
        throwInvalidID_(entity + " ID must be a non-negative Int64 in canonical decimal form", text);
      }
    }
  } // namespace

  OpenSwathLibraryIDNormalizer::SourceIDMapping OpenSwathLibraryIDNormalizer::normalizeSourceIDs(OpenSwath::LightTargetedExperiment& exp)
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

    // Assign canonical precursor IDs from the complete source compound set before constructing
    // the operational subset. This mirrors the persistent PQP convention: an unreferenced source
    // compound may later be omitted from the LightTargetedExperiment, but its assigned ID remains
    // reserved so surviving precursor IDs are stable across direct source loading and source->PQP
    // round-trips. Sparse canonical precursor IDs are therefore expected and valid.
    std::vector<std::string> source_compound_ids;
    source_compound_ids.reserve(exp.compounds.size());
    for (const auto& compound : exp.compounds)
    {
      source_compound_ids.emplace_back(compound.id);
    }
    std::sort(source_compound_ids.begin(), source_compound_ids.end());

    std::unordered_map<std::string, std::string> precursor_id_map;
    precursor_id_map.reserve(source_compound_ids.size());
    for (Size i = 0; i < source_compound_ids.size(); ++i)
    {
      precursor_id_map.emplace(source_compound_ids[i], StringUtils::toStr(static_cast<int64_t>(i)));
    }

    SourceIDMapping source_ids;
    source_ids.precursor_source_to_canonical = precursor_id_map;
    source_ids.precursor_canonical_to_source.reserve(precursor_id_map.size());
    for (const auto& [source_id, canonical_id] : precursor_id_map)
    {
      source_ids.precursor_canonical_to_source.emplace(canonical_id, source_id);
    }

    // Only compounds referenced by transitions are retained in the operational OpenSWATH
    // experiment. Their IDs are not compressed after filtering; they keep the values assigned
    // from the complete source compound set above.
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
      if (!referenced_compound_ids.contains(source_id))
      {
        continue;
      }

      const auto precursor_it = precursor_id_map.find(source_id);
      if (precursor_it == precursor_id_map.end())
      {
        // All source compounds were inserted into precursor_id_map above.
        throwInvalidID_("Source compound is missing from the canonical precursor map", source_id);
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

      const std::string source_transition_id = transition.transition_name;
      if (!transition.getDecoy() &&
          (hasConventionalDecoyPrefix_(source_ref) || hasConventionalDecoyPrefix_(source_transition_id)))
      {
        transition.setDecoy(true);
      }

      const std::string canonical_transition_id = StringUtils::toStr(static_cast<int64_t>(i));
      source_ids.transition_canonical_to_source.emplace(canonical_transition_id, source_transition_id);
      transition.peptide_ref = precursor_it->second;
      transition.transition_name = canonical_transition_id;
    }

    validateCanonicalIDs(normalized);
    exp = std::move(normalized);
    return source_ids;
  }

  void OpenSwathLibraryIDNormalizer::materializeDecoyPrefix(OpenSwath::LightTargetedExperiment& exp,
                                                             const SourceIDMapping& source_ids,
                                                             const std::string& decoy_prefix)
  {
    if (decoy_prefix.empty())
    {
      return;
    }

    std::unordered_set<std::string> decoy_precursor_ids;
    decoy_precursor_ids.reserve(source_ids.precursor_source_to_canonical.size());
    for (const auto& [source_id, canonical_id] : source_ids.precursor_source_to_canonical)
    {
      if (source_id.rfind(decoy_prefix, 0) == 0)
      {
        decoy_precursor_ids.insert(canonical_id);
      }
    }

    for (auto& transition : exp.transitions)
    {
      bool source_transition_is_decoy = false;
      const auto transition_source_it = source_ids.transition_canonical_to_source.find(transition.transition_name);
      if (transition_source_it != source_ids.transition_canonical_to_source.end())
      {
        source_transition_is_decoy = transition_source_it->second.rfind(decoy_prefix, 0) == 0;
      }

      if (decoy_precursor_ids.contains(transition.peptide_ref) || source_transition_is_decoy)
      {
        transition.setDecoy(true);
      }
    }
  }

  std::optional<std::string> OpenSwathLibraryIDNormalizer::canonicalTargetForDecoyPrecursor(
    const std::string& decoy_ref,
    const SourceIDMapping& source_ids,
    const std::string& decoy_prefix)
  {
    if (decoy_prefix.empty())
    {
      return std::nullopt;
    }

    // Canonical path: recover the source precursor ID, strip the configured prefix
    // there, then map the target source ID back into the canonical operational domain.
    const auto source_it = source_ids.precursor_canonical_to_source.find(decoy_ref);
    if (source_it != source_ids.precursor_canonical_to_source.end())
    {
      const std::string& source_decoy_id = source_it->second;
      if (source_decoy_id.rfind(decoy_prefix, 0) != 0)
      {
        return std::nullopt;
      }

      const std::string target_source_id = source_decoy_id.substr(decoy_prefix.size());
      const auto target_it = source_ids.precursor_source_to_canonical.find(target_source_id);
      if (target_it == source_ids.precursor_source_to_canonical.end())
      {
        return std::nullopt;
      }
      return target_it->second;
    }

    // Compatibility for callers that still operate directly on source IDs and therefore
    // do not have a provenance map. Do not use this fallback when a mapping is present,
    // because a numeric-looking source ID could otherwise collide with a canonical ID.
    if (source_ids.precursor_canonical_to_source.empty() &&
        source_ids.precursor_source_to_canonical.empty() &&
        decoy_ref.rfind(decoy_prefix, 0) == 0)
    {
      return decoy_ref.substr(decoy_prefix.size());
    }

    return std::nullopt;
  }

  bool OpenSwathLibraryIDNormalizer::hasCanonicalIDFormat(const OpenSwath::LightTargetedExperiment& exp) noexcept
  {
    for (const auto& compound : exp.compounds)
    {
      if (!hasCanonicalIDFormat_(compound.id))
      {
        return false;
      }
    }

    for (const auto& transition : exp.transitions)
    {
      if (!hasCanonicalIDFormat_(transition.transition_name) ||
          !hasCanonicalIDFormat_(transition.peptide_ref))
      {
        return false;
      }
    }
    return true;
  }

  bool OpenSwathLibraryIDNormalizer::hasCanonicalIDs(const OpenSwath::LightTargetedExperiment& exp)
  {
    if (!hasCanonicalIDFormat(exp))
    {
      return false;
    }

    std::unordered_set<std::string> precursor_ids;
    precursor_ids.reserve(exp.compounds.size());
    for (const auto& compound : exp.compounds)
    {
      if (!precursor_ids.insert(compound.id).second)
      {
        return false;
      }
    }

    std::unordered_set<std::string> transition_ids;
    transition_ids.reserve(exp.transitions.size());
    for (const auto& transition : exp.transitions)
    {
      if (!transition_ids.insert(transition.transition_name).second ||
          !precursor_ids.contains(transition.peptide_ref))
      {
        return false;
      }
    }
    return true;
  }

  void OpenSwathLibraryIDNormalizer::validateCanonicalIDs(const OpenSwath::LightTargetedExperiment& exp)
  {
    std::unordered_set<std::string> precursor_ids;
    precursor_ids.reserve(exp.compounds.size());

    for (const auto& compound : exp.compounds)
    {
      const std::string id = compound.id;
      validateCanonicalIDFormat_(id, "Precursor");
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
      validateCanonicalIDFormat_(transition_id, "Transition");
      if (!transition_ids.insert(transition_id).second)
      {
        throwInvalidID_("Canonical transition IDs must be unique", transition_id);
      }

      const std::string precursor_ref = transition.peptide_ref;
      validateCanonicalIDFormat_(precursor_ref, "Transition precursor reference");
      if (!precursor_ids.contains(precursor_ref))
      {
        throwInvalidID_("Transition references an unknown canonical precursor ID", precursor_ref);
      }
    }
  }
} // namespace OpenMS
