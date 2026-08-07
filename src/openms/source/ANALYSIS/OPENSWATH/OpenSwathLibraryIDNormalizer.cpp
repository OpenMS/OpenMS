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
    std::vector<std::string> source_compound_ids;
    source_compound_ids.reserve(exp.compounds.size());

    std::unordered_set<std::string> seen_compound_ids;
    seen_compound_ids.reserve(exp.compounds.size());

    for (const auto& compound : exp.compounds)
    {
      const std::string source_id = compound.id;
      if (source_id.empty())
      {
        throwInvalidID_("Source compound ID must not be empty", source_id);
      }
      if (!seen_compound_ids.insert(source_id).second)
      {
        throwInvalidID_("Source compound IDs must be unique", source_id);
      }
      source_compound_ids.push_back(source_id);
    }

    std::sort(source_compound_ids.begin(), source_compound_ids.end());

    std::unordered_map<std::string, std::string> precursor_id_map;
    precursor_id_map.reserve(source_compound_ids.size());
    for (Size i = 0; i < source_compound_ids.size(); ++i)
    {
      precursor_id_map.emplace(source_compound_ids[i], StringUtils::toStr(static_cast<int64_t>(i)));
    }

    for (auto& compound : exp.compounds)
    {
      compound.id = precursor_id_map.at(std::string(compound.id));
    }

    for (Size i = 0; i < exp.transitions.size(); ++i)
    {
      auto& transition = exp.transitions[i];
      const std::string source_ref = transition.peptide_ref;
      if (source_ref.empty())
      {
        throwInvalidID_("Transition precursor reference must not be empty", source_ref);
      }

      const auto precursor_it = precursor_id_map.find(source_ref);
      if (precursor_it == precursor_id_map.end())
      {
        throwInvalidID_("Transition references an unknown source compound ID", source_ref);
      }

      transition.peptide_ref = precursor_it->second;
      transition.transition_name = StringUtils::toStr(static_cast<int64_t>(i));
    }

    validateCanonicalIDs(exp);
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
