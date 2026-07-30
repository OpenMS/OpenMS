// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/IdentifierMSRunMapper.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
  // Static member initialization
  const StringList IdentifierMSRunMapper::empty_stringlist_;

  IdentifierMSRunMapper::IdentifierMSRunMapper(const std::vector<ProteinIdentification>& prot_ids)
  {
    create(prot_ids);
  }

  void IdentifierMSRunMapper::create(const std::vector<ProteinIdentification>& prot_ids)
  {
    identifier_to_msrunpath_.clear();
    runpath_to_identifier_.clear();

    // First pass: build the forward map (identifier -> ms-run-paths). This direction is
    // always unambiguous - even when several runs share the same path - and is all that
    // getPrimaryMSRunPath() needs. Building it completely up front means that if the
    // reverse-map duplicate check below throws, callers that catch the exception still
    // see a fully-populated forward map (i.e. USI resolution degrades cleanly, not in an
    // order-dependent way).
    for (const auto& prot_id : prot_ids)
    {
      StringList ms_run_paths;
      prot_id.getPrimaryMSRunPath(ms_run_paths);
      identifier_to_msrunpath_[prot_id.getIdentifier()] = ms_run_paths;
    }

    // Second pass: build the reverse map (ms-run-paths -> identifier), rejecting ambiguous
    // input where different identifiers map to the same paths.
    for (const auto& prot_id : prot_ids)
    {
      StringList ms_run_paths;
      prot_id.getPrimaryMSRunPath(ms_run_paths);
      const std::string& identifier = prot_id.getIdentifier();

      // Check for duplicate ms_run_paths (different identifiers mapping to same paths)
      const auto it = runpath_to_identifier_.find(ms_run_paths);
      if (it != runpath_to_identifier_.end())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Multiple protein identifications with the same ms-run-path. Check input!",
          ListUtils::concatenate(ms_run_paths, ","));
      }
      runpath_to_identifier_[ms_run_paths] = identifier;
    }
  }

  std::string IdentifierMSRunMapper::getPrimaryMSRunPath(const PeptideIdentification& pepid) const
  {
    const std::string& identifier = pepid.getIdentifier();
    auto it = identifier_to_msrunpath_.find(identifier);
    if (it != identifier_to_msrunpath_.end() && !it->second.empty())
    {
      const StringList& ms_run_paths = it->second;

      // Determine which file to use (default to index 0)
      Size merge_index = 0;
      if (pepid.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
      {
        merge_index = static_cast<Size>(pepid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX));
      }

      if (merge_index < ms_run_paths.size())
      {
        return ms_run_paths[merge_index];
      }
    }

    // Legacy fallback: data read from older idXML files or other formats may annotate
    // the source file directly on the PeptideIdentification via the (deprecated) "base_name"
    // meta value, instead of via the ProteinIdentification's primary MS run path. Honor it so
    // that such files keep resolving to a source run after the removal of base_name accessors.
    if (pepid.metaValueExists(Constants::UserParam::BASE_NAME))
    {
      return StringUtils::toStr(pepid.getMetaValue(Constants::UserParam::BASE_NAME));
    }

    return std::string();
  }

  void IdentifierMSRunMapper::validateMergeIndex(const PeptideIdentification& pep_id, size_t psm_index) const
  {
    const StringList& paths = getMSRunPaths(pep_id.getIdentifier());
    if (paths.size() < 2) { return; }

    const std::string key = std::string(Constants::UserParam::ID_MERGE_INDEX);
    const std::string where = "PSM #" + StringUtils::toStr(psm_index) + " of identification run '"
                            + pep_id.getIdentifier() + "' (" + StringUtils::toStr(paths.size())
                            + " files in the run)";

    if (!pep_id.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Multiple files in a run, but no '" + key + "' in PeptideIdentification found: " + where + ".");
    }

    const DataValue& dv = pep_id.getMetaValue(Constants::UserParam::ID_MERGE_INDEX);
    if (dv.valueType() != DataValue::INT_VALUE)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "'" + key + "' is not an integer: " + where + ".");
    }

    // long long, not Size: DataValue's unsigned conversion throws its own ConversionError on a
    // negative value, which would escape as an undiagnosed exception instead of the message above.
    const long long merge_index = static_cast<long long>(dv);
    if (merge_index < 0 || static_cast<size_t>(merge_index) >= paths.size())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "'" + key + "' = " + StringUtils::toStr(merge_index) + " is out of range: " + where + ".");
    }
  }

  bool IdentifierMSRunMapper::hasIdentifier(const std::string& identifier) const
  {
    return identifier_to_msrunpath_.contains(identifier);
  }

  const std::string& IdentifierMSRunMapper::getIdentifier(const StringList& ms_run_paths) const
  {
    auto it = runpath_to_identifier_.find(ms_run_paths);
    if (it == runpath_to_identifier_.end())
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "MS run paths not found in mapping");
    }
    return it->second;
  }

  bool IdentifierMSRunMapper::empty() const
  {
    return identifier_to_msrunpath_.empty();
  }

  Size IdentifierMSRunMapper::size() const
  {
    return identifier_to_msrunpath_.size();
  }

  const StringList& IdentifierMSRunMapper::getMSRunPaths(const std::string& identifier) const
  {
    auto it = identifier_to_msrunpath_.find(identifier);
    if (it == identifier_to_msrunpath_.end())
    {
      return empty_stringlist_;
    }
    return it->second;
  }

  std::vector<std::string> IdentifierMSRunMapper::getIdentifiers() const
  {
    std::vector<std::string> identifiers;
    identifiers.reserve(identifier_to_msrunpath_.size());
    for (const auto& pair : identifier_to_msrunpath_)
    {
      identifiers.push_back(pair.first);
    }
    return identifiers;
  }

  bool IdentifierMSRunMapper::hasRunPath(const StringList& ms_run_paths) const
  {
    return runpath_to_identifier_.contains(ms_run_paths);
  }

  bool IdentifierMSRunMapper::tryGetIdentifier(const StringList& ms_run_paths, std::string& identifier) const
  {
    auto it = runpath_to_identifier_.find(ms_run_paths);
    if (it == runpath_to_identifier_.end())
    {
      return false;
    }
    identifier = it->second;
    return true;
  }

} // namespace OpenMS
