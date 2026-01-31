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

    for (const auto& prot_id : prot_ids)
    {
      StringList ms_run_paths;
      prot_id.getPrimaryMSRunPath(ms_run_paths);
      const String& identifier = prot_id.getIdentifier();

      identifier_to_msrunpath_[identifier] = ms_run_paths;

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

  String IdentifierMSRunMapper::getPrimaryMSRunPath(const PeptideIdentification& pepid) const
  {
    const String& identifier = pepid.getIdentifier();
    auto it = identifier_to_msrunpath_.find(identifier);
    if (it == identifier_to_msrunpath_.end())
    {
      return String();
    }

    const StringList& ms_run_paths = it->second;
    if (ms_run_paths.empty())
    {
      return String();
    }

    // Determine which file to use (default to index 0)
    Size merge_index = 0;
    if (pepid.metaValueExists(Constants::UserParam::ID_MERGE_INDEX))
    {
      merge_index = static_cast<Size>(pepid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX));
    }

    // Check if index is valid
    if (merge_index >= ms_run_paths.size())
    {
      return String(); // Invalid index
    }

    return ms_run_paths[merge_index];
  }

  bool IdentifierMSRunMapper::hasIdentifier(const String& identifier) const
  {
    return identifier_to_msrunpath_.find(identifier) != identifier_to_msrunpath_.end();
  }

  const String& IdentifierMSRunMapper::getIdentifier(const StringList& ms_run_paths) const
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

  const StringList& IdentifierMSRunMapper::getMSRunPaths(const String& identifier) const
  {
    auto it = identifier_to_msrunpath_.find(identifier);
    if (it == identifier_to_msrunpath_.end())
    {
      return empty_stringlist_;
    }
    return it->second;
  }

  std::vector<String> IdentifierMSRunMapper::getIdentifiers() const
  {
    std::vector<String> identifiers;
    identifiers.reserve(identifier_to_msrunpath_.size());
    for (const auto& pair : identifier_to_msrunpath_)
    {
      identifiers.push_back(pair.first);
    }
    return identifiers;
  }

  bool IdentifierMSRunMapper::hasRunPath(const StringList& ms_run_paths) const
  {
    return runpath_to_identifier_.find(ms_run_paths) != runpath_to_identifier_.end();
  }

  bool IdentifierMSRunMapper::tryGetIdentifier(const StringList& ms_run_paths, String& identifier) const
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
