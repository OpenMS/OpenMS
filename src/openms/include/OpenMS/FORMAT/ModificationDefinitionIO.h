// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <map>
#include <set>
#include <string>
#include <vector>

namespace OpenMS
{
  class ResidueModification;

  /**
    @brief Collect, serialise and re-register the definitions of tool-defined modifications.

    A modification that is not in the shipped vocabularies (a run-time adduct, a custom search space)
    can be written into an identification file only by name, and no other process can resolve that
    name. This class carries the definitions alongside the identifications, as the
    Constants::UserParam::MODIFICATION_DEFINITIONS meta value on ProteinIdentification::SearchParameters
    - the one structure that idXML, featureXML, consensusXML and .idparquet all already serialise, so no
    format needs a schema change.

    Writers call collect() + attach() (or encode()); readers call registerFrom() as soon as the
    SearchParameters are available and BEFORE any peptide sequence is parsed, because
    AASequence::fromString throws on an unknown name.

    What is serialised is decided per modification by isDefinition(): provenance DEFINED and a
    non-empty Id. It is deliberately not "everything in ModificationsDB beyond the startup baseline":
    the database is process-global and accumulates every definition a process has ever loaded, so a
    tool reading many files would otherwise attribute the first file's definitions to an output that
    never references them.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ModificationDefinitionIO
  {
  public:
    /// The selection predicate: a named modification whose definition must travel with the file.
    static bool isDefinition(const ResidueModification* mod);

    /**
      @brief Walk the identifications and collect, per run identifier, the definitions they reference.

      Looks at every hit's N-terminal, residue and C-terminal modifications, plus the names in each
      run's fixed_modifications / variable_modifications (so a defined variable modification that
      happens to occur in no hit keeps its definition). O(hits). Runs referencing no definition do not
      appear in the result.
    */
    static std::map<std::string, std::set<const ResidueModification*>> collect(
      const std::vector<ProteinIdentification>& protein_ids,
      const PeptideIdentificationList& peptide_ids);

    /// Serialise @p defs into one ';'-joined string. Deterministic: records are sorted.
    static std::string encode(const std::set<const ResidueModification*>& defs);

    /**
      @brief Store @p defs on @p sp under Constants::UserParam::MODIFICATION_DEFINITIONS.

      Unions with whatever the meta value already holds rather than overwriting it, so a filter tool
      that drops every hit of a defined modification does not silently drop its definition too.
      A no-op when @p defs is empty and nothing was stored before.
    */
    static void attach(ProteinIdentification::SearchParameters& sp,
                       const std::set<const ResidueModification*>& defs);

    /**
      @brief Register every record in @p records into ModificationsDB.

      A malformed record is logged and skipped - a bad definition degrades to "unresolvable name", the
      pre-existing behaviour, and must not abort the whole file. @return the number of records registered.
    */
    static Size registerFrom(const std::string& records);

    /// Convenience: registerFrom() on the meta value of @p sp, if present. @return records registered.
    static Size registerFrom(const ProteinIdentification::SearchParameters& sp);
  };
} // namespace OpenMS
