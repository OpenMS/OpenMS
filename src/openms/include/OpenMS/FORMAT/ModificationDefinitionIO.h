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
  class ConsensusMap;
  class FeatureMap;
  class ResidueModification;

  /**
    @brief Collect, serialise and re-register the definitions of tool-defined modifications.

    The definitions travel as the Constants::UserParam::MODIFICATION_DEFINITIONS meta value on
    ProteinIdentification::SearchParameters, which idXML, featureXML, consensusXML and .idparquet all
    serialise. Writers call collect() + attach(); readers call registerFrom() before parsing any
    peptide sequence. isDefinition() selects what is serialised: provenance DEFINED and a non-empty Id.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ModificationDefinitionIO
  {
  public:
    /// The selection predicate: a named modification whose definition must travel with the file.
    static bool isDefinition(const ResidueModification* mod);

    /**
      @brief Collect, per run identifier, the definitions the identifications reference.

      Covers every hit's terminal and residue modifications plus each run's fixed/variable modification
      names. Runs referencing no definition do not appear.
    */
    static std::map<std::string, std::set<const ResidueModification*>> collect(
      const std::vector<ProteinIdentification>& protein_ids,
      const PeptideIdentificationList& peptide_ids);

    /// collect() over a FeatureMap: its runs, the identifications of every feature (subordinates included) and the unassigned ones
    static std::map<std::string, std::set<const ResidueModification*>> collect(const FeatureMap& map);

    /// collect() over a ConsensusMap: its runs, the identifications of every consensus feature and the unassigned ones
    static std::map<std::string, std::set<const ResidueModification*>> collect(const ConsensusMap& map);

    /// Serialise @p defs into one ';'-joined string. Deterministic: records are sorted.
    static std::string encode(const std::set<const ResidueModification*>& defs);

    /// Store @p defs on @p sp under Constants::UserParam::MODIFICATION_DEFINITIONS, unioned with
    /// what the meta value already holds. No-op when there is nothing to store.
    static void attach(ProteinIdentification::SearchParameters& sp,
                       const std::set<const ResidueModification*>& defs);

    /// Per run identifier, the value attach() would store: @p defs unioned with what the run's
    /// SearchParameters already carry. Runs with nothing to store are absent.
    static std::map<std::string, std::string> encodeByRun(
      const std::vector<ProteinIdentification>& protein_ids,
      const std::map<std::string, std::set<const ResidueModification*>>& defs);

    /// Register every record in @p records into ModificationsDB; a malformed record is logged and
    /// skipped. @return the number of records registered
    static Size registerFrom(const std::string& records);

    /// Convenience: registerFrom() on the meta value of @p sp, if present. @return records registered.
    static Size registerFrom(const ProteinIdentification::SearchParameters& sp);
  };
} // namespace OpenMS
