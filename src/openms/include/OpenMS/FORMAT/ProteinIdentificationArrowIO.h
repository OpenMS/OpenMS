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
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>

#include <map>
#include <memory>
#include <vector>

// Forward declarations
namespace arrow
{
  class Table;
}

namespace OpenMS
{

/**
  @brief Import and export ProteinIdentification data to/from Apache Arrow format

  This class provides static methods to export and import ProteinIdentification
  data to/from Apache Arrow Tables and Parquet files. Separate tables are
  provided for protein hits, protein groups, and search parameters.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI ProteinIdentificationArrowIO
{
public:
  // ==================== Export methods ====================

  /**
    @brief Export protein hits to Apache Arrow Table

    Each ProteinHit becomes one row with identification, score, and metadata.

    @param[in] protein_identifications Vector of protein identifications
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportProteinsToArrow(
    const std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Export protein hits to Parquet file

    @param[in] protein_identifications Vector of protein identifications
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportProteinsToParquet(
    const std::vector<ProteinIdentification>& protein_identifications,
    const String& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Export protein groups to Apache Arrow Table

    Each ProteinGroup becomes one row with group probability and member accessions.

    @param[in] protein_identifications Vector of protein identifications
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportProteinGroupsToArrow(
    const std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Export protein groups to Parquet file

    @param[in] protein_identifications Vector of protein identifications
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportProteinGroupsToParquet(
    const std::vector<ProteinIdentification>& protein_identifications,
    const String& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Export search parameters to Apache Arrow Table

    Each ProteinIdentification's SearchParameters becomes one row.

    @param[in] protein_identifications Vector of protein identifications
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportSearchParamsToArrow(
    const std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Export search parameters to Parquet file

    @param[in] protein_identifications Vector of protein identifications
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportSearchParamsToParquet(
    const std::vector<ProteinIdentification>& protein_identifications,
    const String& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  // ==================== Import methods ====================

  /**
    @brief Import all three Parquet files and reconstruct ProteinIdentifications

    Reads the three Parquet files and reconstructs a vector of
    ProteinIdentification objects with hits, groups, and search parameters.

    @param[in] proteins_filename Path to proteins Parquet file
    @param[in] protein_groups_filename Path to protein groups Parquet file
    @param[in] search_params_filename Path to search params Parquet file
    @param[out] protein_identifications Reconstructed protein identifications
    @return true on success, false on error
  */
  static bool importFromParquet(
    const String& proteins_filename,
    const String& protein_groups_filename,
    const String& search_params_filename,
    std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Import search parameters from Arrow Table

    Each row becomes a ProteinIdentification shell with run-level metadata
    and SearchParameters populated.

    @param[in] table Arrow Table with search parameters
    @param[out] protein_identifications Reconstructed protein identifications
    @return true on success, false on error
  */
  static bool importSearchParamsFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Import protein hits from Arrow Table

    Adds ProteinHits to matching ProteinIdentifications by run_identifier.
    If no matching ProteinIdentification exists, creates new ones.

    @param[in] table Arrow Table with protein hits
    @param[out] protein_identifications Protein identifications to populate
    @return true on success, false on error
  */
  static bool importProteinsFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Import protein groups from Arrow Table

    Adds ProteinGroups and IndistinguishableProteins to matching
    ProteinIdentifications by run_identifier.

    @param[in] table Arrow Table with protein groups
    @param[out] protein_identifications Protein identifications to populate
    @return true on success, false on error
  */
  static bool importProteinGroupsFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Import search parameters from Parquet file

    @param[in] filename Path to Parquet file
    @param[out] protein_identifications Reconstructed protein identifications
    @return true on success, false on error
  */
  static bool importSearchParamsFromParquet(
    const String& filename,
    std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Import protein hits from Parquet file

    @param[in] filename Path to Parquet file
    @param[out] protein_identifications Protein identifications to populate
    @return true on success, false on error
  */
  static bool importProteinsFromParquet(
    const String& filename,
    std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Import protein groups from Parquet file

    @param[in] filename Path to Parquet file
    @param[out] protein_identifications Protein identifications to populate
    @return true on success, false on error
  */
  static bool importProteinGroupsFromParquet(
    const String& filename,
    std::vector<ProteinIdentification>& protein_identifications);

  // ==================== Identifier handling parity with XML lane ====================

  /**
    @brief Synthesize fresh run identifiers per ProteinIdentification, mirroring IdXMLFile.

    Mirrors IdXMLFile.cpp:530: every load assigns each ProteinIdentification a fresh
    identifier `<search_engine>_<date>_<UniqueIdGenerator>`. The stored identifier
    on disk is informational; the in-memory identifier is regenerated. This is the
    same defense FeatureXMLHandler / ConsensusXMLHandler / IdXMLFile apply on load
    against downstream-collision-after-rip-and-merge scenarios.

    The function mutates @p protein_identifications in place and returns the map
    { stored_id -> synthesized_id } so the caller can apply the same rename to
    each PeptideIdentification collection it owns (FeatureMap has 2: per-feature
    and unassigned; ConsensusMap has 2; PSMArrowIO has 1).

    Edge cases:
      - empty getSearchEngine() falls back to literal "unknown"
      - invalid getDateTime() uses placeholder "1900-01-01T00:00:00"
      - multiple ProtIDs sharing one stored identifier each receive their own
        distinct synthesized identifier; the returned map collapses to the
        last-seen entry. An OPENMS_LOG_WARN is emitted once per such collision.

    @param[in,out] protein_identifications ProtID vector whose identifiers will be replaced
    @return Map from each stored identifier to its synthesized replacement.
  */
  static std::map<String, String> synthesizeRunIdentifiers(
    std::vector<ProteinIdentification>& protein_identifications);

  /**
    @brief Apply a stored->synthesized identifier rename to a PeptideIdentification collection.

    PeptideIdentifications whose stored identifier isn't present in @p rename are
    left untouched. Mirrors the orphan-pep_id semantics of the XML lane (where
    pep_ids that don't match any ProtID retain their stale identifier).

    @param[in] rename Map produced by synthesizeRunIdentifiers
    @param[in,out] pep_ids PeptideIdentification collection to re-stamp
  */
  static void applyRunIdentifierRename(
    const std::map<String, String>& rename,
    PeptideIdentificationList& pep_ids);

  /**
    @brief Reject a ProteinIdentification vector with duplicate identifiers (store-side check).

    Mirrors XMLHandler::checkUniqueIdentifiers_ — throws Exception::InvalidValue
    with the same message text used by the XML lane's fatalError. Called by every
    parquet store entry point before any Arrow builder is allocated, so no
    partial file is created on rejection.

    @param[in] protein_identifications ProtID vector to check
    @throw Exception::InvalidValue when duplicate identifiers are present
  */
  static void checkUniqueIdentifiers(
    const std::vector<ProteinIdentification>& protein_identifications);
};

} // namespace OpenMS
