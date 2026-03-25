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
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>

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
};

} // namespace OpenMS
