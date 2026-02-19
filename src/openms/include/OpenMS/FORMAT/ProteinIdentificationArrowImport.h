// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#ifdef WITH_PARQUET

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

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
  @brief Import ProteinIdentification data from Apache Arrow format

  This class provides static methods to import ProteinIdentification
  data from Apache Arrow Tables and Parquet files. Separate methods are
  provided for protein hits, protein groups, and search parameters.

  This is the reverse of ProteinIdentificationArrowExport and enables
  full round-trip serialization.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI ProteinIdentificationArrowImport
{
public:
  /**
    @brief Import all three Parquet files and reconstruct ProteinIdentifications

    Reads search_params, proteins, and protein_groups files in order,
    combining them into a single vector of ProteinIdentification objects.

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

#endif // WITH_PARQUET
