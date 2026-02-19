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
  @brief Export ProteinIdentification data to Apache Arrow format

  This class provides static methods to export ProteinIdentification
  data to Apache Arrow Tables and Parquet files. Separate tables are
  provided for protein hits, protein groups, and search parameters.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI ProteinIdentificationArrowExport
{
public:
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
};

} // namespace OpenMS

#endif // WITH_PARQUET
