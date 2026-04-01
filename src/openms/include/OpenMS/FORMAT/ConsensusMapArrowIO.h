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
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>

#include <memory>

// Forward declarations
namespace arrow
{
  class Table;
}

namespace OpenMS
{

/**
  @brief Import and export ConsensusMap data to/from Apache Arrow format

  This class provides static methods to export and import ConsensusMap
  data to/from Apache Arrow Tables and Parquet files. Separate tables are
  provided for consensus features (with their handles and metadata) and for
  peptide spectrum matches (PSMs) associated with features.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI ConsensusMapArrowIO
{
public:
  // ==================== Export methods ====================

  /**
    @brief Export consensus features to Apache Arrow Table

    Each ConsensusFeature becomes one row with RT, MZ, intensity, charge,
    quality, width, nested FeatureHandles, and metadata columns.

    @param[in] cmap The ConsensusMap to export
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportFeaturesToArrow(
    const ConsensusMap& cmap);

  /**
    @brief Export peptide spectrum matches (PSMs) associated with consensus features to Apache Arrow Table

    Each PeptideHit from each PeptideIdentification (both feature-level
    and unassigned) becomes one row.

    @param[in] cmap The ConsensusMap whose identifications to export
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportPSMsToArrow(
    const ConsensusMap& cmap);

  /**
    @brief Export ConsensusMap to a directory of Parquet files

    Writes five Parquet files: consensus_features.parquet, psms.parquet,
    proteins.parquet, protein_groups.parquet, and search_params.parquet
    into the specified directory. Protein-level data is delegated to
    ProteinIdentificationArrowIO. ConsensusMap-level metadata
    (column headers, experiment type, DocumentIdentifier, DataProcessing)
    is stored as file-level key-value metadata in consensus_features.parquet.

    @param[in] cmap The ConsensusMap to export
    @param[in] directory Output directory path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const ConsensusMap& cmap,
    const String& directory,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  // ==================== Import methods ====================

  /**
    @brief Import consensus features from Apache Arrow Table

    Each row becomes a ConsensusFeature with RT, MZ, intensity, charge,
    quality, width, FeatureHandles, and metadata populated.

    @param[in] table Arrow Table with consensus feature data
    @param[out] cmap ConsensusMap to populate
    @return true on success, false on error
  */
  static bool importFeaturesFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    ConsensusMap& cmap);

  /**
    @brief Import PSMs from Apache Arrow Table

    Reconstructs PeptideIdentifications and PeptideHits from the table
    and assigns them to the appropriate consensus features or as unassigned.

    @param[in] table Arrow Table with PSM data
    @param[out] cmap ConsensusMap to populate
    @return true on success, false on error
  */
  static bool importPSMsFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    ConsensusMap& cmap);

  /**
    @brief Import ConsensusMap from a directory of Parquet files

    Reads five Parquet files (consensus_features.parquet, psms.parquet,
    proteins.parquet, protein_groups.parquet, search_params.parquet)
    from the specified directory and reconstructs a complete ConsensusMap
    including FeatureHandles, PSM linkage, protein identifications,
    and ConsensusMap-level metadata.

    @param[in] directory Input directory path containing Parquet files
    @param[out] cmap ConsensusMap to populate
    @return true on success, false on error
  */
  static bool importFromParquet(
    const String& directory,
    ConsensusMap& cmap);
};

} // namespace OpenMS
