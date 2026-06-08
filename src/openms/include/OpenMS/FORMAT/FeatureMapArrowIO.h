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
#include <OpenMS/KERNEL/FeatureMap.h>
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
  @brief Import and export FeatureMap data to/from Apache Arrow format

  This class provides static methods to export and import FeatureMap
  data to/from Apache Arrow Tables and Parquet files. Separate tables are
  provided for features (with their geometry and metadata) and for
  peptide spectrum matches (PSMs) associated with features.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI FeatureMapArrowIO
{
public:
  // ==================== Export methods ====================

  /**
    @brief Export features to Apache Arrow Table

    Each Feature becomes one row with RT, MZ, intensity, charge,
    quality, convex hull geometry, and metadata columns.

    @param[in] feature_map The FeatureMap to export
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportFeaturesToArrow(
    const FeatureMap& feature_map);

  /**
    @brief Export peptide spectrum matches (PSMs) associated with features to Apache Arrow Table

    Each PeptideHit from each PeptideIdentification (both feature-level
    and unassigned) becomes one row.

    @param[in] feature_map The FeatureMap whose identifications to export
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportPSMsToArrow(
    const FeatureMap& feature_map);

  /**
    @brief Export FeatureMap to a directory of Parquet files

    Writes five Parquet files: features.parquet, psms.parquet,
    proteins.parquet, protein_groups.parquet, and search_params.parquet
    into the specified directory. Protein-level data is delegated to
    ProteinIdentificationArrowIO. FeatureMap-level metadata
    (DocumentIdentifier, DataProcessing) is stored as file-level
    key-value metadata in features.parquet.

    @param[in] feature_map The FeatureMap to export
    @param[in] directory Output directory path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const FeatureMap& feature_map,
    const String& directory,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  // ==================== Import methods ====================

  /**
    @brief Import features from Apache Arrow Table

    Each row becomes a Feature with RT, MZ, intensity, charge,
    quality, convex hull geometry, and metadata populated.

    @param[in] table Arrow Table with feature data
    @param[out] feature_map FeatureMap to populate
    @return true on success, false on error
  */
  static bool importFeaturesFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    FeatureMap& feature_map);

  /**
    @brief Import PSMs from Apache Arrow Table

    Reconstructs PeptideIdentifications and PeptideHits from the table
    and assigns them to the appropriate features or as unassigned.

    @param[in] table Arrow Table with PSM data
    @param[out] feature_map FeatureMap to populate
    @return true on success, false on error
  */
  static bool importPSMsFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    FeatureMap& feature_map);

  /**
    @brief Import FeatureMap from a directory of Parquet files

    Reads five Parquet files (features.parquet, psms.parquet,
    proteins.parquet, protein_groups.parquet, search_params.parquet)
    from the specified directory and reconstructs a complete FeatureMap
    including feature hierarchy, PSM linkage, protein identifications,
    and FeatureMap-level metadata.

    @param[in] directory Input directory path containing Parquet files
    @param[out] feature_map FeatureMap to populate
    @return true on success, false on error
  */
  static bool importFromParquet(
    const String& directory,
    FeatureMap& feature_map);
};

} // namespace OpenMS
