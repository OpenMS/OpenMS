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
#include <OpenMS/METADATA/ID/MapPeptideIdentifications.h>
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

    Each PeptideHit from each PeptideIdentification of @p identifications (those of the features of @p feature_map,
    subordinates included, by their paths, and the unassigned ones) becomes one row, which names its feature by its
    unique id.

    @param[in] feature_map The FeatureMap whose features the identifications belong to
    @param[in] identifications The identifications to export, e.g. those of @p feature_map (see
               IdentificationDataConverter::exportMapIdentifications())
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportPSMsToArrow(
    const FeatureMap& feature_map,
    const MapPeptideIdentifications& identifications);

  /// The PSMs of the identifications of @p feature_map (see IdentificationDataConverter::exportMapIdentifications())
  static std::shared_ptr<arrow::Table> exportPSMsToArrow(
    const FeatureMap& feature_map);

  /**
    @brief Export FeatureMap to a directory of Parquet files

    Writes five Parquet files: features.parquet, psms.parquet,
    proteins.parquet, protein_groups.parquet, and search_params.parquet
    into the specified directory. FeatureMap-level metadata
    (DocumentIdentifier, DataProcessing) is stored as file-level
    key-value metadata in features.parquet. The PSM and protein tables stay empty:
    they are the layout of bundles of former versions, which hold peptide identifications.

    If the map owns identification data (FeatureMap::getIdentificationData()) or
    features carry native identification links, the bundle also contains the native
    identification bundle identifications/ (see IdentificationDataFile) and
    identification_links.parquet, which links features by unique ID. Linked features
    need distinct valid unique IDs, and every link must resolve.

    The bundle is written to a temporary sibling directory and published only when
    complete. An existing feature bundle or empty directory is replaced; any other
    existing file or directory is left untouched and the export fails.

    @throw Exception::InvalidValue for unresolved identification links

    @param[in] feature_map The FeatureMap to export
    @param[in] directory Output directory path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const FeatureMap& feature_map,
    const std::string& directory,
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
    and adds them to @p identifications: to the features of @p feature_map
    (by the paths of the features with the unique ids of the rows) or as unassigned.

    @param[in] table Arrow Table with PSM data
    @param[in] feature_map FeatureMap with the features the PSMs belong to
    @param[in,out] identifications The identifications to add to (their proteins give the score directions
                   of PSMs without one)
    @return true on success, false on error
  */
  static bool importPSMsFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    const FeatureMap& feature_map,
    MapPeptideIdentifications& identifications);

  /**
    @brief Import FeatureMap from a directory of Parquet files

    Reads five Parquet files (features.parquet, psms.parquet,
    proteins.parquet, protein_groups.parquet, search_params.parquet)
    from the specified directory and reconstructs a complete FeatureMap
    including feature hierarchy, identifications and FeatureMap-level metadata.
    The identifications are those of identifications/ (see exportToParquet()), or the
    peptide and protein identifications of a bundle of a former version, as identification data.

    @param[in] directory Input directory path containing Parquet files
    @param[out] feature_map FeatureMap to populate
    @return true on success, false on error
  */
  static bool importFromParquet(
    const std::string& directory,
    FeatureMap& feature_map);
};

} // namespace OpenMS
