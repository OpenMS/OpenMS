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

    Each PeptideHit from each PeptideIdentification of @p identifications (those of the consensus features of
    @p cmap by their indices, and the unassigned ones) becomes one row, which names its consensus feature by its
    unique id.

    @param[in] cmap The ConsensusMap whose consensus features the identifications belong to
    @param[in] identifications The identifications to export, e.g. those of @p cmap (see
               IdentificationDataConverter::exportMapIdentifications())
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportPSMsToArrow(
    const ConsensusMap& cmap,
    const MapPeptideIdentifications& identifications);

  /// The PSMs of the identifications of @p cmap (see IdentificationDataConverter::exportMapIdentifications())
  static std::shared_ptr<arrow::Table> exportPSMsToArrow(
    const ConsensusMap& cmap);

  /**
    @brief Export ConsensusMap to a directory of Parquet files

    Writes five Parquet files: consensus_features.parquet, psms.parquet,
    proteins.parquet, protein_groups.parquet, and search_params.parquet
    into the specified directory. ConsensusMap-level metadata
    (column headers, experiment type, DocumentIdentifier, DataProcessing)
    is stored as file-level key-value metadata in consensus_features.parquet.
    The PSM and protein tables stay empty: they are the layout of bundles of
    former versions, which hold peptide identifications.

    If the map owns identification data (ConsensusMap::getIdentificationData()) or
    consensus features carry native identification links, the bundle also contains the
    native identification bundle identifications/ (see IdentificationDataFile) and
    identification_links.parquet, which links consensus features by unique ID.

    The bundle is written to a temporary sibling directory and published only when
    complete. An existing consensus bundle or empty directory is replaced; any other
    existing file or directory is left untouched and the export fails.

    @throw Exception::InvalidValue for unresolved identification links

    @param[in] cmap The ConsensusMap to export
    @param[in] directory Output directory path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const ConsensusMap& cmap,
    const std::string& directory,
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
    and adds them to @p identifications: to the consensus features of @p cmap
    (by the indices of the consensus features with the unique ids of the rows) or as unassigned.
    A run identifier that names none of the protein identifications of @p identifications
    adds a protein identification run without hits (see QPXFile::importFromArrow()).

    @param[in] table Arrow Table with PSM data
    @param[in] cmap ConsensusMap with the consensus features the PSMs belong to
    @param[in,out] identifications The identifications to add to
    @return true on success, false on error
  */
  static bool importPSMsFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    const ConsensusMap& cmap,
    MapPeptideIdentifications& identifications);

  /**
    @brief Import ConsensusMap from a directory of Parquet files

    Reads five Parquet files (consensus_features.parquet, psms.parquet,
    proteins.parquet, protein_groups.parquet, search_params.parquet)
    from the specified directory and reconstructs a complete ConsensusMap
    including FeatureHandles, identifications and ConsensusMap-level metadata.
    The identifications are those of identifications/ (see exportToParquet()), or the
    peptide and protein identifications of a bundle of a former version, as identification data.

    @param[in] directory Input directory path containing Parquet files
    @param[out] cmap ConsensusMap to populate
    @return true on success, false on error
  */
  static bool importFromParquet(
    const std::string& directory,
    ConsensusMap& cmap);
};

} // namespace OpenMS
