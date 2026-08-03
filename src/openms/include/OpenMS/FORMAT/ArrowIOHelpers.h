// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>  // for ParquetWriteConfig

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

// Forward declarations
namespace arrow
{
  class Array;
  class KeyValueMetadata;
  class Table;
}

namespace OpenMS
{
class ConsensusMap;
class MetaInfoInterface;

/**
  @brief Public helpers for writing and concatenating Arrow tables to Parquet files

  TOPP tools link against libOpenMS (which exports these helpers) but not directly
  against Arrow/Parquet. These wrappers keep all Arrow/Parquet API calls inside
  libOpenMS so downstream binaries don't need to import Arrow symbols.

  @ingroup FileIO
*/
namespace ArrowIOHelpers
{
  /**
    @brief Generate a lowercase hyphenated RFC 4122 version-4 UUID string

    Used by QPX Parquet exporters when attaching file metadata.

    @return UUID string, e.g. "550e8400-e29b-41d4-a716-446655440000"
  */
  OPENMS_DLLAPI std::string generateUuidV4();

  /**
    @brief Write an Arrow table to a Parquet file

    @param[in] table The Arrow table to write (must not be null)
    @param[in] filename Output file path
    @param[in] config Parquet writer configuration (compression, row group size, ...)
    @return true on success, false on error (errors are logged)
  */
  OPENMS_DLLAPI bool writeTableToParquet(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Concatenate a vector of Arrow tables and write the result to a Parquet file

    All tables must share the same schema. An empty input vector is a no-op
    (returns true without writing).

    @param[in] tables Vector of Arrow tables to concatenate (must share schema)
    @param[in] filename Output file path
    @param[in] config Parquet writer configuration
    @param[in] metadata Schema key-value metadata for the merged file. Concatenation
               otherwise inherits the first input's metadata, i.e. that table's identity;
               pass a fresh one (see qpxFileMetadata()) so the merged file is its own.
               @c nullptr keeps whatever the concatenation produced.
    @return true on success (or if @p tables is empty), false on error
  */
  OPENMS_DLLAPI bool concatenateAndWriteToParquet(
    const std::vector<std::shared_ptr<arrow::Table>>& tables,
    const std::string& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{},
    const std::shared_ptr<const arrow::KeyValueMetadata>& metadata = nullptr);

  // ---------------------------------------------------------------------------
  // QPX file metadata
  // ---------------------------------------------------------------------------

  /**
    @brief Build the canonical QPX file-level key-value metadata

    Writes the keys defined by the QPX serialization spec: @c qpx_version,
    @c file_type, @c creator, @c software_provider, @c creation_date (ISO 8601),
    @c compression_format and @c uuid, plus any @p extra keys.

    Build this <b>once per output file</b> and reuse the returned object for both the
    writer schema and every batch — each call mints a fresh @c uuid and
    @c creation_date, so calling it per batch would produce mismatched schemas.

    @param[in] file_type QPX view token: @c "psm_file", @c "feature_file" or @c "pg_file"
    @param[in] config Write configuration; supplies @c compression_format
    @param[in] extra Additional keys, e.g. <tt>{{"scan_format", "scan"}}</tt>
    @return The metadata, or @c nullptr if @p config selects a compression QPX does not
            define (LZ4). Callers must treat @c nullptr as a write failure.
  */
  OPENMS_DLLAPI std::shared_ptr<const arrow::KeyValueMetadata> qpxFileMetadata(
    const std::string& file_type,
    const ParquetWriteConfig& config = ParquetWriteConfig{},
    const std::map<std::string, std::string>& extra = {});

  /**
    @brief Classify a spectrum native ID into a QPX @c scan_format token

    @param[in] native_id A spectrum native ID (e.g. <tt>controllerType=0 ... scan=1234</tt>)
    @return @c "index" for @c index= IDs, @c "scan" for other recognized native IDs,
            or @c "" when the convention cannot be determined.

    @see qpxScanFormat(const std::vector<std::string>&)
  */
  OPENMS_DLLAPI std::string qpxScanFormat(const std::string& native_id);

  /**
    @brief Derive a single QPX @c scan_format token for a set of native IDs

    Unrecognized IDs are ignored. Returns @c "" when no ID is recognized, or when the
    inputs disagree — mixed conventions are reported once via the log rather than
    guessed at, so an ambiguous export omits @c scan_format instead of mislabeling it.
  */
  OPENMS_DLLAPI std::string qpxScanFormat(const std::vector<std::string>& native_ids);

  /**
    @brief Reduce an MS run path to the QPX @c run_file_name form

    QPX defines @c run_file_name as the raw data file name *without* extension, and uses it
    as a primary-key component in the psm, feature and pg views (and to join them against
    @c run.parquet). The full name with extension belongs in @c run.file_name.

    @param[in] ms_run_path Source path, e.g. <tt>/data/proj/S1_Frontal_1.mzML</tt>
    @return The stem, e.g. <tt>S1_Frontal_1</tt>; empty input yields an empty result.
  */
  OPENMS_DLLAPI std::string qpxRunFileName(const std::string& ms_run_path);

  /**
    @brief Warn once per export about source files that collapse onto one @c run_file_name

    QPX defines the column as the file name without path or extension, so two distinct paths
    sharing a stem (@c /a/run1.mzML and @c /b/run1.mzML) become indistinguishable as a join
    and partition key. Same-named files in different directories are a legitimate layout and
    the origin is known -- only the spec's representation cannot express it -- so this warns
    rather than failing, matching the policy the PSM exporter already applies.

    @param[in] context Caller name for the log line
    @param[in] ms_run_paths Source paths participating in one QPX collection
    @return true if every distinct path yields a distinct run name
  */
  OPENMS_DLLAPI bool qpxWarnOnRunNameCollisions(
    const std::string& context,
    const std::vector<std::string>& ms_run_paths);


  /**
    @brief QPX @c intensities[].label for a ConsensusMap column header

    The label is a join key: QPX matches @c intensities[].label against the unnested
    @c run.samples[].label of @c run.parquet (docs/spec/views.md), so it must be a
    canonical channel token, not a channel index or a file name.

    Isobaric channels resolve to the reporter name prefixed by the method family, using OpenMS'
    own channel names — @c "TMT126", @c "TMT131" (TMT10-plex channel 10), @c "ITRAQ114".
    Label-free maps resolve to @c "LFQ". Recognizable SILAC modification labels resolve by mass
    class; callers that have a whole ConsensusMap should prefer qpxIntensityLabels(), which can
    distinguish a two-plex second channel (heavy by role) from a three-plex medium channel.

    @param[in] column_label ConsensusMap::ColumnHeader::label. IsobaricChannelExtractor
               writes @c "&lt;methodname&gt;_&lt;channelname&gt;" (e.g. @c "tmt10plex_126");
               ProteomicsLFQ writes @c "label-free".
    @param[in] channel_name The header's @c channel_name meta value (e.g. @c "126"); empty
               when the map is not isobaric.
    @return The label, or @c "" when the isobaric method, non-isobaric vocabulary token, or
            SILAC role cannot be identified — writing a guessed token into a join key is worse
            than writing none, so callers must handle the empty result.
  */
  OPENMS_DLLAPI std::string qpxIntensityLabel(
    const std::string& column_label,
    const std::string& channel_name);

  /**
    @brief Derive canonical QPX intensity labels for every ConsensusMap column

    Unlike the scalar overload, this sees all columns of one source run. It recognizes the
    FeatureFinderMultiplex SILAC shapes and maps them to the active SDRF/QPX vocabulary:
    @c "SILAC light" / @c "SILAC heavy" for two-plex and
    @c "SILAC light" / @c "SILAC medium" / @c "SILAC heavy" for three-plex. The mapping is by
    channel role, so a two-plex @c Arg6 channel is correctly called heavy even though the same
    modification is the medium channel in the standard three-plex.

    @param[in] cmap ConsensusMap whose column headers are labelled
    @return Map index to canonical label. An unrepresentable header maps to @c "" and is logged;
            exporters must refuse it.
  */
  OPENMS_DLLAPI std::map<std::uint64_t, std::string> qpxIntensityLabels(
    const ConsensusMap& cmap);

  /**
    @brief Test whether a string is in the canonical SDRF/QPX intensity-label vocabulary

    Covers every isobaric method supported by OpenMS plus the QPX LFQ, SILAC, mTRAQ, and dimethyl
    plex labels. Intended for write-time value validation, after producer-specific labels have
    been normalized.

    @param[in] label Candidate label
    @return true if @p label is a canonical QPX channel token
  */
  OPENMS_DLLAPI bool qpxIsCanonicalIntensityLabel(const std::string& label);

  // ---------------------------------------------------------------------------
  // Read helpers
  // ---------------------------------------------------------------------------

  /**
    @brief Fetch a named column from a table, combining chunks if needed

    Returns nullptr if the column is missing or contains no chunks. When
    @p required is true, missing columns are logged as errors.
  */
  OPENMS_DLLAPI std::shared_ptr<arrow::Array> getColumn(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& name,
    bool required = true);

  /// Read a string at @p row, or "" if null/out-of-bounds
  OPENMS_DLLAPI std::string getStringValue(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row);

  /// Read a double at @p row, or @p default_val if null
  OPENMS_DLLAPI double getDoubleValue(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    double default_val = 0.0);

  /// Read a float at @p row, or @p default_val if null
  OPENMS_DLLAPI float getFloatValue(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    float default_val = 0.0f);

  /// Read an int32 at @p row, or @p default_val if null
  OPENMS_DLLAPI int32_t getInt32Value(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    int32_t default_val = 0);

  /// Read an int64 at @p row, or @p default_val if null
  OPENMS_DLLAPI int64_t getInt64Value(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    int64_t default_val = 0);

  /// Read a bool at @p row, or @p default_val if null
  OPENMS_DLLAPI bool getBoolValue(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    bool default_val = false);

  /// Whether @p array is null at @p row (or unset)
  OPENMS_DLLAPI bool isNull(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row);

  /**
    @brief Read metavalues from a list<struct{name,value,value_type}> column

    Decodes typed entries (int, double/float, *_list, string) and assigns
    them to @p target. Keys in @p excluded_keys are skipped.
  */
  OPENMS_DLLAPI void readMetaValues(
    const std::shared_ptr<arrow::Array>& array,
    int64_t row,
    MetaInfoInterface& target,
    const std::unordered_set<std::string>& excluded_keys = {});
}

} // namespace OpenMS
