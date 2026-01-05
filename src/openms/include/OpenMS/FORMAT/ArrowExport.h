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
#include <OpenMS/KERNEL/MSExperiment.h>

#include <memory>
#include <vector>
#include <string>

// Arrow C Data Interface types
#include <arrow/c/abi.h>

// Forward declarations - avoid exposing Arrow types in header
namespace arrow
{
  class Table;
}

namespace OpenMS
{

/**
  @brief Format for Arrow export

  Long format: One row per peak, with spectrum metadata repeated.
               Better for SQL-like queries and pandas operations.
               Memory overhead from repeated spectrum indices (~25-40% for 32-bit index).

  SemiWide format: One row per spectrum, with peak data as list arrays.
                   More memory efficient, handles empty spectra naturally.
                   Requires Arrow compute kernels for peak-level operations.
*/
enum class ArrowExportFormat
{
  Long,      ///< One row per peak (default)
  SemiWide   ///< One row per spectrum with list arrays for mz/intensity
};

/**
  @brief Configuration for Arrow export of spectra data

  Allows filtering by MS level, RT range, m/z range, and column selection.

  @note When ms_levels is empty, all MS levels are exported.
  @note When columns is empty, all available columns are exported.
  @note RT and m/z ranges of (0, 0) indicate no filtering.
*/
struct OPENMS_DLLAPI ArrowSpectraExportConfig
{
  /// Export format (Long or SemiWide)
  ArrowExportFormat format = ArrowExportFormat::Long;

  /// MS levels to include (empty = all levels)
  std::vector<UInt> ms_levels;

  /// Minimum RT in seconds (0 = no lower bound)
  double min_rt = 0;

  /// Maximum RT in seconds (0 = no upper bound)
  double max_rt = 0;

  /// Minimum m/z (0 = no lower bound)
  double min_mz = 0;

  /// Maximum m/z (0 = no upper bound)
  double max_mz = 0;

  /// Columns to export (empty = all available columns)
  /// Available columns depend on format:
  ///   Long: mz, intensity, rt, ion_mobility, spectrum_index, ms_level, native_id,
  ///         precursor_mz, precursor_charge, precursor_intensity,
  ///         isolation_lower, isolation_upper
  ///   SemiWide: Same but mz/intensity/ion_mobility are list arrays
  std::vector<std::string> columns;

  /// Include precursor information columns (default: true)
  bool include_precursor_info = true;

  /// Include ion mobility if present in data (default: true)
  bool include_ion_mobility = true;
};


/**
  @brief Export MSExperiment spectra to Apache Arrow Table

  Exports spectra data from an MSExperiment to an Arrow Table.

  Long format schema (one row per peak):
    - mz (float64): Peak m/z value
    - intensity (float32): Peak intensity
    - rt (float32): Retention time in seconds
    - ion_mobility (float32, nullable): Ion mobility value if present
    - spectrum_index (uint32): Index of spectrum in MSExperiment
    - ms_level (uint8): MS level (1, 2, ...)
    - native_id (string): Native spectrum identifier
    - precursor_mz (float64, nullable): Precursor m/z (null for MS1)
    - precursor_charge (int16, nullable): Precursor charge
    - precursor_intensity (float32, nullable): Precursor intensity
    - isolation_lower (float64, nullable): Isolation window lower offset
    - isolation_upper (float64, nullable): Isolation window upper offset

  Semi-wide format schema (one row per spectrum):
    - spectrum_index (uint32): Index of spectrum in MSExperiment
    - rt (float32): Retention time in seconds
    - ms_level (uint8): MS level
    - native_id (string): Native spectrum identifier
    - mz (list<float64>): Array of m/z values
    - intensity (list<float32>): Array of intensity values
    - ion_mobility (list<float32>, nullable): Array of ion mobility values
    - precursor_mz (float64, nullable): Precursor m/z (null for MS1)
    - precursor_charge (int16, nullable): Precursor charge
    - precursor_intensity (float32, nullable): Precursor intensity
    - isolation_lower (float64, nullable): Isolation window lower offset
    - isolation_upper (float64, nullable): Isolation window upper offset

  @param[in] exp The MSExperiment to export
  @param[in] config Export configuration (filtering, format, column selection)
  @return Shared pointer to Arrow Table, or nullptr on error

  @note Errors are logged via OPENMS_LOG_ERROR before returning nullptr.
  @note Empty experiments return a valid table with zero rows and the full schema.

  @ingroup FileIO
*/
OPENMS_DLLAPI std::shared_ptr<arrow::Table> exportSpectraToArrow(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config = ArrowSpectraExportConfig{});


/**
  @brief Get available column names for spectra Arrow export

  Returns the list of column names that would be included in the export
  based on the configuration and the actual data in the experiment.

  @param[in] exp The MSExperiment to analyze
  @param[in] config Export configuration
  @return Vector of available column names

  @ingroup FileIO
*/
OPENMS_DLLAPI std::vector<std::string> getSpectraArrowColumns(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config = ArrowSpectraExportConfig{});


/**
  @brief Configuration for Arrow export of chromatogram data
*/
struct OPENMS_DLLAPI ArrowChromatogramExportConfig
{
  /// Export format (Long or SemiWide)
  ArrowExportFormat format = ArrowExportFormat::Long;

  /// Minimum RT in seconds (0 = no lower bound)
  double min_rt = 0;

  /// Maximum RT in seconds (0 = no upper bound)
  double max_rt = 0;

  /// Columns to export (empty = all available columns)
  std::vector<std::string> columns;
};


/**
  @brief Export MSExperiment chromatograms to Apache Arrow Table

  Exports chromatogram data from an MSExperiment to an Arrow Table.

  Long format schema (one row per time point):
    - rt (float64): Retention time
    - intensity (float32): Signal intensity
    - chromatogram_index (uint32): Index of chromatogram in MSExperiment
    - native_id (string): Chromatogram identifier
    - precursor_mz (float64): Q1 m/z (precursor mass)
    - product_mz (float64): Q3 m/z (product mass for SRM)

  Semi-wide format schema (one row per chromatogram):
    - chromatogram_index (uint32): Index of chromatogram
    - native_id (string): Chromatogram identifier
    - rt (list<float64>): Array of retention times
    - intensity (list<float32>): Array of intensities
    - precursor_mz (float64): Q1 m/z
    - product_mz (float64): Q3 m/z

  @param[in] exp The MSExperiment to export
  @param[in] config Export configuration
  @return Shared pointer to Arrow Table, or nullptr on error

  @ingroup FileIO
*/
OPENMS_DLLAPI std::shared_ptr<arrow::Table> exportChromatogramsToArrow(
  const MSExperiment& exp,
  const ArrowChromatogramExportConfig& config = ArrowChromatogramExportConfig{});


/**
  @brief Get available column names for chromatogram Arrow export

  @param[in] exp The MSExperiment to analyze
  @param[in] config Export configuration
  @return Vector of available column names

  @ingroup FileIO
*/
OPENMS_DLLAPI std::vector<std::string> getChromatogramArrowColumns(
  const MSExperiment& exp,
  const ArrowChromatogramExportConfig& config = ArrowChromatogramExportConfig{});


/**
  @brief Export spectra to Arrow via C Data Interface (zero-copy to Python)

  Exports the Arrow schema and array to C Data Interface format, which allows
  zero-copy transfer to PyArrow via pyarrow.Table._import_from_c().

  @param[in] exp The MSExperiment to export
  @param[in] config Export configuration
  @param[out] out_schema Pointer to ArrowSchema struct (caller must allocate)
  @param[out] out_array Pointer to ArrowArray struct (caller must allocate)
  @return true on success, false on error

  @note The caller is responsible for calling the release callbacks on the
        schema and array when done.
  @note This is primarily intended for Python bindings for zero-copy export.

  @ingroup FileIO
*/
OPENMS_DLLAPI bool exportSpectraToArrowCDataInterface(
  const MSExperiment& exp,
  const ArrowSpectraExportConfig& config,
  ::ArrowSchema* out_schema,
  ::ArrowArray* out_array);


/**
  @brief Export chromatograms to Arrow via C Data Interface (zero-copy to Python)

  @param[in] exp The MSExperiment to export
  @param[in] config Export configuration
  @param[out] out_schema Pointer to ArrowSchema struct
  @param[out] out_array Pointer to ArrowArray struct
  @return true on success, false on error

  @ingroup FileIO
*/
OPENMS_DLLAPI bool exportChromatogramsToArrowCDataInterface(
  const MSExperiment& exp,
  const ArrowChromatogramExportConfig& config,
  ::ArrowSchema* out_schema,
  ::ArrowArray* out_array);

/**
  @brief Configuration for Parquet file writing

  Controls compression, row group size, and other Parquet-specific settings.
  These settings affect file size, read performance, and memory usage.

  Performance guidelines for MS data:
  - ZSTD compression typically achieves 3-5x compression on peak data
  - 128MB row groups balance parallelism and compression efficiency
  - Statistics enable efficient range queries on m/z and RT columns
*/
struct OPENMS_DLLAPI ParquetWriteConfig
{
  /// Compression algorithm
  enum class Compression
  {
    NONE,       ///< No compression (fastest write, largest files)
    SNAPPY,     ///< Fast compression, moderate ratio (widely compatible)
    GZIP,       ///< Good compression, slower (universal compatibility)
    LZ4,        ///< Very fast compression, lower ratio
    ZSTD        ///< Best ratio/speed tradeoff (recommended, default)
  };

  /// Compression algorithm (default: ZSTD for best ratio/speed)
  Compression compression = Compression::ZSTD;

  /// Compression level (interpretation depends on algorithm)
  /// ZSTD: 1-22 (default 3, higher = better ratio, slower)
  /// GZIP: 1-9 (default 6)
  /// LZ4/SNAPPY: ignored
  int compression_level = 3;

  /// Target row group size in bytes (default: 128MB)
  /// Smaller = more parallelism for readers, larger = better compression
  /// For MS data with millions of peaks, 128MB typically gives 1-2M rows per group
  int64_t row_group_size = 128 * 1024 * 1024;

  /// Write column statistics (min/max/null_count) for each row group
  /// Enables predicate pushdown for efficient m/z and RT range queries
  /// Small overhead (~1% file size increase)
  bool write_statistics = true;

  /// Data page size in bytes (default: 1MB)
  /// Affects granularity of reads within a row group
  int64_t data_page_size = 1024 * 1024;
};


/**
  @brief Export MSExperiment spectra to Parquet file

  Exports spectra data to Apache Parquet format, which provides:
  - Columnar storage optimized for analytical queries
  - Efficient compression (typically 3-5x for MS data with ZSTD)
  - Fast partial reads (only requested columns/row groups are loaded)
  - Wide ecosystem support (DuckDB, Polars, pandas, Spark, R arrow)

  The schema matches exportSpectraToArrow() output.

  Performance notes:
  - ZSTD compression gives best size/speed tradeoff for MS data
  - Row group size of 128MB balances parallelism and compression
  - Statistics enable efficient range queries on m/z and RT
  - For very large files (>100M peaks), consider exporting MS levels separately

  @param[in] exp The MSExperiment to export
  @param[in] filename Output file path (.parquet extension recommended)
  @param[in] config Arrow export configuration (filtering, format, columns)
  @param[in] parquet_config Parquet writing options (compression, row groups)
  @return true on success, false on error

  @note Errors are logged via OPENMS_LOG_ERROR before returning false.

  Example:
  @code
  MSExperiment exp;
  // ... load data ...

  // Export with default settings (ZSTD compression, 128MB row groups)
  exportSpectraToParquet(exp, "spectra.parquet");

  // Export MS2 only with maximum compression
  ArrowSpectraExportConfig config;
  config.ms_levels = {2};

  ParquetWriteConfig pq_config;
  pq_config.compression = ParquetWriteConfig::Compression::ZSTD;
  pq_config.compression_level = 9;

  exportSpectraToParquet(exp, "ms2_spectra.parquet", config, pq_config);
  @endcode

  @ingroup FileIO
*/
OPENMS_DLLAPI bool exportSpectraToParquet(
  const MSExperiment& exp,
  const String& filename,
  const ArrowSpectraExportConfig& config = ArrowSpectraExportConfig{},
  const ParquetWriteConfig& parquet_config = ParquetWriteConfig{});


/**
  @brief Export MSExperiment chromatograms to Parquet file

  Exports chromatogram data to Apache Parquet format.
  See exportSpectraToParquet() for details on Parquet benefits and options.

  @param[in] exp The MSExperiment to export
  @param[in] filename Output file path
  @param[in] config Arrow export configuration
  @param[in] parquet_config Parquet writing options
  @return true on success, false on error

  @ingroup FileIO
*/
OPENMS_DLLAPI bool exportChromatogramsToParquet(
  const MSExperiment& exp,
  const String& filename,
  const ArrowChromatogramExportConfig& config = ArrowChromatogramExportConfig{},
  const ParquetWriteConfig& parquet_config = ParquetWriteConfig{});


} // namespace OpenMS

#endif // WITH_PARQUET
