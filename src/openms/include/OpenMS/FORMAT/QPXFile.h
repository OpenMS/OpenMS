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
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
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
  @brief Export PSM (Peptide Spectrum Match) data to Apache Arrow format following QPX PSM schema

  This class provides static methods to export PeptideIdentification/ProteinIdentification
  data to Apache Arrow Tables and Parquet files. The schema follows the QPX (Quantitative
  Proteomics Exchange) PSM format.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI QPXFile
{
public:
  /**
     * @brief Export PSMs to Arrow table using PSMSchema for lossless round-trips.
     *
     * Produces a table with PSMSchema columns (score, score_type, rank, etc.)
     * suitable for FeatureMapArrowIO and ConsensusMapArrowIO round-trips.
     * For QPX exchange format output, use exportPSMsToQPXArrow() instead.
     */
  static std::shared_ptr<arrow::Table> exportToArrow(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications,
    bool export_all_psms = false);

  /**
     * @brief Export PSMs to QPX Parquet eXchange format Arrow table (QPXPSMSchema).
     *
     * Unlike exportToArrow() which produces a PSMSchema table for lossless
     * round-trips, this method produces a QPXPSMSchema table optimized for
     * cross-tool exchange (quantms format).
     *
     * @param protein_identifications  Protein identifications (for file name lookup)
     * @param peptide_identifications  Peptide identifications to export
     * @param export_all_psms  If true, export all PSM hits; if false, only best hit per spectrum
     * @return Arrow table with QPXPSMSchema columns, or nullptr on failure
     */
  static std::shared_ptr<arrow::Table> exportPSMsToQPXArrow(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications,
    bool export_all_psms = false);

  /**
    @brief Export PSM data to Parquet file

    @param[in] protein_identifications Vector of protein identifications
    @param[in] peptide_identifications List of peptide identifications
    @param[in] filename Output file path
    @param[in] export_all_psms If true, export all hits per spectrum (default: false, only best hit)
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications,
    const std::string& filename,
    bool export_all_psms = false,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Write a pre-built QPX PSM Arrow table to a Parquet file

    The table is expected to follow QPXPSMSchema (e.g., from exportPSMsToQPXArrow).
    Attaches QPX file metadata (qpx_version, file_type="psm_file", UUID, creation_date,
    compression_format) before writing. Use this overload when the caller already has the
    table built (e.g., for merged output) to avoid rebuilding it.

    @param[in] table QPX PSM Arrow table (must not be null)
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @param[in] scan_format QPX scan_format token for the source native IDs, from
               ArrowIOHelpers::qpxScanFormat(). A pre-built table no longer carries the
               native IDs, so the caller must supply it; when empty the key is omitted
               rather than guessed.
    @return true on success, false on error
  */
  static bool exportToParquet(
    const std::shared_ptr<arrow::Table>& table,
    const std::string& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{},
    const std::string& scan_format = "");

  /**
    @brief Stream PSMs to a QPX Parquet file in row-batches to cap peak memory.

    Builds and flushes the QPXPSMSchema table in batches through a persistent Parquet
    writer instead of materializing the entire table in memory. Intended for very large
    inputs (e.g. millions of PSMs), where the one-shot exportToParquet would spike memory
    by holding all columns for all rows at once. Output is equivalent to exportToParquet
    (same schema and QPX metadata), but written as multiple Parquet row groups.

    @param[in] protein_identifications Protein identifications (for run-file-name lookup)
    @param[in] peptide_identification_ptrs NON-OWNING pointers to the PSMs to export. The
               caller guarantees they outlive the call and are non-null.
    @param[in] filename Output file path
    @param[in] export_all_psms If true, export all hits per spectrum (default: best hit only)
    @param[in] batch_size Number of PeptideIdentifications materialized into one in-memory
               Arrow table at a time. This is the peak-memory knob ONLY; it does not
               determine row-group count (with @p export_all_psms a single
               PeptideIdentification can emit several rows). 0 is treated as the default.
    @param[in] config Parquet writing options. config.row_group_size is the maximum number
               of rows per Parquet row group (the WriteTable chunk size).
    @param[in] n_threads OpenMP threads used to build each batch's partitions in parallel.
               1 = serial (default, preserves prior behaviour); 0 = auto (all available cores,
               i.e. omp_get_max_threads(), which honours the OMP_NUM_THREADS environment variable);
               N = fixed count. The per-row build dominates export cost, so parallelism here is
               the main speedup. Output is identical in row content and order regardless of
               @p n_threads (contiguous partitions are written in index order). The Parquet write
               itself stays serial. Without OpenMP support the export always runs serially.
               @note On hyper-threaded CPUs, using all logical cores (0) can be slower than the
               physical-core count because the build is memory-bandwidth bound; set OMP_NUM_THREADS
               (or pass N) to the physical-core count for best throughput on such machines.
    @return true on success, false on error (errors are logged)
  */
  static bool exportToParquetStreaming(
    const std::vector<ProteinIdentification>& protein_identifications,
    const std::vector<const PeptideIdentification*>& peptide_identification_ptrs,
    const std::string& filename,
    bool export_all_psms = false,
    size_t batch_size = 1000000,
    const ParquetWriteConfig& config = ParquetWriteConfig{},
    int n_threads = 1);

  /**
    @brief Refuse PSMs of a merged run that carry no usable @c id_merge_index

    Without the index every PSM of a merged identification run resolves to the run's FIRST
    file, so @c run_file_name -- a QPX primary-key component -- would be wrong rather than
    missing.

    @note Preflight: call before the output file is opened. The streaming build runs inside an
          OpenMP region whose exception firewall would flatten a throw into a logged
          @c return @c false, leaving a truncated .parquet behind.

    @param[in] protein_identifications Identification runs supplying @c spectra_data
    @param[in] peptide_identifications The PSMs about to be exported
    @throws Exception::MissingInformation if a PSM of a merged run has no usable index
    @throws Exception::InvalidValue if identification runs share an identifier
  */
  static void requireResolvableMergeIndices(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications);

  /**
    @brief Pointer-based overload of requireResolvableMergeIndices()

    Avoids copying every PSM when a caller already owns the identifications in another container.
    Null pointers are ignored, matching the streaming exporter.

    @param[in] protein_identifications Identification runs supplying @c spectra_data
    @param[in] peptide_identifications The PSMs about to be exported
    @throws Exception::MissingInformation if a PSM of a merged run has no usable index
    @throws Exception::InvalidValue if identification runs share an identifier
  */
  static void requireResolvableMergeIndices(
    const std::vector<ProteinIdentification>& protein_identifications,
    const std::vector<const PeptideIdentification*>& peptide_identifications);

  /**
    @brief Import PSMs from a PSMSchema Arrow table.

    Reads `PSMSchema`-conformant rows and appends `PeptideIdentification`s
    to @p peptide_identifications. Each row's `run_identifier` column links
    PSMs back to the matching `ProteinIdentification` already present in
    @p protein_identifications by run identifier. If no match exists, a
    new `ProteinIdentification` shell is appended.

    A shell also gets the run's ordered @c spectra_data back, rebuilt from the per-PSM
    @c reference_file_name and @c id_merge_index: on a merged run the index is what selects a PSM's
    file, so the whole list has to be restored - a single path would leave every index but 0
    pointing past the end. Indices are preserved, and a file that contributed no PSM stays an empty
    entry rather than being closed up, so IdentifierMSRunMapper resolves each imported PSM to the
    same origin file it had before the export. PSMs the table stores no usable index for are placed
    by their path instead, and their @c id_merge_index is set to match (merged runs only - on a
    single-file run index 0 is the default and none is written). Runs @p protein_identifications
    already covers are left alone: their own @c spectra_data is authoritative.

    @param[in]    table                     PSMSchema Arrow table (must not be null)
    @param[in,out] protein_identifications  Existing protein identifications (used for higher_score_better lookup; new shells appended for unknown run_identifiers)
    @param[in,out] peptide_identifications  Peptide identifications appended to (caller may pass an empty or pre-populated list)
    @return true on success, false on schema mismatch or unrecoverable error (errors are logged)
  */
  static bool importFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    std::vector<ProteinIdentification>& protein_identifications,
    PeptideIdentificationList& peptide_identifications);
};

} // namespace OpenMS
