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
#include <OpenMS/FORMAT/QPXIdentity.h>

#include <memory>
#include <string>

// Forward declarations
namespace arrow
{
  class Table;
}

namespace OpenMS
{

/**
  @brief Export ConsensusMap feature data to Apache Arrow format following QPX feature schema

  This class provides static methods to export ConsensusMap data to Apache Arrow
  Tables and Parquet files. The schema follows the QPX (Quantitative Proteomics
  Exchange) feature format.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI ConsensusMapArrowExport
{
public:
  /**
    @brief Export ConsensusMap to Apache Arrow Table

    Exports consensus features following the QPX feature schema. Each
    ConsensusFeature becomes one row with identification, quantification,
    and protein group information.

    @param[in] cmap The ConsensusMap to export
    @param[out] out_links Optional feature&harr;PSM linkage, collected while the rows are built.
                Pass it to the psm exporter to fill @c psm.feature_id. Producing it here rather
                than in a second pass is what keeps the two directions reciprocal -- qpx
                validates that a PSM pointing at a feature is listed back by that feature's
                @c psm_ids.
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportToArrow(const ConsensusMap& cmap,
                                                     QPXIdentity::FeatureLinks* out_links = nullptr);

  /**
    @brief Verify that every feature can be attributed to a single peptide

    The QPX feature and pg views have singular @c sequence / @c peptidoform / @c charge and
    no ConsensusFeature identifier, so a feature whose identifications disagree on the top
    peptide cannot be represented: exporting one would publish a chosen interpretation while
    silently discarding the alternatives, and the pg view would count evidence contributed by
    one peptide under another's identity. The ambiguity-preserving representation is the
    OpenMS-native @c consensusparquet (ConsensusMapArrowIO), which keeps every hit and the
    consensus feature association, and which @ref TOPP_IDConflictResolver accepts as input.

    Several identifications agreeing on the same modified sequence
    (@c FEATURE_ID_MULTIPLE_SAME) are unambiguous and accepted -- that is the deliberate
    output of IDConflictResolverAlgorithm::resolve() with @p keep_matching.

    @note Must be called as a preflight, before any OpenMP region and before the output file
          is opened. An exception escaping the parallel batch build would be flattened into a
          boolean by the worker's catch-all, and one escaping the region at all is
          std::terminate.

    @param[in] cmap The map about to be exported
    @throws Exception::IllegalArgument if any feature has divergent peptide annotations
  */
  static void requireUnambiguousIdentities(const ConsensusMap& cmap);

  /**
    @brief Refuse a map whose features cannot be attributed to an origin MS run

    A feature's identification names its run through @c id_merge_index into the identification
    run's @c spectra_data. Without a usable index every PSM of a merged run resolves to the
    run's FIRST file, so @c run_file_name would be wrong rather than missing.

    Only the identification the exported row actually uses is validated -- validating every
    attached one refuses a feature whose winning hit resolves perfectly because a sibling,
    hitless or simply not selected, lacks an index that is never read.

    @note Same preflight constraint as requireUnambiguousIdentities(): call before any OpenMP
          region and before the output file is opened.

    @param[in] cmap The map about to be exported
    @throws Exception::MissingInformation if a feature's winning identification belongs to a
           merged run but carries no usable @c id_merge_index
  */
  static void requireResolvableIdRuns(const ConsensusMap& cmap);

  /**
    @brief Export ConsensusMap to Parquet file

    @param[in] cmap The ConsensusMap to export
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @param[out] out_links Optional feature&harr;PSM linkage, see exportToArrow()
    @return true on success, false on error
  */
  static bool exportToParquet(
    const ConsensusMap& cmap,
    const std::string& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{},
    QPXIdentity::FeatureLinks* out_links = nullptr);

  /**
    @brief Stream a ConsensusMap to a Parquet file in row batches (bounded peak memory)

    Functionally equivalent to exportToParquet() but builds and flushes the feature
    table one @p batch_size -sized range at a time through a persistent
    @c parquet::arrow::FileWriter, instead of materializing the whole ~N-row Arrow
    table in memory before a single write. For isobaric data (one consensus feature
    per PSM) N can be in the millions, where the one-shot path's transient peak drives
    the process into swap / OOM; here peak memory stays bounded by one batch.

    Each batch is optionally partitioned and built in parallel with OpenMP and written
    in index order (the Parquet writer stays serial), so the written rows and their order
    are identical to exportToParquet() and deterministic for any thread count; only the
    Parquet row-group layout may differ.

    @param[in] cmap The ConsensusMap to export
    @param[in] filename Output file path
    @param[in] batch_size Consensus features materialized per batch (0 is treated as the default)
    @param[in] config Parquet writing options
    @param[in] n_threads OpenMP threads for the per-batch build: 1 = serial (default),
                         0 = all available cores (honors @c OMP_NUM_THREADS), N = fixed
    @param[out] out_links Optional feature&harr;PSM linkage, see exportToArrow(). Each worker
                fills its own map and they are merged in index order, so the collected linkage
                does not depend on @p n_threads either.
    @return true on success, false on error
  */
  static bool exportToParquetStreaming(
    const ConsensusMap& cmap,
    const std::string& filename,
    size_t batch_size = 1000000,
    const ParquetWriteConfig& config = ParquetWriteConfig{},
    int n_threads = 1,
    QPXIdentity::FeatureLinks* out_links = nullptr);
};

} // namespace OpenMS
