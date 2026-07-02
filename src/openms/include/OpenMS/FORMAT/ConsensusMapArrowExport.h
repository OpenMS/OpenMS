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
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportToArrow(const ConsensusMap& cmap);

  /**
    @brief Export ConsensusMap to Parquet file

    @param[in] cmap The ConsensusMap to export
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const ConsensusMap& cmap,
    const std::string& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Stream a ConsensusMap to a Parquet file in row batches (bounded peak memory)

    Functionally equivalent to exportToParquet() but builds and flushes the feature
    table one @p batch_size -sized range at a time through a persistent
    @c parquet::arrow::FileWriter, instead of materializing the whole ~N-row Arrow
    table in memory before a single write. For isobaric data (one consensus feature
    per PSM) N can be in the millions, where the one-shot path's transient peak drives
    the process into swap / OOM; here peak memory stays bounded by one batch.

    The written rows and their order are identical to exportToParquet(); only the
    Parquet row-group layout may differ.

    @param[in] cmap The ConsensusMap to export
    @param[in] filename Output file path
    @param[in] batch_size Consensus features materialized per batch (0 is treated as the default)
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquetStreaming(
    const ConsensusMap& cmap,
    const std::string& filename,
    size_t batch_size = 1000000,
    const ParquetWriteConfig& config = ParquetWriteConfig{});
};

} // namespace OpenMS
