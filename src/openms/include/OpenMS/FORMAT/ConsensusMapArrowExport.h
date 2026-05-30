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
    const String& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});
};

} // namespace OpenMS
