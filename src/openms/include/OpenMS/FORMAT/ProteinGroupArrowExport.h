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
  @brief Export protein group data to Apache Arrow format following QPX pg schema

  This class provides static methods to export protein group quantification data
  from a ConsensusMap to Apache Arrow Tables and Parquet files. The schema follows
  the QPX (Quantitative Proteomics Exchange) protein group format.

  Protein groups must have quantification annotated via
  PeptideAndProteinQuant::annotateQuantificationsToProteins() before export.

  @experimental This API is experimental and may change in future versions.

  @ingroup FileIO
*/
class OPENMS_DLLAPI ProteinGroupArrowExport
{
public:
  /**
    @brief Export protein group data to Apache Arrow Table

    Exports indistinguishable protein groups following the QPX pg schema.
    One row is emitted per protein group per run file.

    @param[in] cmap The ConsensusMap with annotated protein group quantification
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportToArrow(const ConsensusMap& cmap);

  /**
    @brief Export protein group data to Parquet file

    @param[in] cmap The ConsensusMap with annotated protein group quantification
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

#endif // WITH_PARQUET
