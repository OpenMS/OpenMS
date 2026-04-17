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
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <memory>
#include <string>
#include <vector>

// Forward declarations
namespace arrow
{
  class Table;
}

namespace OpenMS
{
  class ProteinIdentification;
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

  /**
    @brief Export protein group data to Arrow table from identification data (no quantification)

    For search-engine output where no ConsensusMap is available. Populates required
    QPX pg fields (pg_accessions, anchor_protein, run_file_name, is_decoy, peptides)
    and sets quantification columns (intensities, additional_intensities) to null.

    @param[in] protein_identifications Protein identifications with protein groups
    @param[in] peptide_identifications Peptide identifications (for peptide-per-protein counts)
    @return Shared pointer to Arrow Table (empty table if no groups, never nullptr)
  */
  static std::shared_ptr<arrow::Table> exportToArrow(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications);

  /**
    @brief Export protein group data to Parquet file from identification data (no quantification)

    @param[in] protein_identifications Protein identifications with protein groups
    @param[in] peptide_identifications Peptide identifications (for peptide-per-protein counts)
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications,
    const String& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Write a pre-built QPX protein group Arrow table to a Parquet file

    The table is expected to follow QPXPgSchema (e.g., from exportToArrow).
    Attaches QPX file metadata (qpx_version, file_type="pg", UUID, creation_date)
    before writing. Use this overload when the caller already has the table built
    (e.g., for merged output) to avoid rebuilding it.

    @param[in] table QPX pg Arrow table (must not be null)
    @param[in] filename Output file path
    @param[in] config Parquet writing options
    @return true on success, false on error
  */
  static bool exportToParquet(
    const std::shared_ptr<arrow::Table>& table,
    const String& filename,
    const ParquetWriteConfig& config = ParquetWriteConfig{});
};

} // namespace OpenMS
