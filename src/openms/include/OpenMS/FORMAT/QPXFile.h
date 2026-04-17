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
    const String& filename,
    bool export_all_psms = false,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Write a pre-built QPX PSM Arrow table to a Parquet file

    The table is expected to follow QPXPSMSchema (e.g., from exportPSMsToQPXArrow).
    Attaches QPX file metadata (qpx_version, file_type="psm", UUID, creation_date)
    before writing. Use this overload when the caller already has the table built
    (e.g., for merged output) to avoid rebuilding it.

    @param[in] table QPX PSM Arrow table (must not be null)
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
