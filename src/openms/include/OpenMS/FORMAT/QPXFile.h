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
    @brief Export PSM data to Apache Arrow Table

    Exports peptide spectrum matches following the QPX PSM schema. Each
    PeptideHit becomes one row with identification, score, and metadata.

    @param[in] protein_identifications Vector of protein identifications
    @param[in] peptide_identifications List of peptide identifications
    @param[in] export_all_psms If true, export all hits per spectrum (default: false, only best hit)
    @return Shared pointer to Arrow Table, or nullptr on error
  */
  static std::shared_ptr<arrow::Table> exportToArrow(
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
};

} // namespace OpenMS

#endif // WITH_PARQUET
