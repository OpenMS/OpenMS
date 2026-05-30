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
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/MSExperimentArrowExport.h>

#include <vector>

namespace OpenMS
{

/**
  @brief Read and write OpenMS identification data as a parquet bundle (.idparquet).

  An idparquet bundle is a directory containing four parquet files:
    - psms.parquet           (PSMSchema, lossless PeptideIdentification data)
    - proteins.parquet       (ProteinSchema, ProteinHits)
    - protein_groups.parquet (ProteinGroupSchema, indistinguishable groups)
    - search_params.parquet  (SearchParamsSchema, run-level parameters)

  All four files are required for a valid bundle on read.

  @ingroup FileIO
*/
class OPENMS_DLLAPI PSMArrowIO
{
public:
  /**
    @brief Export protein and peptide identifications to an idparquet directory bundle.

    Writes (and overwrites) the four canonical files inside @p dir. Other
    files in @p dir are left untouched. If @p dir does not exist it is
    created (only @p dir itself; the parent must already exist). If @p dir
    exists as a regular file, returns false.

    @param[in] protein_identifications  Protein identifications (run-level metadata, hits, groups, search params)
    @param[in] peptide_identifications  Peptide identifications (PSMs)
    @param[in] dir                       Output directory path
    @param[in] export_all_psms           If true, write all hits per PSM; if false, write best hit per spectrum only (default: true — bundle is intended for lossless round-trip)
    @param[in] config                    Parquet writer configuration
    @return true on success, false on error (errors are logged)
  */
  static bool exportToParquet(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications,
    const String& dir,
    bool export_all_psms = true,
    const ParquetWriteConfig& config = ParquetWriteConfig{});

  /**
    @brief Import protein and peptide identifications from an idparquet directory bundle.

    All four canonical files (psms.parquet, proteins.parquet,
    protein_groups.parquet, search_params.parquet) must be present in
    @p dir. Missing any one is an error.

    @param[in]  dir                       Input directory path
    @param[out] protein_identifications   Populated from the three protein-side parquet files
    @param[out] peptide_identifications   Populated from psms.parquet
    @return true on success, false on error (errors are logged)
  */
  static bool importFromParquet(
    const String& dir,
    std::vector<ProteinIdentification>& protein_identifications,
    PeptideIdentificationList& peptide_identifications);
};

} // namespace OpenMS
