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

#include <memory>
#include <vector>

namespace arrow
{
  class Table;
}

namespace OpenMS
{

/**
  @brief PSM Arrow I/O for the idparquet bundle format.

  Provides export/import of PeptideIdentification data to/from Arrow Tables
  using PSMSchema. Optimized for fast internal round-trips by storing the
  native OpenMS modified sequence string (AASequence::toString()) in the
  sequence column, enabling direct AASequence::fromString() on import
  without ProForma conversion overhead.

  @ingroup FileIO
*/
class OPENMS_DLLAPI PSMParquetIO
{
public:
  /**
   * @brief Export PSMs to Arrow table using PSMSchema.
   *
   * Stores the full modified sequence (AASequence::toString()) in the
   * sequence column for fast import. Also writes the ProForma canonical
   * form in the peptidoform column for external tool interoperability.
   */
  static std::shared_ptr<arrow::Table> exportToArrow(
    const std::vector<ProteinIdentification>& protein_identifications,
    const PeptideIdentificationList& peptide_identifications,
    bool export_all_psms = false);

  /**
   * @brief Import PSMs from Arrow table using PSMSchema.
   *
   * Reads the sequence column directly via AASequence::fromString(),
   * bypassing the slower ProForma parse path.
   */
  static bool importFromArrow(
    const std::shared_ptr<arrow::Table>& table,
    std::vector<ProteinIdentification>& protein_identifications,
    PeptideIdentificationList& peptide_identifications);
};

} // namespace OpenMS
