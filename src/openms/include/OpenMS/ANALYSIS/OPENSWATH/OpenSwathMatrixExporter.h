// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportConfig.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportData.h>
#include <OpenMS/config.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Build and write OpenSWATH quantification matrices.

    Matrix construction follows the same high-level PyProphet logic: select the
    best peakgroup per run / precursor first, then summarize to precursor,
    peptide, protein, or gene level.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathMatrixExporter
  {
  public:
    /// Build a matrix from filtered export rows.
    static OpenSwathQuantMatrix buildMatrix(const std::vector<OpenSwathExportRow>& rows,
                                            const OpenSwathMatrixExportConfig& config);

    /// Write a previously built matrix to TSV or Parquet.
    static void writeMatrix(const String& filename,
                            const OpenSwathQuantMatrix& matrix,
                            const OpenSwathMatrixExportConfig& config);
  };
} // namespace OpenMS
