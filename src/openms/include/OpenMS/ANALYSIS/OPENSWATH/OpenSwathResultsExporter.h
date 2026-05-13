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
    @brief Write filtered OpenSWATH feature results to TSV or Parquet.

    The exporter is file-format agnostic with respect to the input rows. OSW- or
    Parquet-specific reading remains outside this class.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathResultsExporter
  {
  public:
    /// Write @p rows to @p filename using the selected output format.
    /// Empty exports are allowed and result in a header-only TSV or an empty Parquet table.
    static void write(const String& filename,
                      const std::vector<OpenSwathExportRow>& rows,
                      const OpenSwathResultsExportConfig& config);
  };
} // namespace OpenMS
