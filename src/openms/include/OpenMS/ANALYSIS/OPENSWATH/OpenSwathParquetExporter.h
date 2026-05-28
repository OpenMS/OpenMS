// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportData.h>
#include <OpenMS/config.h>

namespace OpenMS
{
  /**
    @brief Write scored OpenSWATH feature and transition tables to Parquet.

    @ingroup TargetedQuantitation
  */
  class OPENMS_DLLAPI OpenSwathParquetExporter
  {
  public:
    /// Write precursor / feature score rows to Parquet.
    /// Empty tables are allowed as long as the schema metadata is available.
    static void writeFeatureScores(const String& filename,
                                   const OpenSwathFeatureScoreTable& table);

    /// Write optional transition score rows to Parquet.
    /// Empty tables are allowed as long as the schema metadata is available.
    static void writeTransitionScores(const String& filename,
                                      const OpenSwathTransitionScoreTable& table);
  };
} // namespace OpenMS
