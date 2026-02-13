// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenSwath
{
  struct LightTargetedExperiment;
}

namespace OpenMS
{

  /**
      @brief Write OpenSwathWorkflow output to an OSW Parquet directory (.oswpq).

      This writer mirrors the OSW SQLite tables in Parquet:
      - runs.parquet
      - features.parquet
      - feature_precursor.parquet
      - feature_transition.parquet (including UIS rows when present)

      The output is a directory or zip archive with the following layout:
      @code
      <output>.oswpq
      ├── metadata.json
      ├── library/
      │   ├── precursors.parquet
      │   └── transitions.parquet
      └── runs/
          ├── runs.parquet
          └── run_id=<id>/
              ├── features.parquet
              ├── feature_precursor.parquet
              └── feature_transition.parquet
      @endcode

      If the output directory already exists, the writer will append a new
      run partition under runs/run_id=<id>/, update runs.parquet, and refresh
      metadata.json. Existing runs are not modified.
  */
  class OPENMS_DLLAPI OpenSwathOSWParquetWriter
  {
  public:
    /// Write an OSW parquet directory or zip archive
    void write(const String& output_path,
               const OpenSwath::LightTargetedExperiment& assay_library,
               const FeatureMap& feature_map,
               UInt64 run_id,
               const String& input_filename,
               bool enable_uis_scoring) const;
  };

} // namespace OpenMS
