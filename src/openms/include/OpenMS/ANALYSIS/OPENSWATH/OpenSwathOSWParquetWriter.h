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
        └── run_id=&lt;id&gt;/
          ├── features.parquet
          ├── feature_precursor.parquet
          └── feature_transition.parquet
      @endcode

  If the output directory already exists, the writer will append a new
  run partition under runs/run_id=&lt;id&gt;/, update runs.parquet, and refresh
      metadata.json. Existing runs are not modified.
  */
  class OPENMS_DLLAPI OpenSwathOSWParquetWriter
  {
  public:
    /**
      @brief Write an OSW parquet directory or zip archive.

      Writes the assay library (precursors + transitions) and scored features
      for a single run into the Parquet directory layout. If the output
      directory already exists, the new run is appended under a fresh
      @c run_id= partition.

      @param[in] output_path       Path to the output directory (.oswpq)
      @param[in] assay_library     The transition library used for extraction
      @param[in] feature_map       Scored features for this run
      @param[in] run_id            Unique identifier for this run
      @param[in] input_filename    Original input filename (stored in metadata)
      @param[in] enable_uis_scoring Whether UIS (identification) transitions are included

      @throws Exception::MissingFeature if built without Parquet support
      @throws Exception::MissingInformation if a feature is missing required meta values
      @throws Exception::InvalidValue if a run_id already exists in the output
    */
    void write(const String& output_path,
               const OpenSwath::LightTargetedExperiment& assay_library,
               const FeatureMap& feature_map,
               UInt64 run_id,
               const String& input_filename,
               bool enable_uis_scoring) const;

    /**
     * @brief Set whether the writer should preserve (append) existing archives.
     * 
     * @param[in] preserve If true, the writer will append to existing .oswpq archives instead of overwriting them.
     */
    void setPreserveExisting(bool preserve);
  private:
    // Mutable so the const write() method may consult this flag.
    mutable bool preserve_existing_ = false;
  };

} // namespace OpenMS
