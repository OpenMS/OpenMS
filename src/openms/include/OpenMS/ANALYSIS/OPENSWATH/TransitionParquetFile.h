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

namespace OpenSwath
{
  struct LightTargetedExperiment;
}

namespace OpenMS
{

  /**
      @brief Read OpenSwath Parquet library input (.oswpq) into LightTargetedExperiment.

      The Parquet library format is a directory container with separate tables for
      precursors and transitions. The reader materializes all rows into
      OpenSwath::LightTargetedExperiment.

      The container layout is:
      @code
      <library>.oswpq
      └── library/
          ├── metadata.json
          ├── precursors.parquet
          └── transitions.parquet
      @endcode

      The metadata file contains an OpenMS metadata block and QC counts. The
      new canonical layout (matching OpenSwathOSWParquetWriter) looks like:
      @code{.json}
      {
        "openms": {
          "schema_version": 1,
          "generator": "TransitionParquetFile",
          "openms_version": "<version>",
          "build_time": "<build_time>",
          "tool": {"name": "OpenSwathWorkflow", "version": "<version>"},
          "counts": {
            "proteins": {"total": 0, "target": 0, "decoy": 0},
            "peptides": {"total": 0, "target": 0, "decoy": 0},
            "precursors": {"total": 0, "target": 0, "decoy": 0},
            "compounds": {"total": 0, "target": 0, "decoy": 0},
            "transitions": {"total": 0, "target": 0, "decoy": 0}
          },
          "fragment_type_counts": {
            "target": {"b": 0, "y": 0, "other": 0},
            "decoy": {"b": 0, "y": 0, "other": 0}
          },
          "charge_counts": {
            "precursor": {"target": {"2": 0, "3": 0}, "decoy": {"2": 0, "3": 0}},
            "transition": {"target": {"1": 0, "2": 0}, "decoy": {"1": 0, "2": 0}}
          }
        }
      }
      @endcode

      Required columns for @c precursors.parquet:
      - precursor_id (int64)
      - precursor_mz (float64)
      - charge (int32)
      - library_rt (float64)
      Optional columns:
      - library_drift_time (float64)
      - traml_id (string)
      - decoy (bool)
      - modified_sequence (string)
      - unmodified_sequence (string)
      - protein_accessions (string or list<string>)

      Required columns for @c transitions.parquet:
      - transition_id (int64)
      - precursor_id (int64)
      - product_mz (float64)
      - charge (int32)
      - type (string)
      - ordinal (int32)
      - detecting (bool)
      - identifying (bool)
      - quantifying (bool)
      - library_intensity (float64)
      - decoy (bool)
      Optional columns:
      - traml_id (string)
      - annotation (string)
  */
  class OPENMS_DLLAPI TransitionParquetFile
  {
  public:
    /// Default constructor
    TransitionParquetFile() = default;

    /// Default destructor
    ~TransitionParquetFile() = default;

    /// Read a .oswpq library directory and populate a LightTargetedExperiment
    void convertParquetToTargetedExperiment(const String& oswpq_dir, OpenSwath::LightTargetedExperiment& targeted_exp) const;

    /// Write a LightTargetedExperiment to a .oswpq library (zip file or directory)
    void convertLightTargetedExperimentToParquet(const String& oswpq_path, const OpenSwath::LightTargetedExperiment& targeted_exp) const;
  };

} // namespace OpenMS
