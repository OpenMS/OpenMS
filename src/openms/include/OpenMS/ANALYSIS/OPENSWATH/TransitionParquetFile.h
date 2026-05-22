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

      @ingroup FileIO
  */
  class OPENMS_DLLAPI TransitionParquetFile
  {
  public:
    /// Default constructor (the class is stateless; API consists of the two convert* methods)
    TransitionParquetFile() = default;

    /// Default destructor
    ~TransitionParquetFile() = default;

    /**
      @brief Read a @c .oswpq library and populate a @ref OpenSwath::LightTargetedExperiment.

      Opens @c library/precursors.parquet and @c library/transitions.parquet from @p oswpq_dir
      (zip archive or extracted directory), validates each table against the
      @ref OSWPrecursorSchema / @ref OSWTransitionSchema in subset mode (extra columns are
      tolerated), and materialises the rows into @p targeted_exp.

      @p targeted_exp is **reset** to an empty @ref OpenSwath::LightTargetedExperiment at the
      start of the call — pre-existing contents are discarded rather than appended to.

      @param[in]  oswpq_dir    Path to a @c .oswpq archive or to an already-extracted directory containing @c library/.
      @param[out] targeted_exp Populated targeted experiment; cleared before being filled.
      @throws Exception::MissingInformation If a required parquet entry (precursors / transitions) cannot be located inside @p oswpq_dir.
      @throws Exception::InvalidValue If a loaded parquet table fails schema validation against the precursor / transition schema, or if any other required row-level field is missing during materialisation.
    */
    void convertParquetToTargetedExperiment(const String& oswpq_dir, OpenSwath::LightTargetedExperiment& targeted_exp) const;

    /**
      @brief Write a @ref OpenSwath::LightTargetedExperiment to a @c .oswpq library.

      Output target depends on @p oswpq_path: if it points at an existing directory, the
      parquet files are written directly under @c \<oswpq_path\>/library/. Otherwise the writer
      stages the layout in a temporary directory, then assembles a zip archive at
      @p oswpq_path via a @c .tmp staging archive that is renamed into place once complete.

      @param[in] oswpq_path  Destination — existing directory or zip-file path.
      @param[in] targeted_exp Library to serialise.
      @throws Exception::FileNotWritable If a generated file inside the staging area cannot be written.
      @throws Exception::MissingInformation If required per-row data (e.g. a transition lacking a peptide ref) is missing from @p targeted_exp.
      @throws Exception::InvalidValue If row-level invariants are violated (e.g. duplicate precursor ids, schema-incompatible values).
    */
    void convertLightTargetedExperimentToParquet(const String& oswpq_path, const OpenSwath::LightTargetedExperiment& targeted_exp) const;
  };

} // namespace OpenMS
