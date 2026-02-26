
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
#include <map>
#include <vector>

namespace OpenMS
{
  /**
    @brief Reader for OpenSwath OSW Parquet output.

    This class reads the Parquet output layout produced by
    OpenSwathOSWParquetWriter (library/*.parquet, runs/runs.parquet and
    per-run features.parquet) and exposes a flat table of feature rows that
    combine feature-level scores with precursor/run metadata. The output is
    suitable for downstream scoring workflows (for example pyProphet).

    The reader focuses on a compact, commonly-used column set (feature id,
    run id, precursor id, RT, MS2 area/apex, precursor charge/decoy and the
    number of detecting transitions). Additional columns can be added easily
    if required.

    @ingroup Analysis
  */
  class OPENMS_DLLAPI OpenSwathOSWParquetReader
  {
  public:
    /// Single extracted row combining feature + precursor + run metadata
    struct Row
    {
      int64_t feature_id = 0;
      int64_t run_id = 0;
      int64_t precursor_id = 0;
      double exp_rt = 0.0;
      double ms2_area_intensity = 0.0;
      double ms2_total_area_intensity = 0.0;
      double ms2_apex_intensity = 0.0;
      int precursor_charge = 0;
      bool decoy = false;
      int64_t transition_count = 0;
      String group_id;
    };

  /// MS2 row declared at file scope for pybind/pyx friendliness
  struct OpenSwathOSWParquetReaderRowMS2
    {
      int64_t feature_id = 0;
      int64_t run_id = 0;
      int64_t precursor_id = 0;
      double exp_rt = 0.0;
      int precursor_charge = 0;
      bool decoy = false;
      int64_t transition_count = 0;
      String group_id;

      // Dynamic score maps (keyed by column name without prefix)
      std::map<String, double> ms2_scores; // keys are e.g. "var_ms2_dotprod_score"
      std::map<String, double> ms1_scores; // keys are e.g. "var_ms1_mi_score"
    };

    /// Default constructor
    OpenSwathOSWParquetReader() = default;

    /**
      @brief Convenience constructor that loads from the given oswpq path

      This constructor calls load(oswpq_dir) so the object is ready to use
      after construction. It is provided for Python ergonomics so callers
      can create an instance with a single argument similarly to other
      Parquet helper classes.
    */
    OpenSwathOSWParquetReader(const String& oswpq_dir);

    /**
      @brief Load and extract rows from an OSW Parquet directory or .oswpq archive

      @param[in] oswpq_dir  Path to the unzipped directory or .oswpq archive
    */
  void load(const String& oswpq_dir);

  /// Return the originally provided oswpq path (may be empty)
  const String& oswpqPath() const { return oswpq_dir_; }

    /// Return extracted rows
    const std::vector<Row>& rows() const { return rows_; }

    /**
      @brief Write extracted rows to a simple CSV file (no quoting)

      The CSV contains a header and the columns:
      feature_id,run_id,precursor_id,exp_rt,ms2_area,ms2_total_area,ms2_apex,
      precursor_charge,decoy,transition_count
    */
    void writeCSV(const String& filename) const;

    /**
      @brief Extract MS2-level feature rows across all runs.

      This method reads per-run features.parquet files and returns a vector of
      RowMS2 containing feature identifiers, RT, precursor metadata and maps of
      MS2 (and optionally MS1) score columns. The returned vector is sorted by
      run_id, precursor_id and exp_rt to match the SQL ordering used in the
      sqlite-based extractor.

      @param[out] out_rows   Vector filled with extracted rows
      @param[in] level       "ms2" (default) or "ms1ms2" to also include MS1 scores
      @param[in] main_score  Optional main score name to be used downstream
    */
  void fetchMS2Features(const String& oswpq_dir, std::vector<OpenSwathOSWParquetReaderRowMS2>& out_rows, const String& level = "ms2", const String& main_score = "") const;

    /**
      @brief Write a flattened Parquet file of the MS2 (and optionally MS1) feature table.

      The output Parquet file contains one row per feature and expands all
      score columns discovered in the per-run features.parquet files into
      top-level columns. This method is intended to be used by Python
      bindings which can then load the Parquet file using pandas/pyarrow
      to obtain an identical table shape to the sqlite-based extractor.

      @param[in] oswpq_dir    Path to the unzipped directory or .oswpq archive
      @param[in] output_path  Path where the parquet file will be written
      @param[in] level        "ms2" (default) or "ms1ms2" to include MS1 score columns
    */
    void writeMS2FeaturesParquet(const String& oswpq_dir, const String& output_path, const String& level = "ms2") const;
  private:
    std::vector<Row> rows_;
    // store last-loaded path so Python-side code can call fetch methods without re-supplying the path
    String oswpq_dir_;
  };

} // namespace OpenMS
