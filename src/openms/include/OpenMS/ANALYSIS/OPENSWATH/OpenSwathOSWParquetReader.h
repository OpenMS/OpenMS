
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

    /// Result container for fetchPeakGroupFeatures: column-oriented (SOA) layout
    /// Contains discovered MS2 and optional MS1 score columns alongside core feature columns.
    struct PeakGroupFeatureScoresResult
    {
      // Core per-feature columns (all length N)
      std::vector<int64_t> feature_id;
      std::vector<int64_t> run_id;
      std::vector<int64_t> precursor_id;
      std::vector<double> exp_rt;
      std::vector<int> precursor_charge;
      std::vector<bool> decoy;
      std::vector<int64_t> transition_count;
      std::vector<String> group_id;

      // Discovered MS2 score columns and their column vectors (each vector length N)
      // ms2_columns[i] corresponds to ms2_values[i]
      std::vector<String> ms2_columns;
      std::vector<std::vector<double>> ms2_values;

      // Discovered MS1 score columns and their column vectors (only present if requested)
      std::vector<String> ms1_columns;
      std::vector<std::vector<double>> ms1_values;
    };

    /// Result container for transition-level features (SOA)
    struct TransitionFeaturesResult
    {
      // Core per-transition columns (all length N)
      std::vector<int64_t> feature_id;
      std::vector<int64_t> run_id;
      std::vector<int64_t> precursor_id;
      std::vector<double> exp_rt;
      std::vector<int> precursor_charge;

      std::vector<int64_t> transition_id;
      std::vector<int> product_charge;
      std::vector<bool> decoy; // transition-level decoy flag

      // basic transition-level peak metrics
      std::vector<double> area_intensity;
      std::vector<double> total_area_intensity;
      std::vector<double> apex_intensity;
      std::vector<double> apex_rt;
      std::vector<double> rt_fwhm;
      std::vector<double> masserror_ppm;
      std::vector<double> total_mi;

      // Discovered transition-level var_ columns and their column vectors
      std::vector<String> transition_var_columns;
      std::vector<std::vector<double>> transition_var_values;

      // Group id (run_feature_precursor_transition style)
      std::vector<String> group_id;
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
    PeakGroupFeatureScoresResult fetchPeakGroupFeatures(const String& oswpq_dir, const String& level = "ms2", const String& main_score = "") const;

    /// Extract transition-level feature rows across all runs (SOA)
    /// Reads per-run feature_transition.parquet and joins to features/library to
    /// provide transition-level metrics alongside precursor/feature metadata.
    TransitionFeaturesResult fetchTransitionFeatures(const String& oswpq_dir) const;

  private:
    std::vector<Row> rows_;
    // store last-loaded path so Python-side code can call fetch methods without re-supplying the path
    String oswpq_dir_;
  };

} // namespace OpenMS
