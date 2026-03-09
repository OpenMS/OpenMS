
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
    OpenSwathOSWParquetWriter (library/precursors.parquet, runs/runs.parquet and
    per-run features.parquet) and exposes a flat table of feature rows that
    combine feature-level scores with precursor/run metadata. This output is
    meant for downstream scoring workflows (PyProphet) or for simple single table exports.
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

    /// Result container for fetchPeakGroupFeatures
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

    /// Result container for transition-level features
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

  /// @brief Result container for an unscored table
  ///
  /// An "unscored" table contains feature-level columns but does not include
  /// discriminant scores or FDR estimates. Provides many columns (feature,
  /// precursor, run, MS1/MS2 metrics and discovered feature score columns).
    struct UnscoredResult
    {
      std::vector<int64_t> id_run;
      std::vector<int64_t> id_peptide; 
      std::vector<int64_t> transition_group_id; 
      std::vector<bool> decoy;
      std::vector<int64_t> run_id;
      std::vector<String> filename;
      std::vector<double> RT;
      std::vector<double> assay_rt;   // FEATURE.EXP_RT - FEATURE.DELTA_RT
      std::vector<double> delta_rt;   // FEATURE.DELTA_RT
      std::vector<double> assay_RT;   // PRECURSOR.LIBRARY_RT
      std::vector<double> delta_RT;   // FEATURE.NORM_RT - PRECURSOR.LIBRARY_RT
      std::vector<int64_t> id;        // FEATURE.ID
      std::vector<int> Charge;        // PRECURSOR.CHARGE
      std::vector<double> mz;         // PRECURSOR.PRECURSOR_MZ
      std::vector<double> Intensity;  // FEATURE_MS2.AREA_INTENSITY

      // aggregated MS1 metrics
      std::vector<double> aggr_prec_Peak_Area;
      std::vector<double> aggr_prec_Peak_Apex;

      std::vector<double> leftWidth;
      std::vector<double> rightWidth;

      // ion-mobility columns (may be NULL)
      std::vector<double> EXP_IM;
      std::vector<double> IM_leftWidth;
      std::vector<double> IM_rightWidth;

      // Discovered score columns (var_ms1_ and var_ms2_) and their values
      std::vector<String> ms2_columns;
      std::vector<std::vector<double>> ms2_values;
      std::vector<String> ms1_columns;
      std::vector<std::vector<double>> ms1_values;
    };

    /// Default constructor
    OpenSwathOSWParquetReader() = default;

    /**
      @brief Convenience constructor that loads from the given oswpq path

      This constructor calls load(oswpq_dir) so the object is ready to use
      after construction. It is provided for Python ergonomics so callers
      can create an instance with a single argument similarly to other
      Parquet helper classes.

      @param[in] oswpq_dir  Path to the unzipped OSW Parquet directory or a
                           .oswpq archive (zip) that will be read.
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

      This method reads per-run `features.parquet` files and returns a
      PeakGroupFeatureScoresResult containing feature identifiers, retention
      times, precursor metadata and discovered MS2 (and optionally MS1)
      score columns. The returned result is sorted by run_id, precursor_id and
      exp_rt to match the ordering used by the sqlite-based extractor.

      @param[in] oswpq_dir  Path to the unzipped OSW Parquet directory or a
                           .oswpq archive (zip) that will be read.
      @param[in] level       "ms2" (default) or "ms1ms2" to also include MS1 scores
      @param[in] main_score  Optional main score name to be used downstream. If provided and present among discovered MS2 columns, the specified
             column will be placed first in `ms2_columns`/`ms2_values`
             to ease downstream usage where a primary score is expected.
      @return PeakGroupFeatureScoresResult populated with discovered columns
              and core feature fields.
    */
    PeakGroupFeatureScoresResult fetchPeakGroupFeatures(const String& oswpq_dir, const String& level = "ms2", const String& main_score = "") const;

    /**
      @brief Extract transition-level feature rows across all runs (SOA)

      Reads per-run `feature_transition.parquet` files and joins them with
      the per-run `features.parquet` and library tables to produce a
      column-oriented (structure-of-arrays) result. The returned
      TransitionFeaturesResult contains per-transition columns including:

      - feature identifiers and run/precursor metadata: `feature_id`,
        `run_id`, `precursor_id`, `exp_rt`, `precursor_charge`
      - transition identifiers and properties: `transition_id`, `product_charge`,
        transition-level `decoy` flag
      - basic transition peak metrics: `area_intensity`, `total_area_intensity`,
        `apex_intensity`, `apex_rt`, `rt_fwhm`, `masserror_ppm`, `total_mi`
      - discovered transition-level score columns (e.g. `var_ms2_*`) are
        collected into `transition_var_columns` with their per-row values in
        `transition_var_values` (one vector per discovered column)
      - `group_id` strings in the form run_feature_precursor_transition to
        allow easy grouping or joins on the Python side.

      The ordering of rows is deterministic (sorted by run, precursor, feature,
      transition) to match other extractors.

      @param[in] oswpq_dir  Path to the unzipped OSW Parquet directory or
                           a .oswpq archive (zip) that will be read.
      @return TransitionFeaturesResult populated with transition-level columns
    */
    TransitionFeaturesResult fetchTransitionFeatures(const String& oswpq_dir) const;

    /**
      @brief Read an "unscored" table and return a column-oriented result.

      Reads per-run `features.parquet` (and supporting tables) and returns a
      rich set of per-feature columns in a Structure-of-Arrays layout
      convenient for conversion to a pandas.DataFrame. The returned
      UnscoredResult contains the following columns (vectors of length N):

      - feature/run identifiers: `id_run`, `id_peptide` (optional, 0 if unknown),
        `transition_group_id` (precursor id), `decoy`, `run_id`, `filename`
      - retention time fields: `RT` (feature EXP_RT), `assay_rt` (EXP_RT - DELTA_RT),
        `delta_rt` (FEATURE.DELTA_RT), `assay_RT` (precursor/library RT),
        `delta_RT` (norm_rt - library RT)
      - feature identifiers and precursor metadata: `id` (FEATURE.ID), `Charge`, `mz`
      - MS2/MS1 intensity metrics: `Intensity` (FEATURE_MS2.AREA_INTENSITY),
        `aggr_prec_Peak_Area`, `aggr_prec_Peak_Apex`
      - peak boundary widths: `left_width`, `right_width` (may be NaN if absent; available in the returned struct as `leftWidth`/`rightWidth`)
      - optional ion-mobility fields: canonical Parquet names `exp_im`, `exp_im_leftwidth`, `exp_im_rightwidth` (returned as `EXP_IM`/`IM_leftWidth`/`IM_rightWidth` in the struct)
      - discovered score columns: `ms2_columns` / `ms2_values` and
        `ms1_columns` / `ms1_values` (discovered across runs, e.g. `var_ms2_*`)

      Optional columns that are not present in a given Parquet file will be
      populated with default values (NaN for floating fields, empty strings for
      filenames, 0 for integer ids) so all returned vectors have equal length
      and are safe to convert into tabular formats.

      @param[in] oswpq_dir  Path to the unzipped OSW Parquet directory or a
                           .oswpq archive (zip) that will be read.
      @return UnscoredResult containing the assembled column vectors.
    */
    UnscoredResult fetchUnscoredData(const String& oswpq_dir) const;

  private:
    std::vector<Row> rows_;
    // store last-loaded path so Python-side code can call fetch methods without re-supplying the path
    String oswpq_dir_;
  };

} // namespace OpenMS
