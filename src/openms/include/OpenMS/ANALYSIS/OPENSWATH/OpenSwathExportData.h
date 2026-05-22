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
#include <OpenMS/config.h>

#include <optional>
#include <vector>

namespace OpenMS
{
  /**
    @brief Output file format for OpenSWATH export tables.

    @ingroup TargetedQuantitation
  */
  enum class OpenSwathExportFileFormat
  {
    TSV,
    Parquet
  };

  /**
    @brief IPF handling mode for OpenSWATH exports.

    @ingroup TargetedQuantitation
  */
  enum class OpenSwathIPFExportMode
  {
    /// Report results on peptidoform-level
    Peptidoform,
    /// Augment OpenSWATH results with IPF scores
    Augmented,
    /// Ignore IPF scores, report standard OpenSWATH results only
    Disable
  };

  /**
    @brief Quantification matrix level for OpenSWATH matrix exports.

    @ingroup TargetedQuantitation
  */
  enum class OpenSwathMatrixLevel
  {
    Precursor,
    Peptide,
    Protein,
    Gene
  };

  /**
    @brief Normalization method for OpenSWATH matrix exports.

    @ingroup TargetedQuantitation
  */
  enum class OpenSwathMatrixNormalization
  {
    None,
    Median,
    MedianMedian
  };

  inline String toString(const OpenSwathExportFileFormat format)
  {
    switch (format)
    {
      case OpenSwathExportFileFormat::TSV: return "tsv";
      case OpenSwathExportFileFormat::Parquet: return "parquet";
    }
    return "";
  }

  inline String toString(const OpenSwathIPFExportMode mode)
  {
    switch (mode)
    {
      case OpenSwathIPFExportMode::Peptidoform: return "peptidoform";
      case OpenSwathIPFExportMode::Augmented: return "augmented";
      case OpenSwathIPFExportMode::Disable: return "disable";
    }
    return "";
  }

  inline String toString(const OpenSwathMatrixLevel level)
  {
    switch (level)
    {
      case OpenSwathMatrixLevel::Precursor: return "precursor";
      case OpenSwathMatrixLevel::Peptide: return "peptide";
      case OpenSwathMatrixLevel::Protein: return "protein";
      case OpenSwathMatrixLevel::Gene: return "gene";
    }
    return "";
  }

  inline String toString(const OpenSwathMatrixNormalization normalization)
  {
    switch (normalization)
    {
      case OpenSwathMatrixNormalization::None: return "none";
      case OpenSwathMatrixNormalization::Median: return "median";
      case OpenSwathMatrixNormalization::MedianMedian: return "medianmedian";
    }
    return "";
  }

  /**
    @brief Filtered feature-level row used by results and matrix exports.

    This row intentionally stores only the columns needed by the user-facing
    results and matrix exports. The wider Parquet scoring export uses dedicated
    structs below because it may also carry dynamically discovered FEATURE_MS1 /
    FEATURE_MS2 / FEATURE_TRANSITION columns.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathExportRow
  {
    Int64 run_id = -1;
    String filename;
    String run_name;
    Int64 feature_id = -1;
    Int64 peptide_id = -1;
    Int64 precursor_id = -1;
    String transition_group_id;
    bool decoy = false;
    String sequence;
    String full_peptide_name;
    String protein_name;
    String gene_name;
    Int32 charge = 0;
    double mz = 0.0;
    double rt = 0.0;
    double assay_rt = 0.0;
    double delta_rt = 0.0;
    double irt = 0.0;
    double assay_irt = 0.0;
    double delta_irt = 0.0;
    double intensity = 0.0;
    std::optional<double> aggr_prec_peak_area;
    std::optional<double> aggr_prec_peak_apex;
    double left_width = 0.0;
    double right_width = 0.0;
    std::optional<double> exp_im;
    std::optional<double> im_left_width;
    std::optional<double> im_right_width;
    std::optional<double> ms1_pep;
    std::optional<double> ms2_pep;
    std::optional<double> precursor_pep;
    std::optional<double> ipf_pep;
    Int32 peak_group_rank = 0;
    double d_score = 0.0;
    double m_score = 0.0;
    std::optional<double> pep;
    std::optional<double> ms2_m_score;
    std::optional<double> peptide_run_specific_qvalue;
    std::optional<double> peptide_experiment_wide_qvalue;
    std::optional<double> peptide_global_qvalue;
    std::optional<double> protein_run_specific_qvalue;
    std::optional<double> protein_experiment_wide_qvalue;
    std::optional<double> protein_global_qvalue;
    std::optional<double> gene_run_specific_qvalue;
    std::optional<double> gene_experiment_wide_qvalue;
    std::optional<double> gene_global_qvalue;
    std::optional<Int64> alignment_group_id;
    std::optional<Int64> alignment_reference_feature_id;
    std::optional<double> alignment_reference_rt;
    std::optional<double> alignment_pep;
    std::optional<double> alignment_qvalue;
    bool from_alignment = false;
    String aggr_peak_area;
    String aggr_peak_apex;
    String aggr_fragment_annotation;
    String ipf_full_peptide_name;
    std::optional<double> ipf_precursor_peakgroup_pep;
    std::optional<double> ipf_peptidoform_pep;
    std::optional<double> ipf_peptidoform_m_score;
  };

  /**
    @brief Feature-level scored row for Parquet export.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathFeatureScoreRow
  {
    Int64 protein_id = -1;
    Int64 peptide_id = -1;
    std::optional<Int64> ipf_peptide_id;
    Int64 precursor_id = -1;
    String protein_accession;
    String unmodified_sequence;
    String modified_sequence;
    String precursor_traml_id;
    String precursor_group_label;
    double precursor_mz = 0.0;
    Int32 precursor_charge = 0;
    std::optional<double> precursor_library_intensity;
    std::optional<double> precursor_library_rt;
    std::optional<double> precursor_library_drift_time;
    std::optional<Int64> gene_id;
    String gene_name;
    std::optional<bool> gene_decoy;
    bool protein_decoy = false;
    bool peptide_decoy = false;
    bool precursor_decoy = false;
    Int64 run_id = -1;
    String filename;
    Int64 feature_id = -1;
    double exp_rt = 0.0;
    std::optional<double> exp_im;
    double norm_rt = 0.0;
    double delta_rt = 0.0;
    double left_width = 0.0;
    double right_width = 0.0;
    std::optional<double> im_left_width;
    std::optional<double> im_right_width;
    std::vector<double> feature_ms1_values;
    std::vector<double> feature_ms2_values;
    std::optional<double> score_ms1_score;
    std::optional<Int32> score_ms1_rank;
    std::optional<double> score_ms1_pvalue;
    std::optional<double> score_ms1_qvalue;
    std::optional<double> score_ms1_pep;
    std::optional<double> score_ms2_score;
    std::optional<Int32> score_ms2_peak_group_rank;
    std::optional<double> score_ms2_pvalue;
    std::optional<double> score_ms2_qvalue;
    std::optional<double> score_ms2_pep;
    std::optional<double> score_ipf_precursor_peakgroup_pep;
    std::optional<double> score_ipf_pep;
    std::optional<double> score_ipf_qvalue;
    std::optional<double> score_peptide_global_score;
    std::optional<double> score_peptide_global_pvalue;
    std::optional<double> score_peptide_global_qvalue;
    std::optional<double> score_peptide_global_pep;
    std::optional<double> score_peptide_experiment_wide_score;
    std::optional<double> score_peptide_experiment_wide_pvalue;
    std::optional<double> score_peptide_experiment_wide_qvalue;
    std::optional<double> score_peptide_experiment_wide_pep;
    std::optional<double> score_peptide_run_specific_score;
    std::optional<double> score_peptide_run_specific_pvalue;
    std::optional<double> score_peptide_run_specific_qvalue;
    std::optional<double> score_peptide_run_specific_pep;
    std::optional<double> score_protein_global_score;
    std::optional<double> score_protein_global_pvalue;
    std::optional<double> score_protein_global_qvalue;
    std::optional<double> score_protein_global_pep;
    std::optional<double> score_protein_experiment_wide_score;
    std::optional<double> score_protein_experiment_wide_pvalue;
    std::optional<double> score_protein_experiment_wide_qvalue;
    std::optional<double> score_protein_experiment_wide_pep;
    std::optional<double> score_protein_run_specific_score;
    std::optional<double> score_protein_run_specific_pvalue;
    std::optional<double> score_protein_run_specific_qvalue;
    std::optional<double> score_protein_run_specific_pep;
    std::optional<double> score_gene_global_score;
    std::optional<double> score_gene_global_pvalue;
    std::optional<double> score_gene_global_qvalue;
    std::optional<double> score_gene_global_pep;
    std::optional<double> score_gene_experiment_wide_score;
    std::optional<double> score_gene_experiment_wide_pvalue;
    std::optional<double> score_gene_experiment_wide_qvalue;
    std::optional<double> score_gene_experiment_wide_pep;
    std::optional<double> score_gene_run_specific_score;
    std::optional<double> score_gene_run_specific_pvalue;
    std::optional<double> score_gene_run_specific_qvalue;
    std::optional<double> score_gene_run_specific_pep;
    std::optional<Int64> alignment_group_id;
    std::optional<Int64> alignment_reference_feature_id;
    std::optional<double> alignment_reference_rt;
    std::optional<double> alignment_pep;
    std::optional<double> alignment_qvalue;
    bool from_alignment = false;
  };

  /**
    @brief Feature-level Parquet export table with dynamically discovered FEATURE_MS1 / FEATURE_MS2 columns.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathFeatureScoreTable
  {
    std::vector<String> feature_ms1_column_names;
    std::vector<String> feature_ms2_column_names;
    std::vector<OpenSwathFeatureScoreRow> rows;
  };

  /**
    @brief Transition-level scored row for optional Parquet export.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathTransitionScoreRow
  {
    std::optional<Int64> run_id;
    std::optional<Int64> ipf_peptide_id;
    Int64 precursor_id = -1;
    Int64 transition_id = -1;
    String transition_traml_id;
    double product_mz = 0.0;
    Int32 transition_charge = 0;
    String transition_type;
    Int32 transition_ordinal = 0;
    String annotation;
    bool transition_detecting = false;
    std::optional<double> transition_library_intensity;
    bool transition_decoy = false;
    std::optional<Int64> feature_id;
    std::vector<std::optional<double>> feature_transition_values;
    std::optional<double> score_transition_score;
    std::optional<Int32> score_transition_rank;
    std::optional<double> score_transition_pvalue;
    std::optional<double> score_transition_qvalue;
    std::optional<double> score_transition_pep;
  };

  /**
    @brief Transition-level Parquet export table with dynamically discovered FEATURE_TRANSITION columns.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathTransitionScoreTable
  {
    std::vector<String> feature_transition_column_names;
    std::vector<OpenSwathTransitionScoreRow> rows;
  };

  /**
    @brief Dense matrix representation used by matrix exports.

    Numeric values are stored by row with missing values represented as
    std::nullopt to preserve sparse run coverage in TSV and Parquet output.

    @ingroup TargetedQuantitation
  */
  struct OPENMS_DLLAPI OpenSwathQuantMatrix
  {
    std::vector<String> identifier_column_names;
    std::vector<std::vector<String>> identifier_rows;
    std::vector<String> sample_column_names;
    std::vector<std::vector<std::optional<double>>> values;
  };
} // namespace OpenMS
