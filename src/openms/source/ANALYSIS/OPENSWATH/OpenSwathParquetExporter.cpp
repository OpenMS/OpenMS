// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathParquetExporter.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <arrow/api.h>

namespace OpenMS
{
  namespace
  {
    void appendOptionalDouble_(arrow::DoubleBuilder& builder, const std::optional<double>& value)
    {
      if (value.has_value()) builder.Append(*value);
      else builder.AppendNull();
    }

    void appendOptionalInt64_(arrow::Int64Builder& builder, const std::optional<Int64>& value)
    {
      if (value.has_value()) builder.Append(*value);
      else builder.AppendNull();
    }

    void appendOptionalBool_(arrow::BooleanBuilder& builder, const std::optional<bool>& value)
    {
      if (value.has_value()) builder.Append(*value);
      else builder.AppendNull();
    }

    void appendString_(arrow::StringBuilder& builder, const String& value)
    {
      if (value.empty()) builder.AppendNull();
      else builder.Append(value.c_str());
    }
  }

  void OpenSwathParquetExporter::writeFeatureScores(const String& filename,
                                                    const OpenSwathFeatureScoreTable& table)
  {
    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;

    arrow::Int64Builder protein_id_b, peptide_id_b, precursor_id_b, run_id_b, feature_id_b;
    arrow::StringBuilder prot_b, seq_b, fullpep_b, traml_b, group_label_b, file_b, gene_name_b;
    arrow::DoubleBuilder precursor_mz_b, lib_intensity_b, lib_rt_b, lib_dt_b, exp_rt_b, exp_im_b, norm_rt_b, delta_rt_b, left_b, right_b, im_left_b, im_right_b;
    arrow::Int64Builder precursor_charge_b;
    arrow::BooleanBuilder protein_decoy_b, peptide_decoy_b, precursor_decoy_b;
    arrow::Int64Builder ipf_peptide_id_b, gene_id_b, alignment_group_id_b, alignment_reference_id_b;
    arrow::BooleanBuilder gene_decoy_b, from_alignment_b;
    arrow::DoubleBuilder score_ms1_score_b, score_ms1_pvalue_b, score_ms1_qvalue_b, score_ms1_pep_b;
    arrow::Int64Builder score_ms1_rank_b, score_ms2_rank_b;
    arrow::DoubleBuilder score_ms2_score_b, score_ms2_pvalue_b, score_ms2_qvalue_b, score_ms2_pep_b, score_ipf_prec_b, score_ipf_pep_b, score_ipf_qvalue_b;
    arrow::DoubleBuilder pep_global_score_b, pep_global_pvalue_b, pep_global_qvalue_b, pep_global_pep_b;
    arrow::DoubleBuilder pep_ew_score_b, pep_ew_pvalue_b, pep_ew_qvalue_b, pep_ew_pep_b;
    arrow::DoubleBuilder pep_rs_score_b, pep_rs_pvalue_b, pep_rs_qvalue_b, pep_rs_pep_b;
    arrow::DoubleBuilder prot_global_score_b, prot_global_pvalue_b, prot_global_qvalue_b, prot_global_pep_b;
    arrow::DoubleBuilder prot_ew_score_b, prot_ew_pvalue_b, prot_ew_qvalue_b, prot_ew_pep_b;
    arrow::DoubleBuilder prot_rs_score_b, prot_rs_pvalue_b, prot_rs_qvalue_b, prot_rs_pep_b;
    arrow::DoubleBuilder gene_global_score_b, gene_global_pvalue_b, gene_global_qvalue_b, gene_global_pep_b;
    arrow::DoubleBuilder gene_ew_score_b, gene_ew_pvalue_b, gene_ew_qvalue_b, gene_ew_pep_b;
    arrow::DoubleBuilder gene_rs_score_b, gene_rs_pvalue_b, gene_rs_qvalue_b, gene_rs_pep_b;
    arrow::DoubleBuilder align_rt_b, align_pep_b, align_qvalue_b;

    std::vector<arrow::DoubleBuilder> ms1_builders(table.feature_ms1_column_names.size());
    std::vector<arrow::DoubleBuilder> ms2_builders(table.feature_ms2_column_names.size());

    for (const auto& row : table.rows)
    {
      protein_id_b.Append(row.protein_id);
      peptide_id_b.Append(row.peptide_id);
      appendOptionalInt64_(ipf_peptide_id_b, row.ipf_peptide_id);
      precursor_id_b.Append(row.precursor_id);
      appendString_(prot_b, row.protein_accession);
      appendString_(seq_b, row.unmodified_sequence);
      appendString_(fullpep_b, row.modified_sequence);
      appendString_(traml_b, row.precursor_traml_id);
      appendString_(group_label_b, row.precursor_group_label);
      precursor_mz_b.Append(row.precursor_mz);
      precursor_charge_b.Append(row.precursor_charge);
      appendOptionalDouble_(lib_intensity_b, row.precursor_library_intensity);
      appendOptionalDouble_(lib_rt_b, row.precursor_library_rt);
      appendOptionalDouble_(lib_dt_b, row.precursor_library_drift_time);
      appendOptionalInt64_(gene_id_b, row.gene_id);
      appendString_(gene_name_b, row.gene_name);
      appendOptionalBool_(gene_decoy_b, row.gene_decoy);
      protein_decoy_b.Append(row.protein_decoy);
      peptide_decoy_b.Append(row.peptide_decoy);
      precursor_decoy_b.Append(row.precursor_decoy);
      run_id_b.Append(row.run_id);
      appendString_(file_b, row.filename);
      feature_id_b.Append(row.feature_id);
      exp_rt_b.Append(row.exp_rt);
      appendOptionalDouble_(exp_im_b, row.exp_im);
      norm_rt_b.Append(row.norm_rt);
      delta_rt_b.Append(row.delta_rt);
      left_b.Append(row.left_width);
      right_b.Append(row.right_width);
      appendOptionalDouble_(im_left_b, row.im_left_width);
      appendOptionalDouble_(im_right_b, row.im_right_width);
      for (Size i = 0; i < table.feature_ms1_column_names.size(); ++i)
      {
        if (i < row.feature_ms1_values.size()) ms1_builders[i].Append(row.feature_ms1_values[i]);
        else ms1_builders[i].AppendNull();
      }
      for (Size i = 0; i < table.feature_ms2_column_names.size(); ++i)
      {
        if (i < row.feature_ms2_values.size()) ms2_builders[i].Append(row.feature_ms2_values[i]);
        else ms2_builders[i].AppendNull();
      }
      appendOptionalDouble_(score_ms1_score_b, row.score_ms1_score);
      if (row.score_ms1_rank.has_value()) score_ms1_rank_b.Append(*row.score_ms1_rank);
      else score_ms1_rank_b.AppendNull();
      appendOptionalDouble_(score_ms1_pvalue_b, row.score_ms1_pvalue);
      appendOptionalDouble_(score_ms1_qvalue_b, row.score_ms1_qvalue);
      appendOptionalDouble_(score_ms1_pep_b, row.score_ms1_pep);
      appendOptionalDouble_(score_ms2_score_b, row.score_ms2_score);
      if (row.score_ms2_peak_group_rank.has_value()) score_ms2_rank_b.Append(*row.score_ms2_peak_group_rank);
      else score_ms2_rank_b.AppendNull();
      appendOptionalDouble_(score_ms2_pvalue_b, row.score_ms2_pvalue);
      appendOptionalDouble_(score_ms2_qvalue_b, row.score_ms2_qvalue);
      appendOptionalDouble_(score_ms2_pep_b, row.score_ms2_pep);
      appendOptionalDouble_(score_ipf_prec_b, row.score_ipf_precursor_peakgroup_pep);
      appendOptionalDouble_(score_ipf_pep_b, row.score_ipf_pep);
      appendOptionalDouble_(score_ipf_qvalue_b, row.score_ipf_qvalue);
      appendOptionalDouble_(pep_global_score_b, row.score_peptide_global_score);
      appendOptionalDouble_(pep_global_pvalue_b, row.score_peptide_global_pvalue);
      appendOptionalDouble_(pep_global_qvalue_b, row.score_peptide_global_qvalue);
      appendOptionalDouble_(pep_global_pep_b, row.score_peptide_global_pep);
      appendOptionalDouble_(pep_ew_score_b, row.score_peptide_experiment_wide_score);
      appendOptionalDouble_(pep_ew_pvalue_b, row.score_peptide_experiment_wide_pvalue);
      appendOptionalDouble_(pep_ew_qvalue_b, row.score_peptide_experiment_wide_qvalue);
      appendOptionalDouble_(pep_ew_pep_b, row.score_peptide_experiment_wide_pep);
      appendOptionalDouble_(pep_rs_score_b, row.score_peptide_run_specific_score);
      appendOptionalDouble_(pep_rs_pvalue_b, row.score_peptide_run_specific_pvalue);
      appendOptionalDouble_(pep_rs_qvalue_b, row.score_peptide_run_specific_qvalue);
      appendOptionalDouble_(pep_rs_pep_b, row.score_peptide_run_specific_pep);
      appendOptionalDouble_(prot_global_score_b, row.score_protein_global_score);
      appendOptionalDouble_(prot_global_pvalue_b, row.score_protein_global_pvalue);
      appendOptionalDouble_(prot_global_qvalue_b, row.score_protein_global_qvalue);
      appendOptionalDouble_(prot_global_pep_b, row.score_protein_global_pep);
      appendOptionalDouble_(prot_ew_score_b, row.score_protein_experiment_wide_score);
      appendOptionalDouble_(prot_ew_pvalue_b, row.score_protein_experiment_wide_pvalue);
      appendOptionalDouble_(prot_ew_qvalue_b, row.score_protein_experiment_wide_qvalue);
      appendOptionalDouble_(prot_ew_pep_b, row.score_protein_experiment_wide_pep);
      appendOptionalDouble_(prot_rs_score_b, row.score_protein_run_specific_score);
      appendOptionalDouble_(prot_rs_pvalue_b, row.score_protein_run_specific_pvalue);
      appendOptionalDouble_(prot_rs_qvalue_b, row.score_protein_run_specific_qvalue);
      appendOptionalDouble_(prot_rs_pep_b, row.score_protein_run_specific_pep);
      appendOptionalDouble_(gene_global_score_b, row.score_gene_global_score);
      appendOptionalDouble_(gene_global_pvalue_b, row.score_gene_global_pvalue);
      appendOptionalDouble_(gene_global_qvalue_b, row.score_gene_global_qvalue);
      appendOptionalDouble_(gene_global_pep_b, row.score_gene_global_pep);
      appendOptionalDouble_(gene_ew_score_b, row.score_gene_experiment_wide_score);
      appendOptionalDouble_(gene_ew_pvalue_b, row.score_gene_experiment_wide_pvalue);
      appendOptionalDouble_(gene_ew_qvalue_b, row.score_gene_experiment_wide_qvalue);
      appendOptionalDouble_(gene_ew_pep_b, row.score_gene_experiment_wide_pep);
      appendOptionalDouble_(gene_rs_score_b, row.score_gene_run_specific_score);
      appendOptionalDouble_(gene_rs_pvalue_b, row.score_gene_run_specific_pvalue);
      appendOptionalDouble_(gene_rs_qvalue_b, row.score_gene_run_specific_qvalue);
      appendOptionalDouble_(gene_rs_pep_b, row.score_gene_run_specific_pep);
      appendOptionalInt64_(alignment_group_id_b, row.alignment_group_id);
      appendOptionalInt64_(alignment_reference_id_b, row.alignment_reference_feature_id);
      appendOptionalDouble_(align_rt_b, row.alignment_reference_rt);
      appendOptionalDouble_(align_pep_b, row.alignment_pep);
      appendOptionalDouble_(align_qvalue_b, row.alignment_qvalue);
      from_alignment_b.Append(row.from_alignment);
    }

    auto pushArray = [&](const char* name, const std::shared_ptr<arrow::DataType>& type, std::shared_ptr<arrow::Array> array)
    {
      fields.push_back(arrow::field(name, type));
      arrays.push_back(std::move(array));
    };

    pushArray("PROTEIN_ID", arrow::int64(), protein_id_b.Finish().ValueOrDie());
    pushArray("PEPTIDE_ID", arrow::int64(), peptide_id_b.Finish().ValueOrDie());
    pushArray("IPF_PEPTIDE_ID", arrow::int64(), ipf_peptide_id_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_ID", arrow::int64(), precursor_id_b.Finish().ValueOrDie());
    pushArray("PROTEIN_ACCESSION", arrow::utf8(), prot_b.Finish().ValueOrDie());
    pushArray("UNMODIFIED_SEQUENCE", arrow::utf8(), seq_b.Finish().ValueOrDie());
    pushArray("MODIFIED_SEQUENCE", arrow::utf8(), fullpep_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_TRAML_ID", arrow::utf8(), traml_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_GROUP_LABEL", arrow::utf8(), group_label_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_MZ", arrow::float64(), precursor_mz_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_CHARGE", arrow::int64(), precursor_charge_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_LIBRARY_INTENSITY", arrow::float64(), lib_intensity_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_LIBRARY_RT", arrow::float64(), lib_rt_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_LIBRARY_DRIFT_TIME", arrow::float64(), lib_dt_b.Finish().ValueOrDie());
    pushArray("GENE_ID", arrow::int64(), gene_id_b.Finish().ValueOrDie());
    pushArray("GENE_NAME", arrow::utf8(), gene_name_b.Finish().ValueOrDie());
    pushArray("GENE_DECOY", arrow::boolean(), gene_decoy_b.Finish().ValueOrDie());
    pushArray("PROTEIN_DECOY", arrow::boolean(), protein_decoy_b.Finish().ValueOrDie());
    pushArray("PEPTIDE_DECOY", arrow::boolean(), peptide_decoy_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_DECOY", arrow::boolean(), precursor_decoy_b.Finish().ValueOrDie());
    pushArray("RUN_ID", arrow::int64(), run_id_b.Finish().ValueOrDie());
    pushArray("FILENAME", arrow::utf8(), file_b.Finish().ValueOrDie());
    pushArray("FEATURE_ID", arrow::int64(), feature_id_b.Finish().ValueOrDie());
    pushArray("EXP_RT", arrow::float64(), exp_rt_b.Finish().ValueOrDie());
    pushArray("EXP_IM", arrow::float64(), exp_im_b.Finish().ValueOrDie());
    pushArray("NORM_RT", arrow::float64(), norm_rt_b.Finish().ValueOrDie());
    pushArray("DELTA_RT", arrow::float64(), delta_rt_b.Finish().ValueOrDie());
    pushArray("LEFT_WIDTH", arrow::float64(), left_b.Finish().ValueOrDie());
    pushArray("RIGHT_WIDTH", arrow::float64(), right_b.Finish().ValueOrDie());
    pushArray("IM_leftWidth", arrow::float64(), im_left_b.Finish().ValueOrDie());
    pushArray("IM_rightWidth", arrow::float64(), im_right_b.Finish().ValueOrDie());
    for (Size i = 0; i < table.feature_ms1_column_names.size(); ++i)
    {
      pushArray(("FEATURE_MS1_" + table.feature_ms1_column_names[i]).c_str(), arrow::float64(), ms1_builders[i].Finish().ValueOrDie());
    }
    for (Size i = 0; i < table.feature_ms2_column_names.size(); ++i)
    {
      pushArray(("FEATURE_MS2_" + table.feature_ms2_column_names[i]).c_str(), arrow::float64(), ms2_builders[i].Finish().ValueOrDie());
    }
    pushArray("SCORE_MS1_SCORE", arrow::float64(), score_ms1_score_b.Finish().ValueOrDie());
    pushArray("SCORE_MS1_RANK", arrow::int64(), score_ms1_rank_b.Finish().ValueOrDie());
    pushArray("SCORE_MS1_PVALUE", arrow::float64(), score_ms1_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_MS1_QVALUE", arrow::float64(), score_ms1_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_MS1_PEP", arrow::float64(), score_ms1_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_MS2_SCORE", arrow::float64(), score_ms2_score_b.Finish().ValueOrDie());
    pushArray("SCORE_MS2_PEAK_GROUP_RANK", arrow::int64(), score_ms2_rank_b.Finish().ValueOrDie());
    pushArray("SCORE_MS2_PVALUE", arrow::float64(), score_ms2_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_MS2_QVALUE", arrow::float64(), score_ms2_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_MS2_PEP", arrow::float64(), score_ms2_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_IPF_PRECURSOR_PEAKGROUP_PEP", arrow::float64(), score_ipf_prec_b.Finish().ValueOrDie());
    pushArray("SCORE_IPF_PEP", arrow::float64(), score_ipf_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_IPF_QVALUE", arrow::float64(), score_ipf_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_GLOBAL_SCORE", arrow::float64(), pep_global_score_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_GLOBAL_PVALUE", arrow::float64(), pep_global_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_GLOBAL_QVALUE", arrow::float64(), pep_global_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_GLOBAL_PEP", arrow::float64(), pep_global_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_EXPERIMENT_WIDE_SCORE", arrow::float64(), pep_ew_score_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_EXPERIMENT_WIDE_PVALUE", arrow::float64(), pep_ew_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_EXPERIMENT_WIDE_QVALUE", arrow::float64(), pep_ew_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_EXPERIMENT_WIDE_PEP", arrow::float64(), pep_ew_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_RUN_SPECIFIC_SCORE", arrow::float64(), pep_rs_score_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_RUN_SPECIFIC_PVALUE", arrow::float64(), pep_rs_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_RUN_SPECIFIC_QVALUE", arrow::float64(), pep_rs_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PEPTIDE_RUN_SPECIFIC_PEP", arrow::float64(), pep_rs_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_GLOBAL_SCORE", arrow::float64(), prot_global_score_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_GLOBAL_PVALUE", arrow::float64(), prot_global_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_GLOBAL_QVALUE", arrow::float64(), prot_global_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_GLOBAL_PEP", arrow::float64(), prot_global_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_EXPERIMENT_WIDE_SCORE", arrow::float64(), prot_ew_score_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_EXPERIMENT_WIDE_PVALUE", arrow::float64(), prot_ew_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_EXPERIMENT_WIDE_QVALUE", arrow::float64(), prot_ew_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_EXPERIMENT_WIDE_PEP", arrow::float64(), prot_ew_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_RUN_SPECIFIC_SCORE", arrow::float64(), prot_rs_score_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_RUN_SPECIFIC_PVALUE", arrow::float64(), prot_rs_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_RUN_SPECIFIC_QVALUE", arrow::float64(), prot_rs_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_PROTEIN_RUN_SPECIFIC_PEP", arrow::float64(), prot_rs_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_GLOBAL_SCORE", arrow::float64(), gene_global_score_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_GLOBAL_PVALUE", arrow::float64(), gene_global_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_GLOBAL_QVALUE", arrow::float64(), gene_global_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_GLOBAL_PEP", arrow::float64(), gene_global_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_EXPERIMENT_WIDE_SCORE", arrow::float64(), gene_ew_score_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_EXPERIMENT_WIDE_PVALUE", arrow::float64(), gene_ew_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_EXPERIMENT_WIDE_QVALUE", arrow::float64(), gene_ew_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_EXPERIMENT_WIDE_PEP", arrow::float64(), gene_ew_pep_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_RUN_SPECIFIC_SCORE", arrow::float64(), gene_rs_score_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_RUN_SPECIFIC_PVALUE", arrow::float64(), gene_rs_pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_RUN_SPECIFIC_QVALUE", arrow::float64(), gene_rs_qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_GENE_RUN_SPECIFIC_PEP", arrow::float64(), gene_rs_pep_b.Finish().ValueOrDie());
    pushArray("alignment_group_id", arrow::int64(), alignment_group_id_b.Finish().ValueOrDie());
    pushArray("alignment_reference_feature_id", arrow::int64(), alignment_reference_id_b.Finish().ValueOrDie());
    pushArray("alignment_reference_rt", arrow::float64(), align_rt_b.Finish().ValueOrDie());
    pushArray("alignment_pep", arrow::float64(), align_pep_b.Finish().ValueOrDie());
    pushArray("alignment_qvalue", arrow::float64(), align_qvalue_b.Finish().ValueOrDie());
    pushArray("from_alignment", arrow::boolean(), from_alignment_b.Finish().ValueOrDie());

    const auto arrow_table = arrow::Table::Make(arrow::schema(fields), arrays);
    if (!ArrowIOHelpers::writeTableToParquet(arrow_table, filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  void OpenSwathParquetExporter::writeTransitionScores(const String& filename,
                                                       const OpenSwathTransitionScoreTable& table)
  {
    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;

    arrow::Int64Builder run_id_b, ipf_peptide_id_b, precursor_id_b, transition_id_b, charge_b, ordinal_b, feature_id_b, score_rank_b;
    arrow::StringBuilder traml_b, type_b, annotation_b;
    arrow::DoubleBuilder product_mz_b, lib_intensity_b, score_b, pvalue_b, qvalue_b, pep_b;
    arrow::BooleanBuilder detecting_b, decoy_b;

    std::vector<arrow::DoubleBuilder> ft_builders(table.feature_transition_column_names.size());

    for (const auto& row : table.rows)
    {
      appendOptionalInt64_(run_id_b, row.run_id);
      appendOptionalInt64_(ipf_peptide_id_b, row.ipf_peptide_id);
      precursor_id_b.Append(row.precursor_id);
      transition_id_b.Append(row.transition_id);
      appendString_(traml_b, row.transition_traml_id);
      product_mz_b.Append(row.product_mz);
      charge_b.Append(row.transition_charge);
      appendString_(type_b, row.transition_type);
      ordinal_b.Append(row.transition_ordinal);
      appendString_(annotation_b, row.annotation);
      detecting_b.Append(row.transition_detecting);
      appendOptionalDouble_(lib_intensity_b, row.transition_library_intensity);
      decoy_b.Append(row.transition_decoy);
      appendOptionalInt64_(feature_id_b, row.feature_id);
      for (Size i = 0; i < table.feature_transition_column_names.size(); ++i)
      {
        if (i < row.feature_transition_values.size() && row.feature_transition_values[i].has_value())
        {
          ft_builders[i].Append(*row.feature_transition_values[i]);
        }
        else
        {
          ft_builders[i].AppendNull();
        }
      }
      appendOptionalDouble_(score_b, row.score_transition_score);
      if (row.score_transition_rank.has_value()) score_rank_b.Append(*row.score_transition_rank);
      else score_rank_b.AppendNull();
      appendOptionalDouble_(pvalue_b, row.score_transition_pvalue);
      appendOptionalDouble_(qvalue_b, row.score_transition_qvalue);
      appendOptionalDouble_(pep_b, row.score_transition_pep);
    }

    auto pushArray = [&](const char* name, const std::shared_ptr<arrow::DataType>& type, std::shared_ptr<arrow::Array> array)
    {
      fields.push_back(arrow::field(name, type));
      arrays.push_back(std::move(array));
    };

    pushArray("RUN_ID", arrow::int64(), run_id_b.Finish().ValueOrDie());
    pushArray("IPF_PEPTIDE_ID", arrow::int64(), ipf_peptide_id_b.Finish().ValueOrDie());
    pushArray("PRECURSOR_ID", arrow::int64(), precursor_id_b.Finish().ValueOrDie());
    pushArray("TRANSITION_ID", arrow::int64(), transition_id_b.Finish().ValueOrDie());
    pushArray("TRANSITION_TRAML_ID", arrow::utf8(), traml_b.Finish().ValueOrDie());
    pushArray("PRODUCT_MZ", arrow::float64(), product_mz_b.Finish().ValueOrDie());
    pushArray("TRANSITION_CHARGE", arrow::int64(), charge_b.Finish().ValueOrDie());
    pushArray("TRANSITION_TYPE", arrow::utf8(), type_b.Finish().ValueOrDie());
    pushArray("TRANSITION_ORDINAL", arrow::int64(), ordinal_b.Finish().ValueOrDie());
    pushArray("ANNOTATION", arrow::utf8(), annotation_b.Finish().ValueOrDie());
    pushArray("TRANSITION_DETECTING", arrow::boolean(), detecting_b.Finish().ValueOrDie());
    pushArray("TRANSITION_LIBRARY_INTENSITY", arrow::float64(), lib_intensity_b.Finish().ValueOrDie());
    pushArray("TRANSITION_DECOY", arrow::boolean(), decoy_b.Finish().ValueOrDie());
    pushArray("FEATURE_ID", arrow::int64(), feature_id_b.Finish().ValueOrDie());
    for (Size i = 0; i < table.feature_transition_column_names.size(); ++i)
    {
      pushArray(("FEATURE_TRANSITION_" + table.feature_transition_column_names[i]).c_str(), arrow::float64(), ft_builders[i].Finish().ValueOrDie());
    }
    pushArray("SCORE_TRANSITION_SCORE", arrow::float64(), score_b.Finish().ValueOrDie());
    pushArray("SCORE_TRANSITION_RANK", arrow::int64(), score_rank_b.Finish().ValueOrDie());
    pushArray("SCORE_TRANSITION_PVALUE", arrow::float64(), pvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_TRANSITION_QVALUE", arrow::float64(), qvalue_b.Finish().ValueOrDie());
    pushArray("SCORE_TRANSITION_PEP", arrow::float64(), pep_b.Finish().ValueOrDie());

    const auto arrow_table = arrow::Table::Make(arrow::schema(fields), arrays);
    if (!ArrowIOHelpers::writeTableToParquet(arrow_table, filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }
} // namespace OpenMS
