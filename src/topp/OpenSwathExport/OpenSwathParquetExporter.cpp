// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenSwathParquetExporter.h"

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/FORMAT/ParquetFile.h>

#include <arrow/api.h>

namespace OpenMS
{
  namespace
  {
    void appendOptionalDouble_(arrow::DoubleBuilder& builder, const std::optional<double>& value, const char* column)
    {
      if (value.has_value()) ParquetFile::appendOrThrow(builder.Append(*value), column);
      else ParquetFile::appendOrThrow(builder.AppendNull(), column);
    }

    void appendOptionalInt64_(arrow::Int64Builder& builder, const std::optional<Int64>& value, const char* column)
    {
      if (value.has_value()) ParquetFile::appendOrThrow(builder.Append(*value), column);
      else ParquetFile::appendOrThrow(builder.AppendNull(), column);
    }

    void appendOptionalBool_(arrow::BooleanBuilder& builder, const std::optional<bool>& value, const char* column)
    {
      if (value.has_value()) ParquetFile::appendOrThrow(builder.Append(*value), column);
      else ParquetFile::appendOrThrow(builder.AppendNull(), column);
    }

    void appendString_(arrow::StringBuilder& builder, const std::string& value, const char* column)
    {
      if (value.empty()) ParquetFile::appendOrThrow(builder.AppendNull(), column);
      else ParquetFile::appendOrThrow(builder.Append(value.c_str()), column);
    }
  }

  void OpenSwathParquetExporter::writeFeatureScores(const std::string& filename,
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
      ParquetFile::appendOrThrow(protein_id_b.Append(row.protein_id), "PROTEIN_ID");
      ParquetFile::appendOrThrow(peptide_id_b.Append(row.peptide_id), "PEPTIDE_ID");
      appendOptionalInt64_(ipf_peptide_id_b, row.ipf_peptide_id, "IPF_PEPTIDE_ID");
      ParquetFile::appendOrThrow(precursor_id_b.Append(row.precursor_id), "PRECURSOR_ID");
      appendString_(prot_b, row.protein_accession, "PROTEIN_ACCESSION");
      appendString_(seq_b, row.unmodified_sequence, "UNMODIFIED_SEQUENCE");
      appendString_(fullpep_b, row.modified_sequence, "MODIFIED_SEQUENCE");
      appendString_(traml_b, row.precursor_traml_id, "PRECURSOR_TRAML_ID");
      appendString_(group_label_b, row.precursor_group_label, "PRECURSOR_GROUP_LABEL");
      ParquetFile::appendOrThrow(precursor_mz_b.Append(row.precursor_mz), "PRECURSOR_MZ");
      ParquetFile::appendOrThrow(precursor_charge_b.Append(row.precursor_charge), "PRECURSOR_CHARGE");
      appendOptionalDouble_(lib_intensity_b, row.precursor_library_intensity, "PRECURSOR_LIBRARY_INTENSITY");
      appendOptionalDouble_(lib_rt_b, row.precursor_library_rt, "PRECURSOR_LIBRARY_RT");
      appendOptionalDouble_(lib_dt_b, row.precursor_library_drift_time, "PRECURSOR_LIBRARY_DRIFT_TIME");
      appendOptionalInt64_(gene_id_b, row.gene_id, "GENE_ID");
      appendString_(gene_name_b, row.gene_name, "GENE_NAME");
      appendOptionalBool_(gene_decoy_b, row.gene_decoy, "GENE_DECOY");
      ParquetFile::appendOrThrow(protein_decoy_b.Append(row.protein_decoy), "PROTEIN_DECOY");
      ParquetFile::appendOrThrow(peptide_decoy_b.Append(row.peptide_decoy), "PEPTIDE_DECOY");
      ParquetFile::appendOrThrow(precursor_decoy_b.Append(row.precursor_decoy), "PRECURSOR_DECOY");
      ParquetFile::appendOrThrow(run_id_b.Append(row.run_id), "RUN_ID");
      appendString_(file_b, row.filename, "FILENAME");
      ParquetFile::appendOrThrow(feature_id_b.Append(row.feature_id), "FEATURE_ID");
      ParquetFile::appendOrThrow(exp_rt_b.Append(row.exp_rt), "EXP_RT");
      appendOptionalDouble_(exp_im_b, row.exp_im, "EXP_IM");
      ParquetFile::appendOrThrow(norm_rt_b.Append(row.norm_rt), "NORM_RT");
      ParquetFile::appendOrThrow(delta_rt_b.Append(row.delta_rt), "DELTA_RT");
      ParquetFile::appendOrThrow(left_b.Append(row.left_width), "LEFT_WIDTH");
      ParquetFile::appendOrThrow(right_b.Append(row.right_width), "RIGHT_WIDTH");
      appendOptionalDouble_(im_left_b, row.im_left_width, "IM_leftWidth");
      appendOptionalDouble_(im_right_b, row.im_right_width, "IM_rightWidth");
      for (Size i = 0; i < table.feature_ms1_column_names.size(); ++i)
      {
        const std::string column = "FEATURE_MS1_" + table.feature_ms1_column_names[i];
        if (i < row.feature_ms1_values.size()) ParquetFile::appendOrThrow(ms1_builders[i].Append(row.feature_ms1_values[i]), column.c_str());
        else ParquetFile::appendOrThrow(ms1_builders[i].AppendNull(), column.c_str());
      }
      for (Size i = 0; i < table.feature_ms2_column_names.size(); ++i)
      {
        const std::string column = "FEATURE_MS2_" + table.feature_ms2_column_names[i];
        if (i < row.feature_ms2_values.size()) ParquetFile::appendOrThrow(ms2_builders[i].Append(row.feature_ms2_values[i]), column.c_str());
        else ParquetFile::appendOrThrow(ms2_builders[i].AppendNull(), column.c_str());
      }
      appendOptionalDouble_(score_ms1_score_b, row.score_ms1_score, "SCORE_MS1_SCORE");
      if (row.score_ms1_rank.has_value()) ParquetFile::appendOrThrow(score_ms1_rank_b.Append(*row.score_ms1_rank), "SCORE_MS1_RANK");
      else ParquetFile::appendOrThrow(score_ms1_rank_b.AppendNull(), "SCORE_MS1_RANK");
      appendOptionalDouble_(score_ms1_pvalue_b, row.score_ms1_pvalue, "SCORE_MS1_PVALUE");
      appendOptionalDouble_(score_ms1_qvalue_b, row.score_ms1_qvalue, "SCORE_MS1_QVALUE");
      appendOptionalDouble_(score_ms1_pep_b, row.score_ms1_pep, "SCORE_MS1_PEP");
      appendOptionalDouble_(score_ms2_score_b, row.score_ms2_score, "SCORE_MS2_SCORE");
      if (row.score_ms2_peak_group_rank.has_value()) ParquetFile::appendOrThrow(score_ms2_rank_b.Append(*row.score_ms2_peak_group_rank), "SCORE_MS2_PEAK_GROUP_RANK");
      else ParquetFile::appendOrThrow(score_ms2_rank_b.AppendNull(), "SCORE_MS2_PEAK_GROUP_RANK");
      appendOptionalDouble_(score_ms2_pvalue_b, row.score_ms2_pvalue, "SCORE_MS2_PVALUE");
      appendOptionalDouble_(score_ms2_qvalue_b, row.score_ms2_qvalue, "SCORE_MS2_QVALUE");
      appendOptionalDouble_(score_ms2_pep_b, row.score_ms2_pep, "SCORE_MS2_PEP");
      appendOptionalDouble_(score_ipf_prec_b, row.score_ipf_precursor_peakgroup_pep, "SCORE_IPF_PRECURSOR_PEAKGROUP_PEP");
      appendOptionalDouble_(score_ipf_pep_b, row.score_ipf_pep, "SCORE_IPF_PEP");
      appendOptionalDouble_(score_ipf_qvalue_b, row.score_ipf_qvalue, "SCORE_IPF_QVALUE");
      appendOptionalDouble_(pep_global_score_b, row.score_peptide_global_score, "SCORE_PEPTIDE_GLOBAL_SCORE");
      appendOptionalDouble_(pep_global_pvalue_b, row.score_peptide_global_pvalue, "SCORE_PEPTIDE_GLOBAL_PVALUE");
      appendOptionalDouble_(pep_global_qvalue_b, row.score_peptide_global_qvalue, "SCORE_PEPTIDE_GLOBAL_QVALUE");
      appendOptionalDouble_(pep_global_pep_b, row.score_peptide_global_pep, "SCORE_PEPTIDE_GLOBAL_PEP");
      appendOptionalDouble_(pep_ew_score_b, row.score_peptide_experiment_wide_score, "SCORE_PEPTIDE_EXPERIMENT_WIDE_SCORE");
      appendOptionalDouble_(pep_ew_pvalue_b, row.score_peptide_experiment_wide_pvalue, "SCORE_PEPTIDE_EXPERIMENT_WIDE_PVALUE");
      appendOptionalDouble_(pep_ew_qvalue_b, row.score_peptide_experiment_wide_qvalue, "SCORE_PEPTIDE_EXPERIMENT_WIDE_QVALUE");
      appendOptionalDouble_(pep_ew_pep_b, row.score_peptide_experiment_wide_pep, "SCORE_PEPTIDE_EXPERIMENT_WIDE_PEP");
      appendOptionalDouble_(pep_rs_score_b, row.score_peptide_run_specific_score, "SCORE_PEPTIDE_RUN_SPECIFIC_SCORE");
      appendOptionalDouble_(pep_rs_pvalue_b, row.score_peptide_run_specific_pvalue, "SCORE_PEPTIDE_RUN_SPECIFIC_PVALUE");
      appendOptionalDouble_(pep_rs_qvalue_b, row.score_peptide_run_specific_qvalue, "SCORE_PEPTIDE_RUN_SPECIFIC_QVALUE");
      appendOptionalDouble_(pep_rs_pep_b, row.score_peptide_run_specific_pep, "SCORE_PEPTIDE_RUN_SPECIFIC_PEP");
      appendOptionalDouble_(prot_global_score_b, row.score_protein_global_score, "SCORE_PROTEIN_GLOBAL_SCORE");
      appendOptionalDouble_(prot_global_pvalue_b, row.score_protein_global_pvalue, "SCORE_PROTEIN_GLOBAL_PVALUE");
      appendOptionalDouble_(prot_global_qvalue_b, row.score_protein_global_qvalue, "SCORE_PROTEIN_GLOBAL_QVALUE");
      appendOptionalDouble_(prot_global_pep_b, row.score_protein_global_pep, "SCORE_PROTEIN_GLOBAL_PEP");
      appendOptionalDouble_(prot_ew_score_b, row.score_protein_experiment_wide_score, "SCORE_PROTEIN_EXPERIMENT_WIDE_SCORE");
      appendOptionalDouble_(prot_ew_pvalue_b, row.score_protein_experiment_wide_pvalue, "SCORE_PROTEIN_EXPERIMENT_WIDE_PVALUE");
      appendOptionalDouble_(prot_ew_qvalue_b, row.score_protein_experiment_wide_qvalue, "SCORE_PROTEIN_EXPERIMENT_WIDE_QVALUE");
      appendOptionalDouble_(prot_ew_pep_b, row.score_protein_experiment_wide_pep, "SCORE_PROTEIN_EXPERIMENT_WIDE_PEP");
      appendOptionalDouble_(prot_rs_score_b, row.score_protein_run_specific_score, "SCORE_PROTEIN_RUN_SPECIFIC_SCORE");
      appendOptionalDouble_(prot_rs_pvalue_b, row.score_protein_run_specific_pvalue, "SCORE_PROTEIN_RUN_SPECIFIC_PVALUE");
      appendOptionalDouble_(prot_rs_qvalue_b, row.score_protein_run_specific_qvalue, "SCORE_PROTEIN_RUN_SPECIFIC_QVALUE");
      appendOptionalDouble_(prot_rs_pep_b, row.score_protein_run_specific_pep, "SCORE_PROTEIN_RUN_SPECIFIC_PEP");
      appendOptionalDouble_(gene_global_score_b, row.score_gene_global_score, "SCORE_GENE_GLOBAL_SCORE");
      appendOptionalDouble_(gene_global_pvalue_b, row.score_gene_global_pvalue, "SCORE_GENE_GLOBAL_PVALUE");
      appendOptionalDouble_(gene_global_qvalue_b, row.score_gene_global_qvalue, "SCORE_GENE_GLOBAL_QVALUE");
      appendOptionalDouble_(gene_global_pep_b, row.score_gene_global_pep, "SCORE_GENE_GLOBAL_PEP");
      appendOptionalDouble_(gene_ew_score_b, row.score_gene_experiment_wide_score, "SCORE_GENE_EXPERIMENT_WIDE_SCORE");
      appendOptionalDouble_(gene_ew_pvalue_b, row.score_gene_experiment_wide_pvalue, "SCORE_GENE_EXPERIMENT_WIDE_PVALUE");
      appendOptionalDouble_(gene_ew_qvalue_b, row.score_gene_experiment_wide_qvalue, "SCORE_GENE_EXPERIMENT_WIDE_QVALUE");
      appendOptionalDouble_(gene_ew_pep_b, row.score_gene_experiment_wide_pep, "SCORE_GENE_EXPERIMENT_WIDE_PEP");
      appendOptionalDouble_(gene_rs_score_b, row.score_gene_run_specific_score, "SCORE_GENE_RUN_SPECIFIC_SCORE");
      appendOptionalDouble_(gene_rs_pvalue_b, row.score_gene_run_specific_pvalue, "SCORE_GENE_RUN_SPECIFIC_PVALUE");
      appendOptionalDouble_(gene_rs_qvalue_b, row.score_gene_run_specific_qvalue, "SCORE_GENE_RUN_SPECIFIC_QVALUE");
      appendOptionalDouble_(gene_rs_pep_b, row.score_gene_run_specific_pep, "SCORE_GENE_RUN_SPECIFIC_PEP");
      appendOptionalInt64_(alignment_group_id_b, row.alignment_group_id, "alignment_group_id");
      appendOptionalInt64_(alignment_reference_id_b, row.alignment_reference_feature_id, "alignment_reference_feature_id");
      appendOptionalDouble_(align_rt_b, row.alignment_reference_rt, "alignment_reference_rt");
      appendOptionalDouble_(align_pep_b, row.alignment_pep, "alignment_pep");
      appendOptionalDouble_(align_qvalue_b, row.alignment_qvalue, "alignment_qvalue");
      ParquetFile::appendOrThrow(from_alignment_b.Append(row.from_alignment), "from_alignment");
    }

    auto pushArray = [&](const char* name, const std::shared_ptr<arrow::DataType>& type, auto& builder)
    {
      fields.push_back(arrow::field(name, type));
      arrays.push_back(ParquetFile::finishArray(builder, name));
    };

    pushArray("PROTEIN_ID", arrow::int64(), protein_id_b);
    pushArray("PEPTIDE_ID", arrow::int64(), peptide_id_b);
    pushArray("IPF_PEPTIDE_ID", arrow::int64(), ipf_peptide_id_b);
    pushArray("PRECURSOR_ID", arrow::int64(), precursor_id_b);
    pushArray("PROTEIN_ACCESSION", arrow::utf8(), prot_b);
    pushArray("UNMODIFIED_SEQUENCE", arrow::utf8(), seq_b);
    pushArray("MODIFIED_SEQUENCE", arrow::utf8(), fullpep_b);
    pushArray("PRECURSOR_TRAML_ID", arrow::utf8(), traml_b);
    pushArray("PRECURSOR_GROUP_LABEL", arrow::utf8(), group_label_b);
    pushArray("PRECURSOR_MZ", arrow::float64(), precursor_mz_b);
    pushArray("PRECURSOR_CHARGE", arrow::int64(), precursor_charge_b);
    pushArray("PRECURSOR_LIBRARY_INTENSITY", arrow::float64(), lib_intensity_b);
    pushArray("PRECURSOR_LIBRARY_RT", arrow::float64(), lib_rt_b);
    pushArray("PRECURSOR_LIBRARY_DRIFT_TIME", arrow::float64(), lib_dt_b);
    pushArray("GENE_ID", arrow::int64(), gene_id_b);
    pushArray("GENE_NAME", arrow::utf8(), gene_name_b);
    pushArray("GENE_DECOY", arrow::boolean(), gene_decoy_b);
    pushArray("PROTEIN_DECOY", arrow::boolean(), protein_decoy_b);
    pushArray("PEPTIDE_DECOY", arrow::boolean(), peptide_decoy_b);
    pushArray("PRECURSOR_DECOY", arrow::boolean(), precursor_decoy_b);
    pushArray("RUN_ID", arrow::int64(), run_id_b);
    pushArray("FILENAME", arrow::utf8(), file_b);
    pushArray("FEATURE_ID", arrow::int64(), feature_id_b);
    pushArray("EXP_RT", arrow::float64(), exp_rt_b);
    pushArray("EXP_IM", arrow::float64(), exp_im_b);
    pushArray("NORM_RT", arrow::float64(), norm_rt_b);
    pushArray("DELTA_RT", arrow::float64(), delta_rt_b);
    pushArray("LEFT_WIDTH", arrow::float64(), left_b);
    pushArray("RIGHT_WIDTH", arrow::float64(), right_b);
    pushArray("IM_leftWidth", arrow::float64(), im_left_b);
    pushArray("IM_rightWidth", arrow::float64(), im_right_b);
    for (Size i = 0; i < table.feature_ms1_column_names.size(); ++i)
    {
      pushArray(("FEATURE_MS1_" + table.feature_ms1_column_names[i]).c_str(), arrow::float64(), ms1_builders[i]);
    }
    for (Size i = 0; i < table.feature_ms2_column_names.size(); ++i)
    {
      pushArray(("FEATURE_MS2_" + table.feature_ms2_column_names[i]).c_str(), arrow::float64(), ms2_builders[i]);
    }
    pushArray("SCORE_MS1_SCORE", arrow::float64(), score_ms1_score_b);
    pushArray("SCORE_MS1_RANK", arrow::int64(), score_ms1_rank_b);
    pushArray("SCORE_MS1_PVALUE", arrow::float64(), score_ms1_pvalue_b);
    pushArray("SCORE_MS1_QVALUE", arrow::float64(), score_ms1_qvalue_b);
    pushArray("SCORE_MS1_PEP", arrow::float64(), score_ms1_pep_b);
    pushArray("SCORE_MS2_SCORE", arrow::float64(), score_ms2_score_b);
    pushArray("SCORE_MS2_PEAK_GROUP_RANK", arrow::int64(), score_ms2_rank_b);
    pushArray("SCORE_MS2_PVALUE", arrow::float64(), score_ms2_pvalue_b);
    pushArray("SCORE_MS2_QVALUE", arrow::float64(), score_ms2_qvalue_b);
    pushArray("SCORE_MS2_PEP", arrow::float64(), score_ms2_pep_b);
    pushArray("SCORE_IPF_PRECURSOR_PEAKGROUP_PEP", arrow::float64(), score_ipf_prec_b);
    pushArray("SCORE_IPF_PEP", arrow::float64(), score_ipf_pep_b);
    pushArray("SCORE_IPF_QVALUE", arrow::float64(), score_ipf_qvalue_b);
    pushArray("SCORE_PEPTIDE_GLOBAL_SCORE", arrow::float64(), pep_global_score_b);
    pushArray("SCORE_PEPTIDE_GLOBAL_PVALUE", arrow::float64(), pep_global_pvalue_b);
    pushArray("SCORE_PEPTIDE_GLOBAL_QVALUE", arrow::float64(), pep_global_qvalue_b);
    pushArray("SCORE_PEPTIDE_GLOBAL_PEP", arrow::float64(), pep_global_pep_b);
    pushArray("SCORE_PEPTIDE_EXPERIMENT_WIDE_SCORE", arrow::float64(), pep_ew_score_b);
    pushArray("SCORE_PEPTIDE_EXPERIMENT_WIDE_PVALUE", arrow::float64(), pep_ew_pvalue_b);
    pushArray("SCORE_PEPTIDE_EXPERIMENT_WIDE_QVALUE", arrow::float64(), pep_ew_qvalue_b);
    pushArray("SCORE_PEPTIDE_EXPERIMENT_WIDE_PEP", arrow::float64(), pep_ew_pep_b);
    pushArray("SCORE_PEPTIDE_RUN_SPECIFIC_SCORE", arrow::float64(), pep_rs_score_b);
    pushArray("SCORE_PEPTIDE_RUN_SPECIFIC_PVALUE", arrow::float64(), pep_rs_pvalue_b);
    pushArray("SCORE_PEPTIDE_RUN_SPECIFIC_QVALUE", arrow::float64(), pep_rs_qvalue_b);
    pushArray("SCORE_PEPTIDE_RUN_SPECIFIC_PEP", arrow::float64(), pep_rs_pep_b);
    pushArray("SCORE_PROTEIN_GLOBAL_SCORE", arrow::float64(), prot_global_score_b);
    pushArray("SCORE_PROTEIN_GLOBAL_PVALUE", arrow::float64(), prot_global_pvalue_b);
    pushArray("SCORE_PROTEIN_GLOBAL_QVALUE", arrow::float64(), prot_global_qvalue_b);
    pushArray("SCORE_PROTEIN_GLOBAL_PEP", arrow::float64(), prot_global_pep_b);
    pushArray("SCORE_PROTEIN_EXPERIMENT_WIDE_SCORE", arrow::float64(), prot_ew_score_b);
    pushArray("SCORE_PROTEIN_EXPERIMENT_WIDE_PVALUE", arrow::float64(), prot_ew_pvalue_b);
    pushArray("SCORE_PROTEIN_EXPERIMENT_WIDE_QVALUE", arrow::float64(), prot_ew_qvalue_b);
    pushArray("SCORE_PROTEIN_EXPERIMENT_WIDE_PEP", arrow::float64(), prot_ew_pep_b);
    pushArray("SCORE_PROTEIN_RUN_SPECIFIC_SCORE", arrow::float64(), prot_rs_score_b);
    pushArray("SCORE_PROTEIN_RUN_SPECIFIC_PVALUE", arrow::float64(), prot_rs_pvalue_b);
    pushArray("SCORE_PROTEIN_RUN_SPECIFIC_QVALUE", arrow::float64(), prot_rs_qvalue_b);
    pushArray("SCORE_PROTEIN_RUN_SPECIFIC_PEP", arrow::float64(), prot_rs_pep_b);
    pushArray("SCORE_GENE_GLOBAL_SCORE", arrow::float64(), gene_global_score_b);
    pushArray("SCORE_GENE_GLOBAL_PVALUE", arrow::float64(), gene_global_pvalue_b);
    pushArray("SCORE_GENE_GLOBAL_QVALUE", arrow::float64(), gene_global_qvalue_b);
    pushArray("SCORE_GENE_GLOBAL_PEP", arrow::float64(), gene_global_pep_b);
    pushArray("SCORE_GENE_EXPERIMENT_WIDE_SCORE", arrow::float64(), gene_ew_score_b);
    pushArray("SCORE_GENE_EXPERIMENT_WIDE_PVALUE", arrow::float64(), gene_ew_pvalue_b);
    pushArray("SCORE_GENE_EXPERIMENT_WIDE_QVALUE", arrow::float64(), gene_ew_qvalue_b);
    pushArray("SCORE_GENE_EXPERIMENT_WIDE_PEP", arrow::float64(), gene_ew_pep_b);
    pushArray("SCORE_GENE_RUN_SPECIFIC_SCORE", arrow::float64(), gene_rs_score_b);
    pushArray("SCORE_GENE_RUN_SPECIFIC_PVALUE", arrow::float64(), gene_rs_pvalue_b);
    pushArray("SCORE_GENE_RUN_SPECIFIC_QVALUE", arrow::float64(), gene_rs_qvalue_b);
    pushArray("SCORE_GENE_RUN_SPECIFIC_PEP", arrow::float64(), gene_rs_pep_b);
    pushArray("alignment_group_id", arrow::int64(), alignment_group_id_b);
    pushArray("alignment_reference_feature_id", arrow::int64(), alignment_reference_id_b);
    pushArray("alignment_reference_rt", arrow::float64(), align_rt_b);
    pushArray("alignment_pep", arrow::float64(), align_pep_b);
    pushArray("alignment_qvalue", arrow::float64(), align_qvalue_b);
    pushArray("from_alignment", arrow::boolean(), from_alignment_b);

    const auto arrow_table = arrow::Table::Make(arrow::schema(fields), arrays);
    if (!ArrowIOHelpers::writeTableToParquet(arrow_table, filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }

  void OpenSwathParquetExporter::writeTransitionScores(const std::string& filename,
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
      appendOptionalInt64_(run_id_b, row.run_id, "RUN_ID");
      appendOptionalInt64_(ipf_peptide_id_b, row.ipf_peptide_id, "IPF_PEPTIDE_ID");
      ParquetFile::appendOrThrow(precursor_id_b.Append(row.precursor_id), "PRECURSOR_ID");
      ParquetFile::appendOrThrow(transition_id_b.Append(row.transition_id), "TRANSITION_ID");
      appendString_(traml_b, row.transition_traml_id, "TRANSITION_TRAML_ID");
      ParquetFile::appendOrThrow(product_mz_b.Append(row.product_mz), "PRODUCT_MZ");
      ParquetFile::appendOrThrow(charge_b.Append(row.transition_charge), "TRANSITION_CHARGE");
      appendString_(type_b, row.transition_type, "TRANSITION_TYPE");
      ParquetFile::appendOrThrow(ordinal_b.Append(row.transition_ordinal), "TRANSITION_ORDINAL");
      appendString_(annotation_b, row.annotation, "ANNOTATION");
      ParquetFile::appendOrThrow(detecting_b.Append(row.transition_detecting), "TRANSITION_DETECTING");
      appendOptionalDouble_(lib_intensity_b, row.transition_library_intensity, "TRANSITION_LIBRARY_INTENSITY");
      ParquetFile::appendOrThrow(decoy_b.Append(row.transition_decoy), "TRANSITION_DECOY");
      appendOptionalInt64_(feature_id_b, row.feature_id, "FEATURE_ID");
      for (Size i = 0; i < table.feature_transition_column_names.size(); ++i)
      {
        const std::string column = "FEATURE_TRANSITION_" + table.feature_transition_column_names[i];
        if (i < row.feature_transition_values.size() && row.feature_transition_values[i].has_value())
        {
          ParquetFile::appendOrThrow(ft_builders[i].Append(*row.feature_transition_values[i]), column.c_str());
        }
        else
        {
          ParquetFile::appendOrThrow(ft_builders[i].AppendNull(), column.c_str());
        }
      }
      appendOptionalDouble_(score_b, row.score_transition_score, "SCORE_TRANSITION_SCORE");
      if (row.score_transition_rank.has_value()) ParquetFile::appendOrThrow(score_rank_b.Append(*row.score_transition_rank), "SCORE_TRANSITION_RANK");
      else ParquetFile::appendOrThrow(score_rank_b.AppendNull(), "SCORE_TRANSITION_RANK");
      appendOptionalDouble_(pvalue_b, row.score_transition_pvalue, "SCORE_TRANSITION_PVALUE");
      appendOptionalDouble_(qvalue_b, row.score_transition_qvalue, "SCORE_TRANSITION_QVALUE");
      appendOptionalDouble_(pep_b, row.score_transition_pep, "SCORE_TRANSITION_PEP");
    }

    auto pushArray = [&](const char* name, const std::shared_ptr<arrow::DataType>& type, auto& builder)
    {
      fields.push_back(arrow::field(name, type));
      arrays.push_back(ParquetFile::finishArray(builder, name));
    };

    pushArray("RUN_ID", arrow::int64(), run_id_b);
    pushArray("IPF_PEPTIDE_ID", arrow::int64(), ipf_peptide_id_b);
    pushArray("PRECURSOR_ID", arrow::int64(), precursor_id_b);
    pushArray("TRANSITION_ID", arrow::int64(), transition_id_b);
    pushArray("TRANSITION_TRAML_ID", arrow::utf8(), traml_b);
    pushArray("PRODUCT_MZ", arrow::float64(), product_mz_b);
    pushArray("TRANSITION_CHARGE", arrow::int64(), charge_b);
    pushArray("TRANSITION_TYPE", arrow::utf8(), type_b);
    pushArray("TRANSITION_ORDINAL", arrow::int64(), ordinal_b);
    pushArray("ANNOTATION", arrow::utf8(), annotation_b);
    pushArray("TRANSITION_DETECTING", arrow::boolean(), detecting_b);
    pushArray("TRANSITION_LIBRARY_INTENSITY", arrow::float64(), lib_intensity_b);
    pushArray("TRANSITION_DECOY", arrow::boolean(), decoy_b);
    pushArray("FEATURE_ID", arrow::int64(), feature_id_b);
    for (Size i = 0; i < table.feature_transition_column_names.size(); ++i)
    {
      pushArray(("FEATURE_TRANSITION_" + table.feature_transition_column_names[i]).c_str(), arrow::float64(), ft_builders[i]);
    }
    pushArray("SCORE_TRANSITION_SCORE", arrow::float64(), score_b);
    pushArray("SCORE_TRANSITION_RANK", arrow::int64(), score_rank_b);
    pushArray("SCORE_TRANSITION_PVALUE", arrow::float64(), pvalue_b);
    pushArray("SCORE_TRANSITION_QVALUE", arrow::float64(), qvalue_b);
    pushArray("SCORE_TRANSITION_PEP", arrow::float64(), pep_b);

    const auto arrow_table = arrow::Table::Make(arrow::schema(fields), arrays);
    if (!ArrowIOHelpers::writeTableToParquet(arrow_table, filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }
} // namespace OpenMS
