// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathResultsExporter.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>

#include <arrow/api.h>

#include <fstream>

namespace OpenMS
{
  namespace
  {
    String optionalToString_(const std::optional<double>& value)
    {
      return value.has_value() ? String(*value) : String();
    }

    String optionalToString_(const std::optional<Int64>& value)
    {
      return value.has_value() ? String(*value) : String();
    }

    void appendOptionalDouble_(arrow::DoubleBuilder& builder, const std::optional<double>& value)
    {
      if (value.has_value())
      {
        builder.Append(*value);
      }
      else
      {
        builder.AppendNull();
      }
    }

    void appendOptionalInt64_(arrow::Int64Builder& builder, const std::optional<Int64>& value)
    {
      if (value.has_value())
      {
        builder.Append(*value);
      }
      else
      {
        builder.AppendNull();
      }
    }

    void appendString_(arrow::StringBuilder& builder, const String& value)
    {
      if (value.empty())
      {
        builder.AppendNull();
      }
      else
      {
        builder.Append(value.c_str());
      }
    }
  } // namespace

  void OpenSwathResultsExporter::write(const String& filename,
                                       const std::vector<OpenSwathExportRow>& rows,
                                       const OpenSwathResultsExportConfig& config)
  {
    if (config.format == OpenSwathExportFileFormat::TSV)
    {
      std::ofstream os(filename.c_str());
      if (!os)
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }

      os << "run_id\tfilename\trun_name\tfeature_id\tpeptide_id\ttransition_group_id\tprecursor_id\tdecoy\tSequence\tFullPeptideName\tProteinName\tGeneName\tCharge\tmz\tRT\tassay_rt\tdelta_rt\tiRT\tassay_iRT\tdelta_iRT\tIntensity\taggr_prec_Peak_Area\taggr_prec_Peak_Apex\tleftWidth\trightWidth\tEXP_IM\tIM_leftWidth\tIM_rightWidth\tms1_pep\tms2_pep\tprecursor_pep\tipf_pep\tpeak_group_rank\td_score\tm_score\tpep\tms2_m_score\tm_score_peptide_run_specific\tm_score_peptide_experiment_wide\tm_score_peptide_global\tm_score_protein_run_specific\tm_score_protein_experiment_wide\tm_score_protein_global\tm_score_gene_run_specific\tm_score_gene_experiment_wide\tm_score_gene_global\talignment_group_id\talignment_reference_feature_id\talignment_reference_rt\talignment_pep\talignment_qvalue\tfrom_alignment\taggr_Peak_Area\taggr_Peak_Apex\taggr_Fragment_Annotation\tipf_FullUniModPeptideName\tipf_precursor_peakgroup_pep\tipf_peptidoform_pep\tipf_peptidoform_m_score\n";
      for (const auto& row : rows)
      {
        os << row.run_id << '\t'
           << row.filename << '\t'
           << row.run_name << '\t'
           << row.feature_id << '\t'
           << row.peptide_id << '\t'
           << row.transition_group_id << '\t'
           << row.precursor_id << '\t'
           << (row.decoy ? 1 : 0) << '\t'
           << row.sequence << '\t'
           << row.full_peptide_name << '\t'
           << row.protein_name << '\t'
           << row.gene_name << '\t'
           << row.charge << '\t'
           << row.mz << '\t'
           << row.rt << '\t'
           << row.assay_rt << '\t'
           << row.delta_rt << '\t'
           << row.irt << '\t'
           << row.assay_irt << '\t'
           << row.delta_irt << '\t'
           << row.intensity << '\t'
           << optionalToString_(row.aggr_prec_peak_area) << '\t'
           << optionalToString_(row.aggr_prec_peak_apex) << '\t'
           << row.left_width << '\t'
           << row.right_width << '\t'
           << optionalToString_(row.exp_im) << '\t'
           << optionalToString_(row.im_left_width) << '\t'
           << optionalToString_(row.im_right_width) << '\t'
           << optionalToString_(row.ms1_pep) << '\t'
           << optionalToString_(row.ms2_pep) << '\t'
           << optionalToString_(row.precursor_pep) << '\t'
           << optionalToString_(row.ipf_pep) << '\t'
           << row.peak_group_rank << '\t'
           << row.d_score << '\t'
           << row.m_score << '\t'
           << optionalToString_(row.pep) << '\t'
           << optionalToString_(row.ms2_m_score) << '\t'
           << optionalToString_(row.peptide_run_specific_qvalue) << '\t'
           << optionalToString_(row.peptide_experiment_wide_qvalue) << '\t'
           << optionalToString_(row.peptide_global_qvalue) << '\t'
           << optionalToString_(row.protein_run_specific_qvalue) << '\t'
           << optionalToString_(row.protein_experiment_wide_qvalue) << '\t'
           << optionalToString_(row.protein_global_qvalue) << '\t'
           << optionalToString_(row.gene_run_specific_qvalue) << '\t'
           << optionalToString_(row.gene_experiment_wide_qvalue) << '\t'
           << optionalToString_(row.gene_global_qvalue) << '\t'
           << optionalToString_(row.alignment_group_id) << '\t'
           << optionalToString_(row.alignment_reference_feature_id) << '\t'
           << optionalToString_(row.alignment_reference_rt) << '\t'
           << optionalToString_(row.alignment_pep) << '\t'
           << optionalToString_(row.alignment_qvalue) << '\t'
           << (row.from_alignment ? 1 : 0) << '\t'
           << row.aggr_peak_area << '\t'
           << row.aggr_peak_apex << '\t'
           << row.aggr_fragment_annotation << '\t'
           << row.ipf_full_peptide_name << '\t'
           << optionalToString_(row.ipf_precursor_peakgroup_pep) << '\t'
           << optionalToString_(row.ipf_peptidoform_pep) << '\t'
           << optionalToString_(row.ipf_peptidoform_m_score) << '\n';
      }
      return;
    }

    arrow::Int64Builder run_id_b;
    arrow::StringBuilder filename_b, run_name_b, tg_b, seq_b, fullpep_b, prot_b, gene_b;
    arrow::Int64Builder feature_id_b, peptide_id_b, precursor_id_b, charge_b;
    arrow::BooleanBuilder decoy_b, from_alignment_b;
    arrow::DoubleBuilder mz_b, rt_b, assay_rt_b, delta_rt_b, irt_b, assay_irt_b, delta_irt_b, intensity_b;
    arrow::DoubleBuilder left_width_b, right_width_b;
    arrow::DoubleBuilder aggr_prec_area_b, aggr_prec_apex_b, exp_im_b, im_left_b, im_right_b;
    arrow::DoubleBuilder ms1_pep_b, ms2_pep_b, precursor_pep_b, ipf_pep_b, d_score_b, m_score_b, pep_b, ms2_m_score_b;
    arrow::Int64Builder peak_group_rank_b, alignment_group_id_b, alignment_reference_feature_id_b;
    arrow::DoubleBuilder pep_rs_b, pep_ew_b, pep_global_b, prot_rs_b, prot_ew_b, prot_global_b, gene_rs_b, gene_ew_b, gene_global_b;
    arrow::DoubleBuilder alignment_reference_rt_b, alignment_pep_b, alignment_qvalue_b;
    arrow::StringBuilder aggr_peak_area_b, aggr_peak_apex_b, aggr_fragment_b, ipf_fullpep_b;
    arrow::DoubleBuilder ipf_precursor_peakgroup_pep_b, ipf_peptidoform_pep_b, ipf_peptidoform_m_score_b;

    for (const auto& row : rows)
    {
      run_id_b.Append(row.run_id);
      appendString_(filename_b, row.filename);
      appendString_(run_name_b, row.run_name);
      feature_id_b.Append(row.feature_id);
      peptide_id_b.Append(row.peptide_id);
      appendString_(tg_b, row.transition_group_id);
      precursor_id_b.Append(row.precursor_id);
      decoy_b.Append(row.decoy);
      appendString_(seq_b, row.sequence);
      appendString_(fullpep_b, row.full_peptide_name);
      appendString_(prot_b, row.protein_name);
      appendString_(gene_b, row.gene_name);
      charge_b.Append(row.charge);
      mz_b.Append(row.mz);
      rt_b.Append(row.rt);
      assay_rt_b.Append(row.assay_rt);
      delta_rt_b.Append(row.delta_rt);
      irt_b.Append(row.irt);
      assay_irt_b.Append(row.assay_irt);
      delta_irt_b.Append(row.delta_irt);
      intensity_b.Append(row.intensity);
      appendOptionalDouble_(aggr_prec_area_b, row.aggr_prec_peak_area);
      appendOptionalDouble_(aggr_prec_apex_b, row.aggr_prec_peak_apex);
      left_width_b.Append(row.left_width);
      right_width_b.Append(row.right_width);
      appendOptionalDouble_(exp_im_b, row.exp_im);
      appendOptionalDouble_(im_left_b, row.im_left_width);
      appendOptionalDouble_(im_right_b, row.im_right_width);
      appendOptionalDouble_(ms1_pep_b, row.ms1_pep);
      appendOptionalDouble_(ms2_pep_b, row.ms2_pep);
      appendOptionalDouble_(precursor_pep_b, row.precursor_pep);
      appendOptionalDouble_(ipf_pep_b, row.ipf_pep);
      peak_group_rank_b.Append(row.peak_group_rank);
      d_score_b.Append(row.d_score);
      m_score_b.Append(row.m_score);
      appendOptionalDouble_(pep_b, row.pep);
      appendOptionalDouble_(ms2_m_score_b, row.ms2_m_score);
      appendOptionalDouble_(pep_rs_b, row.peptide_run_specific_qvalue);
      appendOptionalDouble_(pep_ew_b, row.peptide_experiment_wide_qvalue);
      appendOptionalDouble_(pep_global_b, row.peptide_global_qvalue);
      appendOptionalDouble_(prot_rs_b, row.protein_run_specific_qvalue);
      appendOptionalDouble_(prot_ew_b, row.protein_experiment_wide_qvalue);
      appendOptionalDouble_(prot_global_b, row.protein_global_qvalue);
      appendOptionalDouble_(gene_rs_b, row.gene_run_specific_qvalue);
      appendOptionalDouble_(gene_ew_b, row.gene_experiment_wide_qvalue);
      appendOptionalDouble_(gene_global_b, row.gene_global_qvalue);
      appendOptionalInt64_(alignment_group_id_b, row.alignment_group_id);
      appendOptionalInt64_(alignment_reference_feature_id_b, row.alignment_reference_feature_id);
      appendOptionalDouble_(alignment_reference_rt_b, row.alignment_reference_rt);
      appendOptionalDouble_(alignment_pep_b, row.alignment_pep);
      appendOptionalDouble_(alignment_qvalue_b, row.alignment_qvalue);
      from_alignment_b.Append(row.from_alignment);
      appendString_(aggr_peak_area_b, row.aggr_peak_area);
      appendString_(aggr_peak_apex_b, row.aggr_peak_apex);
      appendString_(aggr_fragment_b, row.aggr_fragment_annotation);
      appendString_(ipf_fullpep_b, row.ipf_full_peptide_name);
      appendOptionalDouble_(ipf_precursor_peakgroup_pep_b, row.ipf_precursor_peakgroup_pep);
      appendOptionalDouble_(ipf_peptidoform_pep_b, row.ipf_peptidoform_pep);
      appendOptionalDouble_(ipf_peptidoform_m_score_b, row.ipf_peptidoform_m_score);
    }

    auto schema = arrow::schema({
      arrow::field("run_id", arrow::int64()),
      arrow::field("filename", arrow::utf8()),
      arrow::field("run_name", arrow::utf8()),
      arrow::field("feature_id", arrow::int64()),
      arrow::field("peptide_id", arrow::int64()),
      arrow::field("transition_group_id", arrow::utf8()),
      arrow::field("precursor_id", arrow::int64()),
      arrow::field("decoy", arrow::boolean()),
      arrow::field("Sequence", arrow::utf8()),
      arrow::field("FullPeptideName", arrow::utf8()),
      arrow::field("ProteinName", arrow::utf8()),
      arrow::field("GeneName", arrow::utf8()),
      arrow::field("Charge", arrow::int64()),
      arrow::field("mz", arrow::float64()),
      arrow::field("RT", arrow::float64()),
      arrow::field("assay_rt", arrow::float64()),
      arrow::field("delta_rt", arrow::float64()),
      arrow::field("iRT", arrow::float64()),
      arrow::field("assay_iRT", arrow::float64()),
      arrow::field("delta_iRT", arrow::float64()),
      arrow::field("Intensity", arrow::float64()),
      arrow::field("aggr_prec_Peak_Area", arrow::float64()),
      arrow::field("aggr_prec_Peak_Apex", arrow::float64()),
      arrow::field("leftWidth", arrow::float64()),
      arrow::field("rightWidth", arrow::float64()),
      arrow::field("EXP_IM", arrow::float64()),
      arrow::field("IM_leftWidth", arrow::float64()),
      arrow::field("IM_rightWidth", arrow::float64()),
      arrow::field("ms1_pep", arrow::float64()),
      arrow::field("ms2_pep", arrow::float64()),
      arrow::field("precursor_pep", arrow::float64()),
      arrow::field("ipf_pep", arrow::float64()),
      arrow::field("peak_group_rank", arrow::int64()),
      arrow::field("d_score", arrow::float64()),
      arrow::field("m_score", arrow::float64()),
      arrow::field("pep", arrow::float64()),
      arrow::field("ms2_m_score", arrow::float64()),
      arrow::field("m_score_peptide_run_specific", arrow::float64()),
      arrow::field("m_score_peptide_experiment_wide", arrow::float64()),
      arrow::field("m_score_peptide_global", arrow::float64()),
      arrow::field("m_score_protein_run_specific", arrow::float64()),
      arrow::field("m_score_protein_experiment_wide", arrow::float64()),
      arrow::field("m_score_protein_global", arrow::float64()),
      arrow::field("m_score_gene_run_specific", arrow::float64()),
      arrow::field("m_score_gene_experiment_wide", arrow::float64()),
      arrow::field("m_score_gene_global", arrow::float64()),
      arrow::field("alignment_group_id", arrow::int64()),
      arrow::field("alignment_reference_feature_id", arrow::int64()),
      arrow::field("alignment_reference_rt", arrow::float64()),
      arrow::field("alignment_pep", arrow::float64()),
      arrow::field("alignment_qvalue", arrow::float64()),
      arrow::field("from_alignment", arrow::boolean()),
      arrow::field("aggr_Peak_Area", arrow::utf8()),
      arrow::field("aggr_Peak_Apex", arrow::utf8()),
      arrow::field("aggr_Fragment_Annotation", arrow::utf8()),
      arrow::field("ipf_FullUniModPeptideName", arrow::utf8()),
      arrow::field("ipf_precursor_peakgroup_pep", arrow::float64()),
      arrow::field("ipf_peptidoform_pep", arrow::float64()),
      arrow::field("ipf_peptidoform_m_score", arrow::float64())
    });

    auto table = arrow::Table::Make(schema, {
      run_id_b.Finish().ValueOrDie(),
      filename_b.Finish().ValueOrDie(),
      run_name_b.Finish().ValueOrDie(),
      feature_id_b.Finish().ValueOrDie(),
      peptide_id_b.Finish().ValueOrDie(),
      tg_b.Finish().ValueOrDie(),
      precursor_id_b.Finish().ValueOrDie(),
      decoy_b.Finish().ValueOrDie(),
      seq_b.Finish().ValueOrDie(),
      fullpep_b.Finish().ValueOrDie(),
      prot_b.Finish().ValueOrDie(),
      gene_b.Finish().ValueOrDie(),
      charge_b.Finish().ValueOrDie(),
      mz_b.Finish().ValueOrDie(),
      rt_b.Finish().ValueOrDie(),
      assay_rt_b.Finish().ValueOrDie(),
      delta_rt_b.Finish().ValueOrDie(),
      irt_b.Finish().ValueOrDie(),
      assay_irt_b.Finish().ValueOrDie(),
      delta_irt_b.Finish().ValueOrDie(),
      intensity_b.Finish().ValueOrDie(),
      aggr_prec_area_b.Finish().ValueOrDie(),
      aggr_prec_apex_b.Finish().ValueOrDie(),
      left_width_b.Finish().ValueOrDie(),
      right_width_b.Finish().ValueOrDie(),
      exp_im_b.Finish().ValueOrDie(),
      im_left_b.Finish().ValueOrDie(),
      im_right_b.Finish().ValueOrDie(),
      ms1_pep_b.Finish().ValueOrDie(),
      ms2_pep_b.Finish().ValueOrDie(),
      precursor_pep_b.Finish().ValueOrDie(),
      ipf_pep_b.Finish().ValueOrDie(),
      peak_group_rank_b.Finish().ValueOrDie(),
      d_score_b.Finish().ValueOrDie(),
      m_score_b.Finish().ValueOrDie(),
      pep_b.Finish().ValueOrDie(),
      ms2_m_score_b.Finish().ValueOrDie(),
      pep_rs_b.Finish().ValueOrDie(),
      pep_ew_b.Finish().ValueOrDie(),
      pep_global_b.Finish().ValueOrDie(),
      prot_rs_b.Finish().ValueOrDie(),
      prot_ew_b.Finish().ValueOrDie(),
      prot_global_b.Finish().ValueOrDie(),
      gene_rs_b.Finish().ValueOrDie(),
      gene_ew_b.Finish().ValueOrDie(),
      gene_global_b.Finish().ValueOrDie(),
      alignment_group_id_b.Finish().ValueOrDie(),
      alignment_reference_feature_id_b.Finish().ValueOrDie(),
      alignment_reference_rt_b.Finish().ValueOrDie(),
      alignment_pep_b.Finish().ValueOrDie(),
      alignment_qvalue_b.Finish().ValueOrDie(),
      from_alignment_b.Finish().ValueOrDie(),
      aggr_peak_area_b.Finish().ValueOrDie(),
      aggr_peak_apex_b.Finish().ValueOrDie(),
      aggr_fragment_b.Finish().ValueOrDie(),
      ipf_fullpep_b.Finish().ValueOrDie(),
      ipf_precursor_peakgroup_pep_b.Finish().ValueOrDie(),
      ipf_peptidoform_pep_b.Finish().ValueOrDie(),
      ipf_peptidoform_m_score_b.Finish().ValueOrDie()
    });

    if (!ArrowIOHelpers::writeTableToParquet(table, filename))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
  }
} // namespace OpenMS
