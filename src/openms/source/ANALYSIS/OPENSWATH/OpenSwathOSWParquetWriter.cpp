// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h>

#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/SYSTEM/ExternalProcess.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <fstream>
#include <cmath>
#include <unordered_map>
#include <vector>

#ifdef WITH_PARQUET
#ifdef signals
#undef signals
#endif
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#endif

namespace OpenMS
{
  namespace
  {
    using OpenMS::Size;

#ifdef WITH_PARQUET
    void appendOrThrow_(const arrow::Status& status, const char* column)
    {
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("Failed to append value for ") + column, status.ToString());
      }
    }

    template <typename Builder>
    std::shared_ptr<arrow::Array> finishArray_(Builder& builder, const char* name)
    {
      std::shared_ptr<arrow::Array> array;
      auto status = builder.Finish(&array);
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      String("Failed to finish array for ") + name, status.ToString());
      }
      return array;
    }

    void appendOptionalFloat_(arrow::DoubleBuilder& builder, bool has_value, double value, const char* column)
    {
      if (!has_value || !std::isfinite(value))
      {
        appendOrThrow_(builder.AppendNull(), column);
      }
      else
      {
        appendOrThrow_(builder.Append(value), column);
      }
    }

    void appendOptionalInt_(arrow::Int64Builder& builder, bool has_value, int64_t value, const char* column)
    {
      if (!has_value)
      {
        appendOrThrow_(builder.AppendNull(), column);
      }
      else
      {
        appendOrThrow_(builder.Append(value), column);
      }
    }

    bool extractMetaDouble_(const Feature& feature, const std::string& key, double& value)
    {
      if (!feature.metaValueExists(key)) return false;
      const DataValue& meta = feature.getMetaValue(key);
      if (meta.isEmpty()) return false;
      try
      {
        value = meta.toString().toDouble();
        return true;
      }
      catch (Exception::ConversionError&)
      {
        return false;
      }
    }

    bool extractMetaDouble_(const BaseFeature& feature, const std::string& key, double& value)
    {
      if (!feature.metaValueExists(key)) return false;
      const DataValue& meta = feature.getMetaValue(key);
      if (meta.isEmpty()) return false;
      try
      {
        value = meta.toString().toDouble();
        return true;
      }
      catch (Exception::ConversionError&)
      {
        return false;
      }
    }

    std::vector<String> getSeparateScore_(const Feature& feature, const std::string& score_name)
    {
      std::vector<String> separated_scores;

      if (!feature.getMetaValue(score_name).isEmpty())
      {
        if (feature.getMetaValue(score_name).valueType() == DataValue::STRING_LIST)
        {
          separated_scores = feature.getMetaValue(score_name).toStringList();
        }
        else if (feature.getMetaValue(score_name).valueType() == DataValue::INT_LIST)
        {
          std::vector<int> int_scores = feature.getMetaValue(score_name).toIntList();
          for (int score : int_scores) separated_scores.emplace_back(score);
        }
        else if (feature.getMetaValue(score_name).valueType() == DataValue::DOUBLE_LIST)
        {
          std::vector<double> double_scores = feature.getMetaValue(score_name).toDoubleList();
          for (double score : double_scores) separated_scores.emplace_back(score);
        }
        else
        {
          separated_scores.push_back(feature.getMetaValue(score_name).toString());
        }
      }

      return separated_scores;
    }

    bool parseOptionalInt64_(const String& text, int64_t& value)
    {
      try
      {
        value = text.toInt64();
        return true;
      }
      catch (Exception::ConversionError&)
      {
        return false;
      }
    }

    bool parseOptionalDouble_(const String& text, double& value)
    {
      try
      {
        value = text.toDouble();
        return true;
      }
      catch (Exception::ConversionError&)
      {
        return false;
      }
    }

    void appendOptionalFloatFromList_(arrow::DoubleBuilder& builder,
                                      const std::vector<String>& values,
                                      Size index,
                                      const char* column)
    {
      double value = 0.0;
      if (index < values.size() && parseOptionalDouble_(values[index], value))
      {
        appendOptionalFloat_(builder, true, value, column);
      }
      else
      {
        appendOrThrow_(builder.AppendNull(), column);
      }
    }

    void writeParquetTable_(const std::shared_ptr<arrow::Table>& table, const String& filename)
    {
      auto outfile_result = arrow::io::FileOutputStream::Open(std::string(filename));
      if (!outfile_result.ok())
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
      }
      auto outfile = outfile_result.ValueOrDie();
      auto status = parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1024);
      if (!status.ok())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to write parquet table", status.ToString());
      }
    }

    void zipParquetDirectory_(const String& directory_path, const String& output_zip)
    {
      const String output_zip_abs = File::absolutePath(output_zip);
      if (File::exists(output_zip_abs))
      {
        File::remove(output_zip_abs);
      }

      ExternalProcess zip_process;
      QStringList args;
      args << "-r" << "-q" << output_zip_abs.toQString() << ".";
      auto status = zip_process.run("zip", args, directory_path.toQString(), false);
      if (status != ExternalProcess::RETURNSTATE::SUCCESS)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Failed to zip parquet output", output_zip_abs);
      }
    }

    String writeMetadata_(const String& base_dir, UInt64 run_id, const String& input_filename)
    {
      const String metadata_path = base_dir + "/metadata.json";
      std::ofstream out(metadata_path.c_str(), std::ios::out | std::ios::trunc);
      if (!out.is_open())
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, metadata_path);
      }

      out << "{\n"
          << "  \"openms\": {\n"
          << "    \"schema_version\": 1,\n"
          << "    \"generator\": \"OpenSwathOSWParquetWriter\",\n"
          << "    \"run\": {\"id\": " << run_id << ", \"filename\": \"" << input_filename << "\"}\n"
          << "  }\n"
          << "}\n";
      return metadata_path;
    }

#endif // WITH_PARQUET
  } // namespace

  void OpenSwathOSWParquetWriter::write(const String& output_path,
                                        const OpenSwath::LightTargetedExperiment& assay_library,
                                        const FeatureMap& feature_map,
                                        UInt64 run_id,
                                        const String& input_filename,
                                        bool enable_uis_scoring) const
  {
#ifndef WITH_PARQUET
    throw Exception::MissingFeature(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "OpenMS was built without Parquet support");
#else
    const UInt64 run_id_clean = Internal::SqliteHelper::clearSignBit(run_id);
    const bool zip_output = output_path.hasSuffix(".osw_parquet");
    String base_dir = output_path;
    File::TempDir temp_dir;
    if (zip_output)
    {
      base_dir = temp_dir.getPath() + "/osw_parquet";
    }

    if (File::exists(base_dir))
    {
      File::removeDirRecursively(base_dir);
    }
    File::makeDir(base_dir);

    const String library_tmp_dir = base_dir + "/library_tmp";
    const String library_dir = base_dir + "/library";
    const String runs_dir = base_dir + "/runs";
    File::makeDir(library_tmp_dir);
    File::makeDir(library_dir);
    File::makeDir(runs_dir);

    writeMetadata_(base_dir, run_id_clean, input_filename);

    TransitionParquetFile().convertLightTargetedExperimentToParquet(library_tmp_dir, assay_library);
    File::copyDirRecursively(String(library_tmp_dir + "/library").toQString(), library_dir.toQString());
    File::removeDirRecursively(library_tmp_dir);

    const String run_path = runs_dir + "/run_id=" + String(run_id_clean);
    File::makeDir(run_path);

    std::unordered_map<String, int64_t> compound_to_precursor;
    compound_to_precursor.reserve(assay_library.compounds.size());
    int64_t next_precursor_id = 1;
    for (const auto& compound : assay_library.compounds)
    {
      if (compound_to_precursor.find(compound.id) != compound_to_precursor.end())
      {
        continue;
      }

      int64_t precursor_id = 0;
      try
      {
        precursor_id = String(compound.id).toInt64();
      }
      catch (Exception::ConversionError&)
      {
        precursor_id = next_precursor_id++;
      }

      if (precursor_id >= next_precursor_id)
      {
        next_precursor_id = precursor_id + 1;
      }
      compound_to_precursor[compound.id] = precursor_id;
    }

    std::unordered_map<String, int64_t> transition_to_id;
    transition_to_id.reserve(assay_library.transitions.size());
    int64_t transition_id = 1;
    for (const auto& transition : assay_library.transitions)
    {
      transition_to_id[transition.transition_name] = transition_id++;
    }

    {
      arrow::Int64Builder run_id_builder;
      arrow::StringBuilder filename_builder;
      appendOrThrow_(run_id_builder.Append(run_id_clean), "run_id");
      appendOrThrow_(filename_builder.Append(std::string(input_filename)), "filename");

      auto runs_schema = arrow::schema({
        arrow::field("run_id", arrow::int64()),
        arrow::field("filename", arrow::utf8())
      });
      auto runs_table = arrow::Table::Make(runs_schema, {
        finishArray_(run_id_builder, "run_id"),
        finishArray_(filename_builder, "filename")
      });
      writeParquetTable_(runs_table, runs_dir + "/runs.parquet");
    }

    arrow::Int64Builder feature_id_builder;
    arrow::Int64Builder feature_run_id_builder;
    arrow::Int64Builder precursor_id_builder;
    arrow::DoubleBuilder exp_rt_builder;
    arrow::DoubleBuilder exp_im_builder;
    arrow::DoubleBuilder norm_rt_builder;
    arrow::DoubleBuilder delta_rt_builder;
    arrow::DoubleBuilder left_width_builder;
    arrow::DoubleBuilder right_width_builder;
    arrow::DoubleBuilder exp_im_left_builder;
    arrow::DoubleBuilder exp_im_right_builder;

    arrow::Int64Builder ms1_feature_id_builder;
    arrow::Int64Builder ms1_run_id_builder;
    arrow::DoubleBuilder ms1_area_builder;
    arrow::DoubleBuilder ms1_apex_builder;
    arrow::DoubleBuilder ms1_exp_im_builder;
    arrow::DoubleBuilder ms1_delta_im_builder;
    arrow::DoubleBuilder ms1_var_massdev_builder;
    arrow::DoubleBuilder ms1_var_im_ms1_delta_builder;
    arrow::DoubleBuilder ms1_var_mi_builder;
    arrow::DoubleBuilder ms1_var_mi_contrast_builder;
    arrow::DoubleBuilder ms1_var_mi_combined_builder;
    arrow::DoubleBuilder ms1_var_iso_corr_builder;
    arrow::DoubleBuilder ms1_var_iso_overlap_builder;
    arrow::DoubleBuilder ms1_var_xcorr_coelution_builder;
    arrow::DoubleBuilder ms1_var_xcorr_coelution_contrast_builder;
    arrow::DoubleBuilder ms1_var_xcorr_coelution_combined_builder;
    arrow::DoubleBuilder ms1_var_xcorr_shape_builder;
    arrow::DoubleBuilder ms1_var_xcorr_shape_contrast_builder;
    arrow::DoubleBuilder ms1_var_xcorr_shape_combined_builder;

    arrow::Int64Builder ms2_feature_id_builder;
    arrow::Int64Builder ms2_run_id_builder;
    arrow::DoubleBuilder ms2_area_builder;
    arrow::DoubleBuilder ms2_total_area_builder;
    arrow::DoubleBuilder ms2_apex_builder;
    arrow::DoubleBuilder ms2_exp_im_builder;
    arrow::DoubleBuilder ms2_exp_im_left_builder;
    arrow::DoubleBuilder ms2_exp_im_right_builder;
    arrow::DoubleBuilder ms2_delta_im_builder;
    arrow::DoubleBuilder ms2_total_mi_builder;
    arrow::DoubleBuilder ms2_var_bseries_builder;
    arrow::DoubleBuilder ms2_var_dotprod_builder;
    arrow::DoubleBuilder ms2_var_intensity_builder;
    arrow::DoubleBuilder ms2_var_iso_corr_builder;
    arrow::DoubleBuilder ms2_var_iso_overlap_builder;
    arrow::DoubleBuilder ms2_var_library_corr_builder;
    arrow::DoubleBuilder ms2_var_library_dotprod_builder;
    arrow::DoubleBuilder ms2_var_library_manhattan_builder;
    arrow::DoubleBuilder ms2_var_library_rmsd_builder;
    arrow::DoubleBuilder ms2_var_library_rootmeansquare_builder;
    arrow::DoubleBuilder ms2_var_library_sangle_builder;
    arrow::DoubleBuilder ms2_var_log_sn_builder;
    arrow::DoubleBuilder ms2_var_manhattan_builder;
    arrow::DoubleBuilder ms2_var_massdev_builder;
    arrow::DoubleBuilder ms2_var_massdev_weighted_builder;
    arrow::DoubleBuilder ms2_var_mi_builder;
    arrow::DoubleBuilder ms2_var_mi_weighted_builder;
    arrow::DoubleBuilder ms2_var_mi_ratio_builder;
    arrow::DoubleBuilder ms2_var_norm_rt_builder;
    arrow::DoubleBuilder ms2_var_xcorr_coelution_builder;
    arrow::DoubleBuilder ms2_var_xcorr_coelution_weighted_builder;
    arrow::DoubleBuilder ms2_var_xcorr_shape_builder;
    arrow::DoubleBuilder ms2_var_xcorr_shape_weighted_builder;
    arrow::DoubleBuilder ms2_var_yseries_builder;
    arrow::DoubleBuilder ms2_var_elution_model_fit_builder;
    arrow::DoubleBuilder ms2_var_im_xcorr_shape_builder;
    arrow::DoubleBuilder ms2_var_im_xcorr_coelution_builder;
    arrow::DoubleBuilder ms2_var_im_delta_builder;
    arrow::DoubleBuilder ms2_var_im_log_intensity_builder;

    arrow::Int64Builder ft_feature_id_builder;
    arrow::Int64Builder ft_run_id_builder;
    arrow::Int64Builder ft_transition_id_builder;
    arrow::DoubleBuilder ft_area_builder;
    arrow::DoubleBuilder ft_total_area_builder;
    arrow::DoubleBuilder ft_apex_int_builder;
    arrow::DoubleBuilder ft_apex_rt_builder;
    arrow::DoubleBuilder ft_rt_fwhm_builder;
    arrow::DoubleBuilder ft_masserror_builder;
    arrow::DoubleBuilder ft_total_mi_builder;
    arrow::DoubleBuilder ft_var_intensity_builder;
    arrow::DoubleBuilder ft_var_intensity_ratio_builder;
    arrow::DoubleBuilder ft_var_log_intensity_builder;
    arrow::DoubleBuilder ft_var_xcorr_coelution_builder;
    arrow::DoubleBuilder ft_var_xcorr_shape_builder;
    arrow::DoubleBuilder ft_var_log_sn_builder;
    arrow::DoubleBuilder ft_var_massdev_builder;
    arrow::DoubleBuilder ft_var_mi_builder;
    arrow::DoubleBuilder ft_var_mi_ratio_builder;
    arrow::DoubleBuilder ft_var_isotope_corr_builder;
    arrow::DoubleBuilder ft_var_isotope_overlap_builder;
    arrow::DoubleBuilder ft_exp_im_builder;
    arrow::DoubleBuilder ft_exp_im_left_builder;
    arrow::DoubleBuilder ft_exp_im_right_builder;
    arrow::DoubleBuilder ft_delta_im_builder;
    arrow::DoubleBuilder ft_var_im_delta_builder;
    arrow::DoubleBuilder ft_var_im_log_intensity_builder;
    arrow::DoubleBuilder ft_var_im_xcorr_coelution_contrast_builder;
    arrow::DoubleBuilder ft_var_im_xcorr_shape_contrast_builder;
    arrow::DoubleBuilder ft_var_im_xcorr_coelution_combined_builder;
    arrow::DoubleBuilder ft_var_im_xcorr_shape_combined_builder;
    arrow::DoubleBuilder ft_start_position_at_5_builder;
    arrow::DoubleBuilder ft_end_position_at_5_builder;
    arrow::DoubleBuilder ft_start_position_at_10_builder;
    arrow::DoubleBuilder ft_end_position_at_10_builder;
    arrow::DoubleBuilder ft_start_position_at_50_builder;
    arrow::DoubleBuilder ft_end_position_at_50_builder;
    arrow::DoubleBuilder ft_total_width_builder;
    arrow::DoubleBuilder ft_tailing_factor_builder;
    arrow::DoubleBuilder ft_asymmetry_factor_builder;
    arrow::DoubleBuilder ft_slope_of_baseline_builder;
    arrow::DoubleBuilder ft_baseline_delta_2_height_builder;
    arrow::DoubleBuilder ft_points_across_baseline_builder;
    arrow::DoubleBuilder ft_points_across_half_height_builder;

    arrow::Int64Builder fp_feature_id_builder;
    arrow::Int64Builder fp_run_id_builder;
    arrow::Int32Builder fp_isotope_builder;
    arrow::DoubleBuilder fp_area_builder;
    arrow::DoubleBuilder fp_apex_builder;

    for (Size idx = 0; idx < feature_map.size(); ++idx)
    {
      const Feature& feature = feature_map[idx];
      const int64_t feature_id = Internal::SqliteHelper::clearSignBit(feature.getUniqueId());
      if (!feature.metaValueExists("PeptideRef"))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Feature missing PeptideRef meta value");
      }
      const String peptide_ref = feature.getMetaValue("PeptideRef");
      auto precursor_it = compound_to_precursor.find(peptide_ref);
      if (precursor_it == compound_to_precursor.end())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Feature references unknown peptide", peptide_ref);
      }

      appendOrThrow_(feature_id_builder.Append(feature_id), "feature_id");
      appendOrThrow_(feature_run_id_builder.Append(run_id_clean), "run_id");
      appendOrThrow_(precursor_id_builder.Append(precursor_it->second), "precursor_id");
      appendOptionalFloat_(exp_rt_builder, true, feature.getRT(), "exp_rt");

      double value = 0.0;
      appendOptionalFloat_(exp_im_builder, extractMetaDouble_(feature, "im_drift", value), value, "exp_im");
      appendOptionalFloat_(norm_rt_builder, extractMetaDouble_(feature, "norm_RT", value), value, "norm_rt");
      appendOptionalFloat_(delta_rt_builder, extractMetaDouble_(feature, "delta_rt", value), value, "delta_rt");
      appendOptionalFloat_(left_width_builder, extractMetaDouble_(feature, "leftWidth", value), value, "left_width");
      appendOptionalFloat_(right_width_builder, extractMetaDouble_(feature, "rightWidth", value), value, "right_width");
      appendOptionalFloat_(exp_im_left_builder, extractMetaDouble_(feature, "im_drift_left", value), value, "exp_im_leftwidth");
      appendOptionalFloat_(exp_im_right_builder, extractMetaDouble_(feature, "im_drift_right", value), value, "exp_im_rightwidth");

      if (feature.metaValueExists("var_ms1_ppm_diff"))
      {
        appendOrThrow_(ms1_feature_id_builder.Append(feature_id), "feature_id");
        appendOrThrow_(ms1_run_id_builder.Append(run_id_clean), "run_id");
        appendOptionalFloat_(ms1_area_builder, extractMetaDouble_(feature, "ms1_area_intensity", value), value, "area_intensity");
        appendOptionalFloat_(ms1_apex_builder, extractMetaDouble_(feature, "ms1_apex_intensity", value), value, "apex_intensity");
        appendOptionalFloat_(ms1_exp_im_builder, extractMetaDouble_(feature, "im_ms1_drift", value), value, "exp_im");
        appendOptionalFloat_(ms1_delta_im_builder, extractMetaDouble_(feature, "im_ms1_delta", value), value, "delta_im");
        appendOptionalFloat_(ms1_var_massdev_builder, extractMetaDouble_(feature, "var_ms1_ppm_diff", value), value, "var_massdev_score");
        appendOptionalFloat_(ms1_var_im_ms1_delta_builder, extractMetaDouble_(feature, "var_im_ms1_delta_score", value), value, "var_im_ms1_delta_score");
        appendOptionalFloat_(ms1_var_mi_builder, extractMetaDouble_(feature, "var_ms1_mi_score", value), value, "var_mi_score");
        appendOptionalFloat_(ms1_var_mi_contrast_builder, extractMetaDouble_(feature, "var_ms1_mi_contrast_score", value), value, "var_mi_contrast_score");
        appendOptionalFloat_(ms1_var_mi_combined_builder, extractMetaDouble_(feature, "var_ms1_mi_combined_score", value), value, "var_mi_combined_score");
        appendOptionalFloat_(ms1_var_iso_corr_builder, extractMetaDouble_(feature, "var_ms1_isotope_correlation", value), value, "var_isotope_correlation_score");
        appendOptionalFloat_(ms1_var_iso_overlap_builder, extractMetaDouble_(feature, "var_ms1_isotope_overlap", value), value, "var_isotope_overlap_score");
        appendOptionalFloat_(ms1_var_xcorr_coelution_builder, extractMetaDouble_(feature, "var_ms1_xcorr_coelution", value), value, "var_xcorr_coelution");
        appendOptionalFloat_(ms1_var_xcorr_coelution_contrast_builder, extractMetaDouble_(feature, "var_ms1_xcorr_coelution_contrast", value), value, "var_xcorr_coelution_contrast");
        appendOptionalFloat_(ms1_var_xcorr_coelution_combined_builder, extractMetaDouble_(feature, "var_ms1_xcorr_coelution_combined", value), value, "var_xcorr_coelution_combined");
        appendOptionalFloat_(ms1_var_xcorr_shape_builder, extractMetaDouble_(feature, "var_ms1_xcorr_shape", value), value, "var_xcorr_shape");
        appendOptionalFloat_(ms1_var_xcorr_shape_contrast_builder, extractMetaDouble_(feature, "var_ms1_xcorr_shape_contrast", value), value, "var_xcorr_shape_contrast");
        appendOptionalFloat_(ms1_var_xcorr_shape_combined_builder, extractMetaDouble_(feature, "var_ms1_xcorr_shape_combined", value), value, "var_xcorr_shape_combined");
      }

      appendOrThrow_(ms2_feature_id_builder.Append(feature_id), "feature_id");
      appendOrThrow_(ms2_run_id_builder.Append(run_id_clean), "run_id");
      appendOptionalFloat_(ms2_area_builder, true, feature.getIntensity(), "area_intensity");
      appendOptionalFloat_(ms2_total_area_builder, extractMetaDouble_(feature, "total_xic", value), value, "total_area_intensity");
      appendOptionalFloat_(ms2_apex_builder, extractMetaDouble_(feature, "peak_apices_sum", value), value, "apex_intensity");
      appendOptionalFloat_(ms2_exp_im_builder, extractMetaDouble_(feature, "im_drift", value), value, "exp_im");
      appendOptionalFloat_(ms2_exp_im_left_builder, extractMetaDouble_(feature, "im_drift_left", value), value, "exp_im_leftwidth");
      appendOptionalFloat_(ms2_exp_im_right_builder, extractMetaDouble_(feature, "im_drift_right", value), value, "exp_im_rightwidth");
      appendOptionalFloat_(ms2_delta_im_builder, extractMetaDouble_(feature, "im_delta", value), value, "delta_im");
      appendOptionalFloat_(ms2_total_mi_builder, extractMetaDouble_(feature, "total_mi", value), value, "total_mi");
      appendOptionalFloat_(ms2_var_bseries_builder, extractMetaDouble_(feature, "var_bseries_score", value), value, "var_bseries_score");
      appendOptionalFloat_(ms2_var_dotprod_builder, extractMetaDouble_(feature, "var_dotprod_score", value), value, "var_dotprod_score");
      appendOptionalFloat_(ms2_var_intensity_builder, extractMetaDouble_(feature, "var_intensity_score", value), value, "var_intensity_score");
      appendOptionalFloat_(ms2_var_iso_corr_builder, extractMetaDouble_(feature, "var_isotope_correlation_score", value), value, "var_isotope_correlation_score");
      appendOptionalFloat_(ms2_var_iso_overlap_builder, extractMetaDouble_(feature, "var_isotope_overlap_score", value), value, "var_isotope_overlap_score");
      appendOptionalFloat_(ms2_var_library_corr_builder, extractMetaDouble_(feature, "var_library_corr", value), value, "var_library_corr");
      appendOptionalFloat_(ms2_var_library_dotprod_builder, extractMetaDouble_(feature, "var_library_dotprod", value), value, "var_library_dotprod");
      appendOptionalFloat_(ms2_var_library_manhattan_builder, extractMetaDouble_(feature, "var_library_manhattan", value), value, "var_library_manhattan");
      appendOptionalFloat_(ms2_var_library_rmsd_builder, extractMetaDouble_(feature, "var_library_rmsd", value), value, "var_library_rmsd");
      appendOptionalFloat_(ms2_var_library_rootmeansquare_builder, extractMetaDouble_(feature, "var_library_rootmeansquare", value), value, "var_library_rootmeansquare");
      appendOptionalFloat_(ms2_var_library_sangle_builder, extractMetaDouble_(feature, "var_library_sangle", value), value, "var_library_sangle");
      appendOptionalFloat_(ms2_var_log_sn_builder, extractMetaDouble_(feature, "var_log_sn_score", value), value, "var_log_sn_score");
      appendOptionalFloat_(ms2_var_manhattan_builder, extractMetaDouble_(feature, "var_manhatt_score", value), value, "var_manhattan_score");
      appendOptionalFloat_(ms2_var_massdev_builder, extractMetaDouble_(feature, "var_massdev_score", value), value, "var_massdev_score");
      appendOptionalFloat_(ms2_var_massdev_weighted_builder, extractMetaDouble_(feature, "var_massdev_score_weighted", value), value, "var_massdev_score_weighted");
      appendOptionalFloat_(ms2_var_mi_builder, extractMetaDouble_(feature, "var_mi_score", value), value, "var_mi_score");
      appendOptionalFloat_(ms2_var_mi_weighted_builder, extractMetaDouble_(feature, "var_mi_weighted_score", value), value, "var_mi_weighted_score");
      appendOptionalFloat_(ms2_var_mi_ratio_builder, extractMetaDouble_(feature, "var_mi_ratio_score", value), value, "var_mi_ratio_score");
      appendOptionalFloat_(ms2_var_norm_rt_builder, extractMetaDouble_(feature, "var_norm_rt_score", value), value, "var_norm_rt_score");
      appendOptionalFloat_(ms2_var_xcorr_coelution_builder, extractMetaDouble_(feature, "var_xcorr_coelution", value), value, "var_xcorr_coelution");
      appendOptionalFloat_(ms2_var_xcorr_coelution_weighted_builder, extractMetaDouble_(feature, "var_xcorr_coelution_weighted", value), value, "var_xcorr_coelution_weighted");
      appendOptionalFloat_(ms2_var_xcorr_shape_builder, extractMetaDouble_(feature, "var_xcorr_shape", value), value, "var_xcorr_shape");
      appendOptionalFloat_(ms2_var_xcorr_shape_weighted_builder, extractMetaDouble_(feature, "var_xcorr_shape_weighted", value), value, "var_xcorr_shape_weighted");
      appendOptionalFloat_(ms2_var_yseries_builder, extractMetaDouble_(feature, "var_yseries_score", value), value, "var_yseries_score");
      appendOptionalFloat_(ms2_var_elution_model_fit_builder, extractMetaDouble_(feature, "var_elution_model_fit_score", value), value, "var_elution_model_fit_score");
      appendOptionalFloat_(ms2_var_im_xcorr_shape_builder, extractMetaDouble_(feature, "var_im_xcorr_shape", value), value, "var_im_xcorr_shape");
      appendOptionalFloat_(ms2_var_im_xcorr_coelution_builder, extractMetaDouble_(feature, "var_im_xcorr_coelution", value), value, "var_im_xcorr_coelution");
      appendOptionalFloat_(ms2_var_im_delta_builder, extractMetaDouble_(feature, "var_im_delta_score", value), value, "var_im_delta_score");
      appendOptionalFloat_(ms2_var_im_log_intensity_builder, extractMetaDouble_(feature, "im_log_intensity", value), value, "var_im_log_intensity");

      auto masserror_ppm = getSeparateScore_(feature, "masserror_ppm");
      const auto& subordinates = feature.getSubordinates();
      for (Size i = 0; i < subordinates.size(); ++i)
      {
        const auto& sub_it = subordinates[i];
        if (!(sub_it.metaValueExists("FeatureLevel") && sub_it.getMetaValue("FeatureLevel") == "MS2"))
        {
          continue;
        }

        if (!sub_it.metaValueExists("native_id"))
        {
          continue;
        }
        const String native_id = sub_it.getMetaValue("native_id");
        int64_t transition_id_value = 0;
        if (!parseOptionalInt64_(native_id, transition_id_value))
        {
          auto it = transition_to_id.find(native_id);
          if (it == transition_to_id.end())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Transition references unknown id", native_id);
          }
          transition_id_value = it->second;
        }

        appendOrThrow_(ft_feature_id_builder.Append(feature_id), "feature_id");
        appendOrThrow_(ft_run_id_builder.Append(run_id_clean), "run_id");
        appendOrThrow_(ft_transition_id_builder.Append(transition_id_value), "transition_id");
        appendOptionalFloat_(ft_area_builder, true, sub_it.getIntensity(), "area_intensity");
        appendOptionalFloat_(ft_total_area_builder, extractMetaDouble_(sub_it, "total_xic", value), value, "total_area_intensity");
        appendOptionalFloat_(ft_apex_int_builder, extractMetaDouble_(sub_it, "peak_apex_int", value), value, "apex_intensity");
        appendOptionalFloat_(ft_apex_rt_builder, extractMetaDouble_(sub_it, "peak_apex_position", value), value, "apex_rt");
        appendOptionalFloat_(ft_rt_fwhm_builder, extractMetaDouble_(sub_it, "width_at_50", value), value, "rt_fwhm");
        appendOptionalFloatFromList_(ft_masserror_builder, masserror_ppm, i, "masserror_ppm");
        appendOptionalFloat_(ft_total_mi_builder, extractMetaDouble_(sub_it, "total_mi", value), value, "total_mi");
        appendOrThrow_(ft_var_intensity_builder.AppendNull(), "var_intensity_score");
        appendOrThrow_(ft_var_intensity_ratio_builder.AppendNull(), "var_intensity_ratio_score");
        appendOrThrow_(ft_var_log_intensity_builder.AppendNull(), "var_log_intensity");
        appendOrThrow_(ft_var_xcorr_coelution_builder.AppendNull(), "var_xcorr_coelution");
        appendOrThrow_(ft_var_xcorr_shape_builder.AppendNull(), "var_xcorr_shape");
        appendOrThrow_(ft_var_log_sn_builder.AppendNull(), "var_log_sn_score");
        appendOrThrow_(ft_var_massdev_builder.AppendNull(), "var_massdev_score");
        appendOrThrow_(ft_var_mi_builder.AppendNull(), "var_mi_score");
        appendOrThrow_(ft_var_mi_ratio_builder.AppendNull(), "var_mi_ratio_score");
        appendOrThrow_(ft_var_isotope_corr_builder.AppendNull(), "var_isotope_correlation_score");
        appendOrThrow_(ft_var_isotope_overlap_builder.AppendNull(), "var_isotope_overlap_score");
        appendOrThrow_(ft_exp_im_builder.AppendNull(), "exp_im");
        appendOrThrow_(ft_exp_im_left_builder.AppendNull(), "exp_im_leftwidth");
        appendOrThrow_(ft_exp_im_right_builder.AppendNull(), "exp_im_rightwidth");
        appendOrThrow_(ft_delta_im_builder.AppendNull(), "delta_im");
        appendOrThrow_(ft_var_im_delta_builder.AppendNull(), "var_im_delta_score");
        appendOrThrow_(ft_var_im_log_intensity_builder.AppendNull(), "var_im_log_intensity");
        appendOrThrow_(ft_var_im_xcorr_coelution_contrast_builder.AppendNull(), "var_im_xcorr_coelution_contrast");
        appendOrThrow_(ft_var_im_xcorr_shape_contrast_builder.AppendNull(), "var_im_xcorr_shape_contrast");
        appendOrThrow_(ft_var_im_xcorr_coelution_combined_builder.AppendNull(), "var_im_xcorr_coelution_combined");
        appendOrThrow_(ft_var_im_xcorr_shape_combined_builder.AppendNull(), "var_im_xcorr_shape_combined");

        const bool has_peak_shape = sub_it.metaValueExists("start_position_at_5");
        appendOptionalFloat_(ft_start_position_at_5_builder, has_peak_shape && extractMetaDouble_(sub_it, "start_position_at_5", value), value, "start_position_at_5");
        appendOptionalFloat_(ft_end_position_at_5_builder, has_peak_shape && extractMetaDouble_(sub_it, "end_position_at_5", value), value, "end_position_at_5");
        appendOptionalFloat_(ft_start_position_at_10_builder, has_peak_shape && extractMetaDouble_(sub_it, "start_position_at_10", value), value, "start_position_at_10");
        appendOptionalFloat_(ft_end_position_at_10_builder, has_peak_shape && extractMetaDouble_(sub_it, "end_position_at_10", value), value, "end_position_at_10");
        appendOptionalFloat_(ft_start_position_at_50_builder, has_peak_shape && extractMetaDouble_(sub_it, "start_position_at_50", value), value, "start_position_at_50");
        appendOptionalFloat_(ft_end_position_at_50_builder, has_peak_shape && extractMetaDouble_(sub_it, "end_position_at_50", value), value, "end_position_at_50");
        appendOptionalFloat_(ft_total_width_builder, has_peak_shape && extractMetaDouble_(sub_it, "total_width", value), value, "total_width");
        appendOptionalFloat_(ft_tailing_factor_builder, has_peak_shape && extractMetaDouble_(sub_it, "tailing_factor", value), value, "tailing_factor");
        appendOptionalFloat_(ft_asymmetry_factor_builder, has_peak_shape && extractMetaDouble_(sub_it, "asymmetry_factor", value), value, "asymmetry_factor");
        appendOptionalFloat_(ft_slope_of_baseline_builder, has_peak_shape && extractMetaDouble_(sub_it, "slope_of_baseline", value), value, "slope_of_baseline");
        appendOptionalFloat_(ft_baseline_delta_2_height_builder, has_peak_shape && extractMetaDouble_(sub_it, "baseline_delta_2_height", value), value, "baseline_delta_2_height");
        appendOptionalFloat_(ft_points_across_baseline_builder, has_peak_shape && extractMetaDouble_(sub_it, "points_across_baseline", value), value, "points_across_baseline");
        appendOptionalFloat_(ft_points_across_half_height_builder, has_peak_shape && extractMetaDouble_(sub_it, "points_across_half_height", value), value, "points_across_half_height");
      }
      else if (sub_it.metaValueExists("FeatureLevel") && sub_it.getMetaValue("FeatureLevel") == "MS1" && sub_it.getIntensity() > 0.0)
      {
        std::vector<String> precursor_id;
        String(sub_it.getMetaValue("native_id")).split(String("Precursor_i"), precursor_id);
        if (precursor_id.size() < 2)
        {
          continue;
        }
        int64_t isotope_value = 0;
        if (!parseOptionalInt64_(precursor_id[1], isotope_value))
        {
          continue;
        }

        appendOrThrow_(fp_feature_id_builder.Append(feature_id), "feature_id");
        appendOrThrow_(fp_run_id_builder.Append(run_id_clean), "run_id");
        appendOrThrow_(fp_isotope_builder.Append(static_cast<int32_t>(isotope_value)), "isotope");
        appendOptionalFloat_(fp_area_builder, true, sub_it.getIntensity(), "area_intensity");
        appendOptionalFloat_(fp_apex_builder, extractMetaDouble_(sub_it, "peak_apex_int", value), value, "apex_intensity");
      }

      if (enable_uis_scoring)
      {
        auto id_target_transition_names = getSeparateScore_(feature, "id_target_transition_names");
        auto id_target_area_intensity = getSeparateScore_(feature, "id_target_area_intensity");
        auto id_target_total_area_intensity = getSeparateScore_(feature, "id_target_total_area_intensity");
        auto id_target_apex_intensity = getSeparateScore_(feature, "id_target_apex_intensity");
        auto id_target_peak_apex_position = getSeparateScore_(feature, "id_target_peak_apex_position");
        auto id_target_peak_fwhm = getSeparateScore_(feature, "id_target_width_at_50");
        auto id_target_total_mi = getSeparateScore_(feature, "id_target_total_mi");
        auto id_target_intensity_score = getSeparateScore_(feature, "id_target_intensity_score");
        auto id_target_intensity_ratio_score = getSeparateScore_(feature, "id_target_intensity_ratio_score");
        auto id_target_log_intensity = getSeparateScore_(feature, "id_target_ind_log_intensity");
        auto id_target_xcorr_coelution = getSeparateScore_(feature, "id_target_ind_xcorr_coelution");
        auto id_target_xcorr_shape = getSeparateScore_(feature, "id_target_ind_xcorr_shape");
        auto id_target_log_sn_score = getSeparateScore_(feature, "id_target_ind_log_sn_score");
        auto id_target_massdev_score = getSeparateScore_(feature, "id_target_ind_massdev_score");
        auto id_target_mi_score = getSeparateScore_(feature, "id_target_ind_mi_score");
        auto id_target_mi_ratio_score = getSeparateScore_(feature, "id_target_ind_mi_ratio_score");
        auto id_target_isotope_correlation = getSeparateScore_(feature, "id_target_ind_isotope_correlation");
        auto id_target_isotope_overlap = getSeparateScore_(feature, "id_target_ind_isotope_overlap");
        auto id_target_im_drift = getSeparateScore_(feature, "id_target_ind_im_drift");
        auto id_target_im_drift_left = getSeparateScore_(feature, "id_target_ind_im_drift_left");
        auto id_target_im_drift_right = getSeparateScore_(feature, "id_target_ind_im_drift_right");
        auto id_target_im_delta = getSeparateScore_(feature, "id_target_ind_im_delta");
        auto id_target_im_delta_score = getSeparateScore_(feature, "id_target_ind_im_delta_score");
        auto id_target_im_log_intensity = getSeparateScore_(feature, "id_target_ind_im_log_intensity");
        auto id_target_im_contrast_coelution = getSeparateScore_(feature, "id_target_ind_im_contrast_coelution");
        auto id_target_im_contrast_shape = getSeparateScore_(feature, "id_target_ind_im_contrast_shape");
        auto id_target_im_sum_contrast_coelution = getSeparateScore_(feature, "id_target_ind_im_sum_contrast_coelution");
        auto id_target_im_sum_contrast_shape = getSeparateScore_(feature, "id_target_ind_im_sum_contrast_shape");

        auto id_target_ind_start_position_at_5 = getSeparateScore_(feature, "id_target_ind_start_position_at_5");
        const bool enable_target_peak_shape = !id_target_ind_start_position_at_5.empty() && id_target_ind_start_position_at_5[0] != "0";
        auto id_target_ind_end_position_at_5 = getSeparateScore_(feature, "id_target_ind_end_position_at_5");
        auto id_target_ind_start_position_at_10 = getSeparateScore_(feature, "id_target_ind_start_position_at_10");
        auto id_target_ind_end_position_at_10 = getSeparateScore_(feature, "id_target_ind_end_position_at_10");
        auto id_target_ind_start_position_at_50 = getSeparateScore_(feature, "id_target_ind_start_position_at_50");
        auto id_target_ind_end_position_at_50 = getSeparateScore_(feature, "id_target_ind_end_position_at_50");
        auto id_target_ind_total_width = getSeparateScore_(feature, "id_target_ind_total_width");
        auto id_target_ind_tailing_factor = getSeparateScore_(feature, "id_target_ind_tailing_factor");
        auto id_target_ind_asymmetry_factor = getSeparateScore_(feature, "id_target_ind_asymmetry_factor");
        auto id_target_ind_slope_of_baseline = getSeparateScore_(feature, "id_target_ind_slope_of_baseline");
        auto id_target_ind_baseline_delta_2_height = getSeparateScore_(feature, "id_target_ind_baseline_delta_2_height");
        auto id_target_ind_points_across_baseline = getSeparateScore_(feature, "id_target_ind_points_across_baseline");
        auto id_target_ind_points_across_half_height = getSeparateScore_(feature, "id_target_ind_points_across_half_height");

        for (Size i = 0; i < id_target_transition_names.size(); ++i)
        {
          const String transition_name = id_target_transition_names[i];
          auto it = transition_to_id.find(transition_name);
          if (it == transition_to_id.end()) continue;

          appendOrThrow_(ft_feature_id_builder.Append(feature_id), "feature_id");
          appendOrThrow_(ft_run_id_builder.Append(run_id_clean), "run_id");
          appendOrThrow_(ft_transition_id_builder.Append(it->second), "transition_id");
          appendOptionalFloatFromList_(ft_area_builder, id_target_area_intensity, i, "area_intensity");
          appendOptionalFloatFromList_(ft_total_area_builder, id_target_total_area_intensity, i, "total_area_intensity");
          appendOptionalFloatFromList_(ft_apex_int_builder, id_target_apex_intensity, i, "apex_intensity");
          appendOptionalFloatFromList_(ft_apex_rt_builder, id_target_peak_apex_position, i, "apex_rt");
          appendOptionalFloatFromList_(ft_rt_fwhm_builder, id_target_peak_fwhm, i, "rt_fwhm");
          appendOptionalFloatFromList_(ft_masserror_builder, id_target_massdev_score, i, "masserror_ppm");
          appendOptionalFloatFromList_(ft_total_mi_builder, id_target_total_mi, i, "total_mi");
          appendOptionalFloatFromList_(ft_var_intensity_builder, id_target_intensity_score, i, "var_intensity_score");
          appendOptionalFloatFromList_(ft_var_intensity_ratio_builder, id_target_intensity_ratio_score, i, "var_intensity_ratio_score");
          appendOptionalFloatFromList_(ft_var_log_intensity_builder, id_target_log_intensity, i, "var_log_intensity");
          appendOptionalFloatFromList_(ft_var_xcorr_coelution_builder, id_target_xcorr_coelution, i, "var_xcorr_coelution");
          appendOptionalFloatFromList_(ft_var_xcorr_shape_builder, id_target_xcorr_shape, i, "var_xcorr_shape");
          appendOptionalFloatFromList_(ft_var_log_sn_builder, id_target_log_sn_score, i, "var_log_sn_score");
          appendOptionalFloatFromList_(ft_var_massdev_builder, id_target_massdev_score, i, "var_massdev_score");
          appendOptionalFloatFromList_(ft_var_mi_builder, id_target_mi_score, i, "var_mi_score");
          appendOptionalFloatFromList_(ft_var_mi_ratio_builder, id_target_mi_ratio_score, i, "var_mi_ratio_score");
          appendOptionalFloatFromList_(ft_var_isotope_corr_builder, id_target_isotope_correlation, i, "var_isotope_correlation_score");
          appendOptionalFloatFromList_(ft_var_isotope_overlap_builder, id_target_isotope_overlap, i, "var_isotope_overlap_score");
          appendOptionalFloatFromList_(ft_exp_im_builder, id_target_im_drift, i, "exp_im");
          appendOptionalFloatFromList_(ft_exp_im_left_builder, id_target_im_drift_left, i, "exp_im_leftwidth");
          appendOptionalFloatFromList_(ft_exp_im_right_builder, id_target_im_drift_right, i, "exp_im_rightwidth");
          appendOptionalFloatFromList_(ft_delta_im_builder, id_target_im_delta, i, "delta_im");
          appendOptionalFloatFromList_(ft_var_im_delta_builder, id_target_im_delta_score, i, "var_im_delta_score");
          appendOptionalFloatFromList_(ft_var_im_log_intensity_builder, id_target_im_log_intensity, i, "var_im_log_intensity");
          appendOptionalFloatFromList_(ft_var_im_xcorr_coelution_contrast_builder, id_target_im_contrast_coelution, i, "var_im_xcorr_coelution_contrast");
          appendOptionalFloatFromList_(ft_var_im_xcorr_shape_contrast_builder, id_target_im_contrast_shape, i, "var_im_xcorr_shape_contrast");
          appendOptionalFloatFromList_(ft_var_im_xcorr_coelution_combined_builder, id_target_im_sum_contrast_coelution, i, "var_im_xcorr_coelution_combined");
          appendOptionalFloatFromList_(ft_var_im_xcorr_shape_combined_builder, id_target_im_sum_contrast_shape, i, "var_im_xcorr_shape_combined");
          if (enable_target_peak_shape)
          {
            appendOptionalFloatFromList_(ft_start_position_at_5_builder, id_target_ind_start_position_at_5, i, "start_position_at_5");
            appendOptionalFloatFromList_(ft_end_position_at_5_builder, id_target_ind_end_position_at_5, i, "end_position_at_5");
            appendOptionalFloatFromList_(ft_start_position_at_10_builder, id_target_ind_start_position_at_10, i, "start_position_at_10");
            appendOptionalFloatFromList_(ft_end_position_at_10_builder, id_target_ind_end_position_at_10, i, "end_position_at_10");
            appendOptionalFloatFromList_(ft_start_position_at_50_builder, id_target_ind_start_position_at_50, i, "start_position_at_50");
            appendOptionalFloatFromList_(ft_end_position_at_50_builder, id_target_ind_end_position_at_50, i, "end_position_at_50");
            appendOptionalFloatFromList_(ft_total_width_builder, id_target_ind_total_width, i, "total_width");
            appendOptionalFloatFromList_(ft_tailing_factor_builder, id_target_ind_tailing_factor, i, "tailing_factor");
            appendOptionalFloatFromList_(ft_asymmetry_factor_builder, id_target_ind_asymmetry_factor, i, "asymmetry_factor");
            appendOptionalFloatFromList_(ft_slope_of_baseline_builder, id_target_ind_slope_of_baseline, i, "slope_of_baseline");
            appendOptionalFloatFromList_(ft_baseline_delta_2_height_builder, id_target_ind_baseline_delta_2_height, i, "baseline_delta_2_height");
            appendOptionalFloatFromList_(ft_points_across_baseline_builder, id_target_ind_points_across_baseline, i, "points_across_baseline");
            appendOptionalFloatFromList_(ft_points_across_half_height_builder, id_target_ind_points_across_half_height, i, "points_across_half_height");
          }
          else
          {
            appendOrThrow_(ft_start_position_at_5_builder.AppendNull(), "start_position_at_5");
            appendOrThrow_(ft_end_position_at_5_builder.AppendNull(), "end_position_at_5");
            appendOrThrow_(ft_start_position_at_10_builder.AppendNull(), "start_position_at_10");
            appendOrThrow_(ft_end_position_at_10_builder.AppendNull(), "end_position_at_10");
            appendOrThrow_(ft_start_position_at_50_builder.AppendNull(), "start_position_at_50");
            appendOrThrow_(ft_end_position_at_50_builder.AppendNull(), "end_position_at_50");
            appendOrThrow_(ft_total_width_builder.AppendNull(), "total_width");
            appendOrThrow_(ft_tailing_factor_builder.AppendNull(), "tailing_factor");
            appendOrThrow_(ft_asymmetry_factor_builder.AppendNull(), "asymmetry_factor");
            appendOrThrow_(ft_slope_of_baseline_builder.AppendNull(), "slope_of_baseline");
            appendOrThrow_(ft_baseline_delta_2_height_builder.AppendNull(), "baseline_delta_2_height");
            appendOrThrow_(ft_points_across_baseline_builder.AppendNull(), "points_across_baseline");
            appendOrThrow_(ft_points_across_half_height_builder.AppendNull(), "points_across_half_height");
          }
        }

        auto id_decoy_transition_names = getSeparateScore_(feature, "id_decoy_transition_names");
        auto id_decoy_area_intensity = getSeparateScore_(feature, "id_decoy_area_intensity");
        auto id_decoy_total_area_intensity = getSeparateScore_(feature, "id_decoy_total_area_intensity");
        auto id_decoy_apex_intensity = getSeparateScore_(feature, "id_decoy_apex_intensity");
        auto id_decoy_peak_apex_position = getSeparateScore_(feature, "id_decoy_peak_apex_position");
        auto id_decoy_peak_fwhm = getSeparateScore_(feature, "id_decoy_width_at_50");
        auto id_decoy_total_mi = getSeparateScore_(feature, "id_decoy_total_mi");
        auto id_decoy_intensity_score = getSeparateScore_(feature, "id_decoy_intensity_score");
        auto id_decoy_intensity_ratio_score = getSeparateScore_(feature, "id_decoy_intensity_ratio_score");
        auto id_decoy_log_intensity = getSeparateScore_(feature, "id_decoy_ind_log_intensity");
        auto id_decoy_xcorr_coelution = getSeparateScore_(feature, "id_decoy_ind_xcorr_coelution");
        auto id_decoy_xcorr_shape = getSeparateScore_(feature, "id_decoy_ind_xcorr_shape");
        auto id_decoy_log_sn_score = getSeparateScore_(feature, "id_decoy_ind_log_sn_score");
        auto id_decoy_massdev_score = getSeparateScore_(feature, "id_decoy_ind_massdev_score");
        auto id_decoy_mi_score = getSeparateScore_(feature, "id_decoy_ind_mi_score");
        auto id_decoy_mi_ratio_score = getSeparateScore_(feature, "id_decoy_ind_mi_ratio_score");
        auto id_decoy_isotope_correlation = getSeparateScore_(feature, "id_decoy_ind_isotope_correlation");
        auto id_decoy_isotope_overlap = getSeparateScore_(feature, "id_decoy_ind_isotope_overlap");
        auto id_decoy_im_drift = getSeparateScore_(feature, "id_decoy_ind_im_drift");
        auto id_decoy_im_drift_left = getSeparateScore_(feature, "id_decoy_ind_im_drift_left");
        auto id_decoy_im_drift_right = getSeparateScore_(feature, "id_decoy_ind_im_drift_right");
        auto id_decoy_im_delta = getSeparateScore_(feature, "id_decoy_ind_im_delta");
        auto id_decoy_im_delta_score = getSeparateScore_(feature, "id_decoy_ind_im_delta_score");
        auto id_decoy_im_log_intensity = getSeparateScore_(feature, "id_decoy_ind_im_log_intensity");
        auto id_decoy_im_contrast_coelution = getSeparateScore_(feature, "id_decoy_ind_im_contrast_coelution");
        auto id_decoy_im_contrast_shape = getSeparateScore_(feature, "id_decoy_ind_im_contrast_shape");
        auto id_decoy_im_sum_contrast_coelution = getSeparateScore_(feature, "id_decoy_ind_im_sum_contrast_coelution");
        auto id_decoy_im_sum_contrast_shape = getSeparateScore_(feature, "id_decoy_ind_im_sum_contrast_shape");

        auto id_decoy_ind_start_position_at_5 = getSeparateScore_(feature, "id_decoy_ind_start_position_at_5");
        const bool enable_decoy_peak_shape = !id_decoy_ind_start_position_at_5.empty() && id_decoy_ind_start_position_at_5[0] != "0";
        auto id_decoy_ind_end_position_at_5 = getSeparateScore_(feature, "id_decoy_ind_end_position_at_5");
        auto id_decoy_ind_start_position_at_10 = getSeparateScore_(feature, "id_decoy_ind_start_position_at_10");
        auto id_decoy_ind_end_position_at_10 = getSeparateScore_(feature, "id_decoy_ind_end_position_at_10");
        auto id_decoy_ind_start_position_at_50 = getSeparateScore_(feature, "id_decoy_ind_start_position_at_50");
        auto id_decoy_ind_end_position_at_50 = getSeparateScore_(feature, "id_decoy_ind_end_position_at_50");
        auto id_decoy_ind_total_width = getSeparateScore_(feature, "id_decoy_ind_total_width");
        auto id_decoy_ind_tailing_factor = getSeparateScore_(feature, "id_decoy_ind_tailing_factor");
        auto id_decoy_ind_asymmetry_factor = getSeparateScore_(feature, "id_decoy_ind_asymmetry_factor");
        auto id_decoy_ind_slope_of_baseline = getSeparateScore_(feature, "id_decoy_ind_slope_of_baseline");
        auto id_decoy_ind_baseline_delta_2_height = getSeparateScore_(feature, "id_decoy_ind_baseline_delta_2_height");
        auto id_decoy_ind_points_across_baseline = getSeparateScore_(feature, "id_decoy_ind_points_across_baseline");
        auto id_decoy_ind_points_across_half_height = getSeparateScore_(feature, "id_decoy_ind_points_across_half_height");

        for (Size i = 0; i < id_decoy_transition_names.size(); ++i)
        {
          const String transition_name = id_decoy_transition_names[i];
          auto it = transition_to_id.find(transition_name);
          if (it == transition_to_id.end()) continue;

          appendOrThrow_(ft_feature_id_builder.Append(feature_id), "feature_id");
          appendOrThrow_(ft_run_id_builder.Append(run_id_clean), "run_id");
          appendOrThrow_(ft_transition_id_builder.Append(it->second), "transition_id");
          appendOptionalFloatFromList_(ft_area_builder, id_decoy_area_intensity, i, "area_intensity");
          appendOptionalFloatFromList_(ft_total_area_builder, id_decoy_total_area_intensity, i, "total_area_intensity");
          appendOptionalFloatFromList_(ft_apex_int_builder, id_decoy_apex_intensity, i, "apex_intensity");
          appendOptionalFloatFromList_(ft_apex_rt_builder, id_decoy_peak_apex_position, i, "apex_rt");
          appendOptionalFloatFromList_(ft_rt_fwhm_builder, id_decoy_peak_fwhm, i, "rt_fwhm");
          appendOptionalFloatFromList_(ft_masserror_builder, id_decoy_massdev_score, i, "masserror_ppm");
          appendOptionalFloatFromList_(ft_total_mi_builder, id_decoy_total_mi, i, "total_mi");
          appendOptionalFloatFromList_(ft_var_intensity_builder, id_decoy_intensity_score, i, "var_intensity_score");
          appendOptionalFloatFromList_(ft_var_intensity_ratio_builder, id_decoy_intensity_ratio_score, i, "var_intensity_ratio_score");
          appendOptionalFloatFromList_(ft_var_log_intensity_builder, id_decoy_log_intensity, i, "var_log_intensity");
          appendOptionalFloatFromList_(ft_var_xcorr_coelution_builder, id_decoy_xcorr_coelution, i, "var_xcorr_coelution");
          appendOptionalFloatFromList_(ft_var_xcorr_shape_builder, id_decoy_xcorr_shape, i, "var_xcorr_shape");
          appendOptionalFloatFromList_(ft_var_log_sn_builder, id_decoy_log_sn_score, i, "var_log_sn_score");
          appendOptionalFloatFromList_(ft_var_massdev_builder, id_decoy_massdev_score, i, "var_massdev_score");
          appendOptionalFloatFromList_(ft_var_mi_builder, id_decoy_mi_score, i, "var_mi_score");
          appendOptionalFloatFromList_(ft_var_mi_ratio_builder, id_decoy_mi_ratio_score, i, "var_mi_ratio_score");
          appendOptionalFloatFromList_(ft_var_isotope_corr_builder, id_decoy_isotope_correlation, i, "var_isotope_correlation_score");
          appendOptionalFloatFromList_(ft_var_isotope_overlap_builder, id_decoy_isotope_overlap, i, "var_isotope_overlap_score");
          appendOptionalFloatFromList_(ft_exp_im_builder, id_decoy_im_drift, i, "exp_im");
          appendOptionalFloatFromList_(ft_exp_im_left_builder, id_decoy_im_drift_left, i, "exp_im_leftwidth");
          appendOptionalFloatFromList_(ft_exp_im_right_builder, id_decoy_im_drift_right, i, "exp_im_rightwidth");
          appendOptionalFloatFromList_(ft_delta_im_builder, id_decoy_im_delta, i, "delta_im");
          appendOptionalFloatFromList_(ft_var_im_delta_builder, id_decoy_im_delta_score, i, "var_im_delta_score");
          appendOptionalFloatFromList_(ft_var_im_log_intensity_builder, id_decoy_im_log_intensity, i, "var_im_log_intensity");
          appendOptionalFloatFromList_(ft_var_im_xcorr_coelution_contrast_builder, id_decoy_im_contrast_coelution, i, "var_im_xcorr_coelution_contrast");
          appendOptionalFloatFromList_(ft_var_im_xcorr_shape_contrast_builder, id_decoy_im_contrast_shape, i, "var_im_xcorr_shape_contrast");
          appendOptionalFloatFromList_(ft_var_im_xcorr_coelution_combined_builder, id_decoy_im_sum_contrast_coelution, i, "var_im_xcorr_coelution_combined");
          appendOptionalFloatFromList_(ft_var_im_xcorr_shape_combined_builder, id_decoy_im_sum_contrast_shape, i, "var_im_xcorr_shape_combined");
          if (enable_decoy_peak_shape)
          {
            appendOptionalFloatFromList_(ft_start_position_at_5_builder, id_decoy_ind_start_position_at_5, i, "start_position_at_5");
            appendOptionalFloatFromList_(ft_end_position_at_5_builder, id_decoy_ind_end_position_at_5, i, "end_position_at_5");
            appendOptionalFloatFromList_(ft_start_position_at_10_builder, id_decoy_ind_start_position_at_10, i, "start_position_at_10");
            appendOptionalFloatFromList_(ft_end_position_at_10_builder, id_decoy_ind_end_position_at_10, i, "end_position_at_10");
            appendOptionalFloatFromList_(ft_start_position_at_50_builder, id_decoy_ind_start_position_at_50, i, "start_position_at_50");
            appendOptionalFloatFromList_(ft_end_position_at_50_builder, id_decoy_ind_end_position_at_50, i, "end_position_at_50");
            appendOptionalFloatFromList_(ft_total_width_builder, id_decoy_ind_total_width, i, "total_width");
            appendOptionalFloatFromList_(ft_tailing_factor_builder, id_decoy_ind_tailing_factor, i, "tailing_factor");
            appendOptionalFloatFromList_(ft_asymmetry_factor_builder, id_decoy_ind_asymmetry_factor, i, "asymmetry_factor");
            appendOptionalFloatFromList_(ft_slope_of_baseline_builder, id_decoy_ind_slope_of_baseline, i, "slope_of_baseline");
            appendOptionalFloatFromList_(ft_baseline_delta_2_height_builder, id_decoy_ind_baseline_delta_2_height, i, "baseline_delta_2_height");
            appendOptionalFloatFromList_(ft_points_across_baseline_builder, id_decoy_ind_points_across_baseline, i, "points_across_baseline");
            appendOptionalFloatFromList_(ft_points_across_half_height_builder, id_decoy_ind_points_across_half_height, i, "points_across_half_height");
          }
          else
          {
            appendOrThrow_(ft_start_position_at_5_builder.AppendNull(), "start_position_at_5");
            appendOrThrow_(ft_end_position_at_5_builder.AppendNull(), "end_position_at_5");
            appendOrThrow_(ft_start_position_at_10_builder.AppendNull(), "start_position_at_10");
            appendOrThrow_(ft_end_position_at_10_builder.AppendNull(), "end_position_at_10");
            appendOrThrow_(ft_start_position_at_50_builder.AppendNull(), "start_position_at_50");
            appendOrThrow_(ft_end_position_at_50_builder.AppendNull(), "end_position_at_50");
            appendOrThrow_(ft_total_width_builder.AppendNull(), "total_width");
            appendOrThrow_(ft_tailing_factor_builder.AppendNull(), "tailing_factor");
            appendOrThrow_(ft_asymmetry_factor_builder.AppendNull(), "asymmetry_factor");
            appendOrThrow_(ft_slope_of_baseline_builder.AppendNull(), "slope_of_baseline");
            appendOrThrow_(ft_baseline_delta_2_height_builder.AppendNull(), "baseline_delta_2_height");
            appendOrThrow_(ft_points_across_baseline_builder.AppendNull(), "points_across_baseline");
            appendOrThrow_(ft_points_across_half_height_builder.AppendNull(), "points_across_half_height");
          }
        }
      }
    }

    {
      auto features_schema = arrow::schema({
        arrow::field("feature_id", arrow::int64()),
        arrow::field("run_id", arrow::int64()),
        arrow::field("precursor_id", arrow::int64()),
        arrow::field("exp_rt", arrow::float64()),
        arrow::field("exp_im", arrow::float64()),
        arrow::field("norm_rt", arrow::float64()),
        arrow::field("delta_rt", arrow::float64()),
        arrow::field("left_width", arrow::float64()),
        arrow::field("right_width", arrow::float64()),
        arrow::field("exp_im_leftwidth", arrow::float64()),
        arrow::field("exp_im_rightwidth", arrow::float64())
      });
      auto features_table = arrow::Table::Make(features_schema, {
        finishArray_(feature_id_builder, "feature_id"),
        finishArray_(feature_run_id_builder, "run_id"),
        finishArray_(precursor_id_builder, "precursor_id"),
        finishArray_(exp_rt_builder, "exp_rt"),
        finishArray_(exp_im_builder, "exp_im"),
        finishArray_(norm_rt_builder, "norm_rt"),
        finishArray_(delta_rt_builder, "delta_rt"),
        finishArray_(left_width_builder, "left_width"),
        finishArray_(right_width_builder, "right_width"),
        finishArray_(exp_im_left_builder, "exp_im_leftwidth"),
        finishArray_(exp_im_right_builder, "exp_im_rightwidth")
      });
      writeParquetTable_(features_table, run_path + "/features.parquet");
    }

    if (ms1_feature_id_builder.length() > 0)
    {
      auto ms1_schema = arrow::schema({
        arrow::field("feature_id", arrow::int64()),
        arrow::field("run_id", arrow::int64()),
        arrow::field("area_intensity", arrow::float64()),
        arrow::field("apex_intensity", arrow::float64()),
        arrow::field("exp_im", arrow::float64()),
        arrow::field("delta_im", arrow::float64()),
        arrow::field("var_massdev_score", arrow::float64()),
        arrow::field("var_im_ms1_delta_score", arrow::float64()),
        arrow::field("var_mi_score", arrow::float64()),
        arrow::field("var_mi_contrast_score", arrow::float64()),
        arrow::field("var_mi_combined_score", arrow::float64()),
        arrow::field("var_isotope_correlation_score", arrow::float64()),
        arrow::field("var_isotope_overlap_score", arrow::float64()),
        arrow::field("var_xcorr_coelution", arrow::float64()),
        arrow::field("var_xcorr_coelution_contrast", arrow::float64()),
        arrow::field("var_xcorr_coelution_combined", arrow::float64()),
        arrow::field("var_xcorr_shape", arrow::float64()),
        arrow::field("var_xcorr_shape_contrast", arrow::float64()),
        arrow::field("var_xcorr_shape_combined", arrow::float64())
      });
      auto ms1_table = arrow::Table::Make(ms1_schema, {
        finishArray_(ms1_feature_id_builder, "feature_id"),
        finishArray_(ms1_run_id_builder, "run_id"),
        finishArray_(ms1_area_builder, "area_intensity"),
        finishArray_(ms1_apex_builder, "apex_intensity"),
        finishArray_(ms1_exp_im_builder, "exp_im"),
        finishArray_(ms1_delta_im_builder, "delta_im"),
        finishArray_(ms1_var_massdev_builder, "var_massdev_score"),
        finishArray_(ms1_var_im_ms1_delta_builder, "var_im_ms1_delta_score"),
        finishArray_(ms1_var_mi_builder, "var_mi_score"),
        finishArray_(ms1_var_mi_contrast_builder, "var_mi_contrast_score"),
        finishArray_(ms1_var_mi_combined_builder, "var_mi_combined_score"),
        finishArray_(ms1_var_iso_corr_builder, "var_isotope_correlation_score"),
        finishArray_(ms1_var_iso_overlap_builder, "var_isotope_overlap_score"),
        finishArray_(ms1_var_xcorr_coelution_builder, "var_xcorr_coelution"),
        finishArray_(ms1_var_xcorr_coelution_contrast_builder, "var_xcorr_coelution_contrast"),
        finishArray_(ms1_var_xcorr_coelution_combined_builder, "var_xcorr_coelution_combined"),
        finishArray_(ms1_var_xcorr_shape_builder, "var_xcorr_shape"),
        finishArray_(ms1_var_xcorr_shape_contrast_builder, "var_xcorr_shape_contrast"),
        finishArray_(ms1_var_xcorr_shape_combined_builder, "var_xcorr_shape_combined")
      });
      writeParquetTable_(ms1_table, run_path + "/feature_ms1.parquet");
    }

    if (fp_feature_id_builder.length() > 0)
    {
      auto fp_schema = arrow::schema({
        arrow::field("feature_id", arrow::int64()),
        arrow::field("run_id", arrow::int64()),
        arrow::field("isotope", arrow::int32()),
        arrow::field("area_intensity", arrow::float64()),
        arrow::field("apex_intensity", arrow::float64())
      });
      auto fp_table = arrow::Table::Make(fp_schema, {
        finishArray_(fp_feature_id_builder, "feature_id"),
        finishArray_(fp_run_id_builder, "run_id"),
        finishArray_(fp_isotope_builder, "isotope"),
        finishArray_(fp_area_builder, "area_intensity"),
        finishArray_(fp_apex_builder, "apex_intensity")
      });
      writeParquetTable_(fp_table, run_path + "/feature_precursor.parquet");
    }

    {
      auto ms2_schema = arrow::schema({
        arrow::field("feature_id", arrow::int64()),
        arrow::field("run_id", arrow::int64()),
        arrow::field("area_intensity", arrow::float64()),
        arrow::field("total_area_intensity", arrow::float64()),
        arrow::field("apex_intensity", arrow::float64()),
        arrow::field("exp_im", arrow::float64()),
        arrow::field("exp_im_leftwidth", arrow::float64()),
        arrow::field("exp_im_rightwidth", arrow::float64()),
        arrow::field("delta_im", arrow::float64()),
        arrow::field("total_mi", arrow::float64()),
        arrow::field("var_bseries_score", arrow::float64()),
        arrow::field("var_dotprod_score", arrow::float64()),
        arrow::field("var_intensity_score", arrow::float64()),
        arrow::field("var_isotope_correlation_score", arrow::float64()),
        arrow::field("var_isotope_overlap_score", arrow::float64()),
        arrow::field("var_library_corr", arrow::float64()),
        arrow::field("var_library_dotprod", arrow::float64()),
        arrow::field("var_library_manhattan", arrow::float64()),
        arrow::field("var_library_rmsd", arrow::float64()),
        arrow::field("var_library_rootmeansquare", arrow::float64()),
        arrow::field("var_library_sangle", arrow::float64()),
        arrow::field("var_log_sn_score", arrow::float64()),
        arrow::field("var_manhattan_score", arrow::float64()),
        arrow::field("var_massdev_score", arrow::float64()),
        arrow::field("var_massdev_score_weighted", arrow::float64()),
        arrow::field("var_mi_score", arrow::float64()),
        arrow::field("var_mi_weighted_score", arrow::float64()),
        arrow::field("var_mi_ratio_score", arrow::float64()),
        arrow::field("var_norm_rt_score", arrow::float64()),
        arrow::field("var_xcorr_coelution", arrow::float64()),
        arrow::field("var_xcorr_coelution_weighted", arrow::float64()),
        arrow::field("var_xcorr_shape", arrow::float64()),
        arrow::field("var_xcorr_shape_weighted", arrow::float64()),
        arrow::field("var_yseries_score", arrow::float64()),
        arrow::field("var_elution_model_fit_score", arrow::float64()),
        arrow::field("var_im_xcorr_shape", arrow::float64()),
        arrow::field("var_im_xcorr_coelution", arrow::float64()),
        arrow::field("var_im_delta_score", arrow::float64()),
        arrow::field("var_im_log_intensity", arrow::float64())
      });
      auto ms2_table = arrow::Table::Make(ms2_schema, {
        finishArray_(ms2_feature_id_builder, "feature_id"),
        finishArray_(ms2_run_id_builder, "run_id"),
        finishArray_(ms2_area_builder, "area_intensity"),
        finishArray_(ms2_total_area_builder, "total_area_intensity"),
        finishArray_(ms2_apex_builder, "apex_intensity"),
        finishArray_(ms2_exp_im_builder, "exp_im"),
        finishArray_(ms2_exp_im_left_builder, "exp_im_leftwidth"),
        finishArray_(ms2_exp_im_right_builder, "exp_im_rightwidth"),
        finishArray_(ms2_delta_im_builder, "delta_im"),
        finishArray_(ms2_total_mi_builder, "total_mi"),
        finishArray_(ms2_var_bseries_builder, "var_bseries_score"),
        finishArray_(ms2_var_dotprod_builder, "var_dotprod_score"),
        finishArray_(ms2_var_intensity_builder, "var_intensity_score"),
        finishArray_(ms2_var_iso_corr_builder, "var_isotope_correlation_score"),
        finishArray_(ms2_var_iso_overlap_builder, "var_isotope_overlap_score"),
        finishArray_(ms2_var_library_corr_builder, "var_library_corr"),
        finishArray_(ms2_var_library_dotprod_builder, "var_library_dotprod"),
        finishArray_(ms2_var_library_manhattan_builder, "var_library_manhattan"),
        finishArray_(ms2_var_library_rmsd_builder, "var_library_rmsd"),
        finishArray_(ms2_var_library_rootmeansquare_builder, "var_library_rootmeansquare"),
        finishArray_(ms2_var_library_sangle_builder, "var_library_sangle"),
        finishArray_(ms2_var_log_sn_builder, "var_log_sn_score"),
        finishArray_(ms2_var_manhattan_builder, "var_manhattan_score"),
        finishArray_(ms2_var_massdev_builder, "var_massdev_score"),
        finishArray_(ms2_var_massdev_weighted_builder, "var_massdev_score_weighted"),
        finishArray_(ms2_var_mi_builder, "var_mi_score"),
        finishArray_(ms2_var_mi_weighted_builder, "var_mi_weighted_score"),
        finishArray_(ms2_var_mi_ratio_builder, "var_mi_ratio_score"),
        finishArray_(ms2_var_norm_rt_builder, "var_norm_rt_score"),
        finishArray_(ms2_var_xcorr_coelution_builder, "var_xcorr_coelution"),
        finishArray_(ms2_var_xcorr_coelution_weighted_builder, "var_xcorr_coelution_weighted"),
        finishArray_(ms2_var_xcorr_shape_builder, "var_xcorr_shape"),
        finishArray_(ms2_var_xcorr_shape_weighted_builder, "var_xcorr_shape_weighted"),
        finishArray_(ms2_var_yseries_builder, "var_yseries_score"),
        finishArray_(ms2_var_elution_model_fit_builder, "var_elution_model_fit_score"),
        finishArray_(ms2_var_im_xcorr_shape_builder, "var_im_xcorr_shape"),
        finishArray_(ms2_var_im_xcorr_coelution_builder, "var_im_xcorr_coelution"),
        finishArray_(ms2_var_im_delta_builder, "var_im_delta_score"),
        finishArray_(ms2_var_im_log_intensity_builder, "var_im_log_intensity")
      });
      writeParquetTable_(ms2_table, run_path + "/feature_ms2.parquet");
    }

    {
      auto ft_schema = arrow::schema({
        arrow::field("feature_id", arrow::int64()),
        arrow::field("run_id", arrow::int64()),
        arrow::field("transition_id", arrow::int64()),
        arrow::field("area_intensity", arrow::float64()),
        arrow::field("total_area_intensity", arrow::float64()),
        arrow::field("apex_intensity", arrow::float64()),
        arrow::field("apex_rt", arrow::float64()),
        arrow::field("rt_fwhm", arrow::float64()),
        arrow::field("masserror_ppm", arrow::float64()),
        arrow::field("total_mi", arrow::float64()),
        arrow::field("var_intensity_score", arrow::float64()),
        arrow::field("var_intensity_ratio_score", arrow::float64()),
        arrow::field("var_log_intensity", arrow::float64()),
        arrow::field("var_xcorr_coelution", arrow::float64()),
        arrow::field("var_xcorr_shape", arrow::float64()),
        arrow::field("var_log_sn_score", arrow::float64()),
        arrow::field("var_massdev_score", arrow::float64()),
        arrow::field("var_mi_score", arrow::float64()),
        arrow::field("var_mi_ratio_score", arrow::float64()),
        arrow::field("var_isotope_correlation_score", arrow::float64()),
        arrow::field("var_isotope_overlap_score", arrow::float64()),
        arrow::field("exp_im", arrow::float64()),
        arrow::field("exp_im_leftwidth", arrow::float64()),
        arrow::field("exp_im_rightwidth", arrow::float64()),
        arrow::field("delta_im", arrow::float64()),
        arrow::field("var_im_delta_score", arrow::float64()),
        arrow::field("var_im_log_intensity", arrow::float64()),
        arrow::field("var_im_xcorr_coelution_contrast", arrow::float64()),
        arrow::field("var_im_xcorr_shape_contrast", arrow::float64()),
        arrow::field("var_im_xcorr_coelution_combined", arrow::float64()),
        arrow::field("var_im_xcorr_shape_combined", arrow::float64()),
        arrow::field("start_position_at_5", arrow::float64()),
        arrow::field("end_position_at_5", arrow::float64()),
        arrow::field("start_position_at_10", arrow::float64()),
        arrow::field("end_position_at_10", arrow::float64()),
        arrow::field("start_position_at_50", arrow::float64()),
        arrow::field("end_position_at_50", arrow::float64()),
        arrow::field("total_width", arrow::float64()),
        arrow::field("tailing_factor", arrow::float64()),
        arrow::field("asymmetry_factor", arrow::float64()),
        arrow::field("slope_of_baseline", arrow::float64()),
        arrow::field("baseline_delta_2_height", arrow::float64()),
        arrow::field("points_across_baseline", arrow::float64()),
        arrow::field("points_across_half_height", arrow::float64())
      });
      auto ft_table = arrow::Table::Make(ft_schema, {
        finishArray_(ft_feature_id_builder, "feature_id"),
        finishArray_(ft_run_id_builder, "run_id"),
        finishArray_(ft_transition_id_builder, "transition_id"),
        finishArray_(ft_area_builder, "area_intensity"),
        finishArray_(ft_total_area_builder, "total_area_intensity"),
        finishArray_(ft_apex_int_builder, "apex_intensity"),
        finishArray_(ft_apex_rt_builder, "apex_rt"),
        finishArray_(ft_rt_fwhm_builder, "rt_fwhm"),
        finishArray_(ft_masserror_builder, "masserror_ppm"),
        finishArray_(ft_total_mi_builder, "total_mi"),
        finishArray_(ft_var_intensity_builder, "var_intensity_score"),
        finishArray_(ft_var_intensity_ratio_builder, "var_intensity_ratio_score"),
        finishArray_(ft_var_log_intensity_builder, "var_log_intensity"),
        finishArray_(ft_var_xcorr_coelution_builder, "var_xcorr_coelution"),
        finishArray_(ft_var_xcorr_shape_builder, "var_xcorr_shape"),
        finishArray_(ft_var_log_sn_builder, "var_log_sn_score"),
        finishArray_(ft_var_massdev_builder, "var_massdev_score"),
        finishArray_(ft_var_mi_builder, "var_mi_score"),
        finishArray_(ft_var_mi_ratio_builder, "var_mi_ratio_score"),
        finishArray_(ft_var_isotope_corr_builder, "var_isotope_correlation_score"),
        finishArray_(ft_var_isotope_overlap_builder, "var_isotope_overlap_score"),
        finishArray_(ft_exp_im_builder, "exp_im"),
        finishArray_(ft_exp_im_left_builder, "exp_im_leftwidth"),
        finishArray_(ft_exp_im_right_builder, "exp_im_rightwidth"),
        finishArray_(ft_delta_im_builder, "delta_im"),
        finishArray_(ft_var_im_delta_builder, "var_im_delta_score"),
        finishArray_(ft_var_im_log_intensity_builder, "var_im_log_intensity"),
        finishArray_(ft_var_im_xcorr_coelution_contrast_builder, "var_im_xcorr_coelution_contrast"),
        finishArray_(ft_var_im_xcorr_shape_contrast_builder, "var_im_xcorr_shape_contrast"),
        finishArray_(ft_var_im_xcorr_coelution_combined_builder, "var_im_xcorr_coelution_combined"),
        finishArray_(ft_var_im_xcorr_shape_combined_builder, "var_im_xcorr_shape_combined"),
        finishArray_(ft_start_position_at_5_builder, "start_position_at_5"),
        finishArray_(ft_end_position_at_5_builder, "end_position_at_5"),
        finishArray_(ft_start_position_at_10_builder, "start_position_at_10"),
        finishArray_(ft_end_position_at_10_builder, "end_position_at_10"),
        finishArray_(ft_start_position_at_50_builder, "start_position_at_50"),
        finishArray_(ft_end_position_at_50_builder, "end_position_at_50"),
        finishArray_(ft_total_width_builder, "total_width"),
        finishArray_(ft_tailing_factor_builder, "tailing_factor"),
        finishArray_(ft_asymmetry_factor_builder, "asymmetry_factor"),
        finishArray_(ft_slope_of_baseline_builder, "slope_of_baseline"),
        finishArray_(ft_baseline_delta_2_height_builder, "baseline_delta_2_height"),
        finishArray_(ft_points_across_baseline_builder, "points_across_baseline"),
        finishArray_(ft_points_across_half_height_builder, "points_across_half_height")
      });
      writeParquetTable_(ft_table, run_path + "/feature_transition.parquet");
    }

    if (zip_output)
    {
      zipParquetDirectory_(base_dir, output_path);
    }
#endif
  }

} // namespace OpenMS
