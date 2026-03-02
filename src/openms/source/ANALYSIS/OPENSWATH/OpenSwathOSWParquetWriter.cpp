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
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <fstream>
#include <cmath>
#include <future>
#include <limits>
#include <unordered_map>
#include <vector>

#include <QString>

#ifdef WITH_PARQUET
#ifdef signals
#undef signals
#endif
#include <arrow/api.h>
#include <parquet/file_reader.h>
#endif

namespace OpenMS
{
  namespace
  {
    using OpenMS::Size;

#ifdef WITH_PARQUET
    void appendOptionalFloat_(arrow::DoubleBuilder& builder, bool has_value, double value, const char* column)
    {
      if (!has_value || !std::isfinite(value))
      {
        ParquetFile::appendOrThrow(builder.AppendNull(), column);
      }
      else
      {
        ParquetFile::appendOrThrow(builder.Append(value), column);
      }
    }

    void appendOptionalInt_(arrow::Int64Builder& builder, bool has_value, int64_t value, const char* column)
    {
      if (!has_value)
      {
        ParquetFile::appendOrThrow(builder.AppendNull(), column);
      }
      else
      {
        ParquetFile::appendOrThrow(builder.Append(value), column);
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

    void appendFeatureScore_(arrow::DoubleBuilder& builder,
                             const Feature& feature,
                             const std::string& key,
                             const char* column,
                             double fallback_value = std::numeric_limits<double>::quiet_NaN(),
                             bool use_fallback = false)
    {
      double value = 0.0;
      if (extractMetaDouble_(feature, key, value))
      {
        appendOptionalFloat_(builder, true, value, column);
      }
      else if (use_fallback && std::isfinite(fallback_value))
      {
        appendOptionalFloat_(builder, true, fallback_value, column);
      }
      else
      {
        ParquetFile::appendOrThrow(builder.AppendNull(), column);
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
        ParquetFile::appendOrThrow(builder.AppendNull(), column);
      }
    }

    struct RunEntry
    {
      int64_t run_id = 0;
      String filename;
    };

    struct RunCounts
    {
      int64_t features = 0;
      int64_t feature_precursor = 0;
      int64_t feature_transition = 0;
    };

    int64_t getParquetRowCount_(const String& filename)
    {
      return OpenMS::ParquetFile::rowCount(filename);
    }

    std::vector<RunEntry> readRuns_(const String& runs_parquet)
    {
      std::vector<RunEntry> runs;
      if (!File::exists(runs_parquet))
      {
        return runs;
      }

      auto runs_table = ParquetFile::readTable(runs_parquet);
      auto run_id_array = std::static_pointer_cast<arrow::Int64Array>(ParquetFile::getColumn(runs_table, "run_id"));
      auto filename_array = std::static_pointer_cast<arrow::StringArray>(ParquetFile::getColumn(runs_table, "filename"));

      const int64_t num_rows = runs_table->num_rows();
      runs.reserve(num_rows);
      for (int64_t row = 0; row < num_rows; ++row)
      {
        RunEntry entry;
        entry.run_id = run_id_array->Value(row);
        entry.filename = filename_array->GetString(row);
        runs.push_back(entry);
      }

      return runs;
    }

    std::string jsonEscape_(const String& input)
    {
      return OpenMS::ParquetFile::jsonEscape(input);
    }

    void writeMetadata_(const String& base_dir,
                        const std::vector<RunEntry>& runs,
                        const std::vector<RunCounts>& run_counts,
                        const RunCounts& total_counts)
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
          << "    \"openms_version\": \"" << jsonEscape_(VersionInfo::getVersion()) << "\",\n"
          << "    \"build_time\": \"" << jsonEscape_(VersionInfo::getTime()) << "\",\n"
          << "    \"tool\": {\"name\": \"OpenSwathWorkflow\", \"version\": \"" << jsonEscape_(VersionInfo::getVersion()) << "\"},\n"
          << "    \"counts\": {\n"
          << "      \"runs\": " << runs.size() << ",\n"
          << "      \"features\": " << total_counts.features << ",\n"
          << "      \"feature_precursor\": " << total_counts.feature_precursor << ",\n"
          << "      \"feature_transition\": " << total_counts.feature_transition << "\n"
          << "    }\n"
          << "  },\n"
          << "  \"runs\": [\n";

      if (runs.size() != run_counts.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Run count metadata mismatch", String(run_counts.size()));
      }

      for (Size i = 0; i < runs.size(); ++i)
      {
        const RunEntry& run = runs[i];
        const RunCounts& counts = run_counts[i];
        out << "    {\"id\": " << run.run_id
            << ", \"filename\": \"" << jsonEscape_(run.filename) << "\""
            << ", \"counts\": {\"features\": " << counts.features
            << ", \"feature_precursor\": " << counts.feature_precursor
            << ", \"feature_transition\": " << counts.feature_transition << "}}";
        if (i + 1 < runs.size())
        {
          out << ",";
        }
        out << "\n";
      }
      out << "  ]\n"
          << "}\n";
    }

    void waitForWriteTasks_(std::vector<std::future<void>>& tasks)
    {
      for (auto& task : tasks)
      {
        task.get();
      }
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
    const String base_dir = output_path;

    const String library_dir = base_dir + "/library";
    const String runs_dir = base_dir + "/runs";
    if (!File::exists(base_dir))
    {
      File::makeDir(base_dir);
    }
    if (!File::exists(runs_dir))
    {
      File::makeDir(runs_dir);
    }

    const String precursors_path = library_dir + "/precursors.parquet";
    const String transitions_path = library_dir + "/transitions.parquet";
    const bool library_ready = File::exists(precursors_path) && File::exists(transitions_path);
    // If library files exist, perform a basic compatibility check to avoid
    // silently reusing an incompatible library when appending runs. A mismatch
    // between the existing library tables and the provided `assay_library`
    // (counts differ) is a strong signal of incompatibility and can create
    // broken foreign-key relationships between run-level files and the
    // library tables. Fail fast with a clear error message.
    if (library_ready)
    {
      const int64_t existing_precursors = getParquetRowCount_(precursors_path);
      const int64_t existing_transitions = getParquetRowCount_(transitions_path);
      const int64_t expected_precursors = static_cast<int64_t>(assay_library.compounds.size());
      const int64_t expected_transitions = static_cast<int64_t>(assay_library.transitions.size());
      if (existing_precursors != expected_precursors || existing_transitions != expected_transitions)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Existing library at '" + library_dir + "' appears incompatible with provided assay_library. Please rebuild the library or use a different output path.",
                                      String("existing_precursors=") + String(existing_precursors) + ", expected_precursors=" + String(expected_precursors));
      }
    }
    if (!library_ready)
    {
      const String library_tmp_dir = base_dir + "/library_tmp";
      if (File::exists(library_tmp_dir))
      {
        File::removeDirRecursively(library_tmp_dir);
      }
      if (File::exists(library_dir))
      {
        File::removeDirRecursively(library_dir);
      }
      File::makeDir(library_tmp_dir);
      File::makeDir(library_dir);
      TransitionParquetFile().convertLightTargetedExperimentToParquet(library_tmp_dir, assay_library);
      File::copyDirRecursively(String(library_tmp_dir + "/library").toQString(), library_dir.toQString());
      File::removeDirRecursively(library_tmp_dir);
    }

    const String run_path = runs_dir + "/run_id=" + String(run_id_clean);
    if (File::exists(run_path))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Run output directory already exists", run_path);
    }
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

    arrow::DoubleBuilder ms1_area_builder;
    arrow::DoubleBuilder ms1_apex_builder;
    arrow::DoubleBuilder ms1_exp_im_builder;
    arrow::DoubleBuilder ms1_delta_im_builder;
    arrow::DoubleBuilder var_ms1_massdev_builder;
    arrow::DoubleBuilder var_ms1_im_ms1_delta_builder;
    arrow::DoubleBuilder var_ms1_mi_builder;
    arrow::DoubleBuilder var_ms1_mi_contrast_builder;
    arrow::DoubleBuilder var_ms1_mi_combined_builder;
    arrow::DoubleBuilder var_ms1_iso_corr_builder;
    arrow::DoubleBuilder var_ms1_iso_overlap_builder;
    arrow::DoubleBuilder var_ms1_xcorr_coelution_builder;
    arrow::DoubleBuilder var_ms1_xcorr_coelution_contrast_builder;
    arrow::DoubleBuilder var_ms1_xcorr_coelution_combined_builder;
    arrow::DoubleBuilder var_ms1_xcorr_shape_builder;
    arrow::DoubleBuilder var_ms1_xcorr_shape_contrast_builder;
    arrow::DoubleBuilder var_ms1_xcorr_shape_combined_builder;

    arrow::DoubleBuilder ms2_area_builder;
    arrow::DoubleBuilder ms2_total_area_builder;
    arrow::DoubleBuilder ms2_apex_builder;
    arrow::DoubleBuilder ms2_exp_im_builder;
    arrow::DoubleBuilder ms2_exp_im_left_builder;
    arrow::DoubleBuilder ms2_exp_im_right_builder;
    arrow::DoubleBuilder ms2_delta_im_builder;
    arrow::DoubleBuilder ms2_total_mi_builder;
    arrow::DoubleBuilder var_ms2_bseries_builder;
    arrow::DoubleBuilder var_ms2_dotprod_builder;
    arrow::DoubleBuilder var_ms2_intensity_builder;
    arrow::DoubleBuilder var_ms2_iso_corr_builder;
    arrow::DoubleBuilder var_ms2_iso_overlap_builder;
    arrow::DoubleBuilder var_ms2_library_corr_builder;
    arrow::DoubleBuilder var_ms2_library_dotprod_builder;
    arrow::DoubleBuilder var_ms2_library_manhattan_builder;
    arrow::DoubleBuilder var_ms2_library_rmsd_builder;
    arrow::DoubleBuilder var_ms2_library_rootmeansquare_builder;
    arrow::DoubleBuilder var_ms2_library_sangle_builder;
    arrow::DoubleBuilder var_ms2_log_sn_builder;
    arrow::DoubleBuilder var_ms2_manhattan_builder;
    arrow::DoubleBuilder var_ms2_massdev_builder;
    arrow::DoubleBuilder var_ms2_massdev_weighted_builder;
    arrow::DoubleBuilder var_ms2_mi_builder;
    arrow::DoubleBuilder var_ms2_mi_weighted_builder;
    arrow::DoubleBuilder var_ms2_mi_ratio_builder;
    arrow::DoubleBuilder var_ms2_norm_rt_builder;
    arrow::DoubleBuilder var_ms2_xcorr_coelution_builder;
    arrow::DoubleBuilder var_ms2_xcorr_coelution_weighted_builder;
    arrow::DoubleBuilder var_ms2_xcorr_shape_builder;
    arrow::DoubleBuilder var_ms2_xcorr_shape_weighted_builder;
    arrow::DoubleBuilder var_ms2_yseries_builder;
    arrow::DoubleBuilder var_ms2_elution_model_fit_builder;
    arrow::DoubleBuilder var_ms2_im_xcorr_shape_builder;
    arrow::DoubleBuilder var_ms2_im_xcorr_coelution_builder;
    arrow::DoubleBuilder var_ms2_im_delta_builder;
    arrow::DoubleBuilder var_ms2_im_log_intensity_builder;

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

      ParquetFile::appendOrThrow(feature_id_builder.Append(feature_id), "feature_id");
      ParquetFile::appendOrThrow(feature_run_id_builder.Append(run_id_clean), "run_id");
      ParquetFile::appendOrThrow(precursor_id_builder.Append(precursor_it->second), "precursor_id");
      appendOptionalFloat_(exp_rt_builder, true, feature.getRT(), "exp_rt");

      double value = 0.0;
      appendFeatureScore_(exp_im_builder, feature, "im_drift", "exp_im");
      appendFeatureScore_(norm_rt_builder, feature, "norm_RT", "norm_rt", -1.0, true);
      appendFeatureScore_(delta_rt_builder, feature, "delta_rt", "delta_rt", -1.0, true);
      appendFeatureScore_(left_width_builder, feature, "leftWidth", "left_width");
      appendFeatureScore_(right_width_builder, feature, "rightWidth", "right_width");
      appendFeatureScore_(exp_im_left_builder, feature, "im_drift_left", "exp_im_leftwidth");
      appendFeatureScore_(exp_im_right_builder, feature, "im_drift_right", "exp_im_rightwidth");

      const bool has_ms1 = feature.metaValueExists("var_ms1_ppm_diff");
      if (has_ms1)
      {
        appendOptionalFloat_(ms1_area_builder, extractMetaDouble_(feature, "ms1_area_intensity", value), value, "ms1_area_intensity");
        appendOptionalFloat_(ms1_apex_builder, extractMetaDouble_(feature, "ms1_apex_intensity", value), value, "ms1_apex_intensity");
        appendOptionalFloat_(ms1_exp_im_builder, extractMetaDouble_(feature, "im_ms1_drift", value), value, "ms1_exp_im");
        appendOptionalFloat_(ms1_delta_im_builder, extractMetaDouble_(feature, "im_ms1_delta", value), value, "ms1_delta_im");
        appendOptionalFloat_(var_ms1_massdev_builder, extractMetaDouble_(feature, "var_ms1_ppm_diff", value), value, "var_ms1_massdev_score");
        appendOptionalFloat_(var_ms1_im_ms1_delta_builder, extractMetaDouble_(feature, "var_im_ms1_delta_score", value), value, "var_ms1_im_ms1_delta_score");
        appendOptionalFloat_(var_ms1_mi_builder, extractMetaDouble_(feature, "var_ms1_mi_score", value), value, "var_ms1_mi_score");
        appendOptionalFloat_(var_ms1_mi_contrast_builder, extractMetaDouble_(feature, "var_ms1_mi_contrast_score", value), value, "var_ms1_mi_contrast_score");
        appendOptionalFloat_(var_ms1_mi_combined_builder, extractMetaDouble_(feature, "var_ms1_mi_combined_score", value), value, "var_ms1_mi_combined_score");
        appendOptionalFloat_(var_ms1_iso_corr_builder, extractMetaDouble_(feature, "var_ms1_isotope_correlation", value), value, "var_ms1_isotope_correlation_score");
        appendOptionalFloat_(var_ms1_iso_overlap_builder, extractMetaDouble_(feature, "var_ms1_isotope_overlap", value), value, "var_ms1_isotope_overlap_score");
        appendOptionalFloat_(var_ms1_xcorr_coelution_builder, extractMetaDouble_(feature, "var_ms1_xcorr_coelution", value), value, "var_ms1_xcorr_coelution");
        appendOptionalFloat_(var_ms1_xcorr_coelution_contrast_builder, extractMetaDouble_(feature, "var_ms1_xcorr_coelution_contrast", value), value, "var_ms1_xcorr_coelution_contrast");
        appendOptionalFloat_(var_ms1_xcorr_coelution_combined_builder, extractMetaDouble_(feature, "var_ms1_xcorr_coelution_combined", value), value, "var_ms1_xcorr_coelution_combined");
        appendOptionalFloat_(var_ms1_xcorr_shape_builder, extractMetaDouble_(feature, "var_ms1_xcorr_shape", value), value, "var_ms1_xcorr_shape");
        appendOptionalFloat_(var_ms1_xcorr_shape_contrast_builder, extractMetaDouble_(feature, "var_ms1_xcorr_shape_contrast", value), value, "var_ms1_xcorr_shape_contrast");
        appendOptionalFloat_(var_ms1_xcorr_shape_combined_builder, extractMetaDouble_(feature, "var_ms1_xcorr_shape_combined", value), value, "var_ms1_xcorr_shape_combined");
      }
      else
      {
        ParquetFile::appendOrThrow(ms1_area_builder.AppendNull(), "ms1_area_intensity");
        ParquetFile::appendOrThrow(ms1_apex_builder.AppendNull(), "ms1_apex_intensity");
        ParquetFile::appendOrThrow(ms1_exp_im_builder.AppendNull(), "ms1_exp_im");
        ParquetFile::appendOrThrow(ms1_delta_im_builder.AppendNull(), "ms1_delta_im");
        ParquetFile::appendOrThrow(var_ms1_massdev_builder.AppendNull(), "var_ms1_massdev_score");
        ParquetFile::appendOrThrow(var_ms1_im_ms1_delta_builder.AppendNull(), "var_ms1_im_ms1_delta_score");
        ParquetFile::appendOrThrow(var_ms1_mi_builder.AppendNull(), "var_ms1_mi_score");
        ParquetFile::appendOrThrow(var_ms1_mi_contrast_builder.AppendNull(), "var_ms1_mi_contrast_score");
        ParquetFile::appendOrThrow(var_ms1_mi_combined_builder.AppendNull(), "var_ms1_mi_combined_score");
        ParquetFile::appendOrThrow(var_ms1_iso_corr_builder.AppendNull(), "var_ms1_isotope_correlation_score");
        ParquetFile::appendOrThrow(var_ms1_iso_overlap_builder.AppendNull(), "var_ms1_isotope_overlap_score");
        ParquetFile::appendOrThrow(var_ms1_xcorr_coelution_builder.AppendNull(), "var_ms1_xcorr_coelution");
        ParquetFile::appendOrThrow(var_ms1_xcorr_coelution_contrast_builder.AppendNull(), "var_ms1_xcorr_coelution_contrast");
        ParquetFile::appendOrThrow(var_ms1_xcorr_coelution_combined_builder.AppendNull(), "var_ms1_xcorr_coelution_combined");
        ParquetFile::appendOrThrow(var_ms1_xcorr_shape_builder.AppendNull(), "var_ms1_xcorr_shape");
        ParquetFile::appendOrThrow(var_ms1_xcorr_shape_contrast_builder.AppendNull(), "var_ms1_xcorr_shape_contrast");
        ParquetFile::appendOrThrow(var_ms1_xcorr_shape_combined_builder.AppendNull(), "var_ms1_xcorr_shape_combined");
      }

      appendOptionalFloat_(ms2_area_builder, true, feature.getIntensity(), "ms2_area_intensity");
      appendOptionalFloat_(ms2_total_area_builder, extractMetaDouble_(feature, "total_xic", value), value, "ms2_total_area_intensity");
      appendOptionalFloat_(ms2_apex_builder, extractMetaDouble_(feature, "peak_apices_sum", value), value, "ms2_apex_intensity");
      appendOptionalFloat_(ms2_exp_im_builder, extractMetaDouble_(feature, "im_drift", value), value, "ms2_exp_im");
      appendOptionalFloat_(ms2_exp_im_left_builder, extractMetaDouble_(feature, "im_drift_left", value), value, "ms2_exp_im_leftwidth");
      appendOptionalFloat_(ms2_exp_im_right_builder, extractMetaDouble_(feature, "im_drift_right", value), value, "ms2_exp_im_rightwidth");
      appendOptionalFloat_(ms2_delta_im_builder, extractMetaDouble_(feature, "im_delta", value), value, "ms2_delta_im");
      appendOptionalFloat_(ms2_total_mi_builder, extractMetaDouble_(feature, "total_mi", value), value, "ms2_total_mi");
      appendOptionalFloat_(var_ms2_bseries_builder, extractMetaDouble_(feature, "var_bseries_score", value), value, "var_ms2_bseries_score");
      appendOptionalFloat_(var_ms2_dotprod_builder, extractMetaDouble_(feature, "var_dotprod_score", value), value, "var_ms2_dotprod_score");
      appendOptionalFloat_(var_ms2_intensity_builder, extractMetaDouble_(feature, "var_intensity_score", value), value, "var_ms2_intensity_score");
      appendOptionalFloat_(var_ms2_iso_corr_builder, extractMetaDouble_(feature, "var_isotope_correlation_score", value), value, "var_ms2_isotope_correlation_score");
      appendOptionalFloat_(var_ms2_iso_overlap_builder, extractMetaDouble_(feature, "var_isotope_overlap_score", value), value, "var_ms2_isotope_overlap_score");
      appendOptionalFloat_(var_ms2_library_corr_builder, extractMetaDouble_(feature, "var_library_corr", value), value, "var_ms2_library_corr");
      appendOptionalFloat_(var_ms2_library_dotprod_builder, extractMetaDouble_(feature, "var_library_dotprod", value), value, "var_ms2_library_dotprod");
      appendOptionalFloat_(var_ms2_library_manhattan_builder, extractMetaDouble_(feature, "var_library_manhattan", value), value, "var_ms2_library_manhattan");
      appendOptionalFloat_(var_ms2_library_rmsd_builder, extractMetaDouble_(feature, "var_library_rmsd", value), value, "var_ms2_library_rmsd");
      appendOptionalFloat_(var_ms2_library_rootmeansquare_builder, extractMetaDouble_(feature, "var_library_rootmeansquare", value), value, "var_ms2_library_rootmeansquare");
      appendOptionalFloat_(var_ms2_library_sangle_builder, extractMetaDouble_(feature, "var_library_sangle", value), value, "var_ms2_library_sangle");
      appendOptionalFloat_(var_ms2_log_sn_builder, extractMetaDouble_(feature, "var_log_sn_score", value), value, "var_ms2_log_sn_score");
      appendOptionalFloat_(var_ms2_manhattan_builder, extractMetaDouble_(feature, "var_manhatt_score", value), value, "var_ms2_manhattan_score");
      appendOptionalFloat_(var_ms2_massdev_builder, extractMetaDouble_(feature, "var_massdev_score", value), value, "var_ms2_massdev_score");
      appendOptionalFloat_(var_ms2_massdev_weighted_builder, extractMetaDouble_(feature, "var_massdev_score_weighted", value), value, "var_ms2_massdev_score_weighted");
      appendOptionalFloat_(var_ms2_mi_builder, extractMetaDouble_(feature, "var_mi_score", value), value, "var_ms2_mi_score");
      appendOptionalFloat_(var_ms2_mi_weighted_builder, extractMetaDouble_(feature, "var_mi_weighted_score", value), value, "var_ms2_mi_weighted_score");
      appendOptionalFloat_(var_ms2_mi_ratio_builder, extractMetaDouble_(feature, "var_mi_ratio_score", value), value, "var_ms2_mi_ratio_score");
      appendOptionalFloat_(var_ms2_norm_rt_builder, extractMetaDouble_(feature, "var_norm_rt_score", value), value, "var_ms2_norm_rt_score");
      appendOptionalFloat_(var_ms2_xcorr_coelution_builder, extractMetaDouble_(feature, "var_xcorr_coelution", value), value, "var_ms2_xcorr_coelution");
      appendOptionalFloat_(var_ms2_xcorr_coelution_weighted_builder, extractMetaDouble_(feature, "var_xcorr_coelution_weighted", value), value, "var_ms2_xcorr_coelution_weighted");
      appendOptionalFloat_(var_ms2_xcorr_shape_builder, extractMetaDouble_(feature, "var_xcorr_shape", value), value, "var_ms2_xcorr_shape");
      appendOptionalFloat_(var_ms2_xcorr_shape_weighted_builder, extractMetaDouble_(feature, "var_xcorr_shape_weighted", value), value, "var_ms2_xcorr_shape_weighted");
      appendOptionalFloat_(var_ms2_yseries_builder, extractMetaDouble_(feature, "var_yseries_score", value), value, "var_ms2_yseries_score");
      appendOptionalFloat_(var_ms2_elution_model_fit_builder, extractMetaDouble_(feature, "var_elution_model_fit_score", value), value, "var_ms2_elution_model_fit_score");
      appendOptionalFloat_(var_ms2_im_xcorr_shape_builder, extractMetaDouble_(feature, "var_im_xcorr_shape", value), value, "var_ms2_im_xcorr_shape");
      appendOptionalFloat_(var_ms2_im_xcorr_coelution_builder, extractMetaDouble_(feature, "var_im_xcorr_coelution", value), value, "var_ms2_im_xcorr_coelution");
      appendOptionalFloat_(var_ms2_im_delta_builder, extractMetaDouble_(feature, "var_im_delta_score", value), value, "var_ms2_im_delta_score");
      appendOptionalFloat_(var_ms2_im_log_intensity_builder, extractMetaDouble_(feature, "im_log_intensity", value), value, "var_ms2_im_log_intensity");

      auto masserror_ppm = getSeparateScore_(feature, "masserror_ppm");
      const auto& subordinates = feature.getSubordinates();
      for (Size i = 0; i < subordinates.size(); ++i)
      {
        const auto& sub_it = subordinates[i];
        if (sub_it.metaValueExists("FeatureLevel") && sub_it.getMetaValue("FeatureLevel") == "MS2")
        {
          if (!sub_it.metaValueExists("native_id"))
          {
            continue;
          }
          const String native_id = sub_it.getMetaValue("native_id");
          int64_t transition_id_value = 0;
          // Prefer explicit name lookup in the library mapping. Only use the
          // numeric fallback when the native_id is numeric and the id falls
          // inside the valid id domain (1..transition_id-1). This prevents
          // accidental acceptance of arbitrary numeric strings that are not
          // actually present in the library table.
          auto it = transition_to_id.find(native_id);
          if (it != transition_to_id.end())
          {
            transition_id_value = it->second;
          }
          else if (parseOptionalInt64_(native_id, transition_id_value))
          {
            // numeric fallback is only valid if present in library id domain
            if (transition_id_value <= 0 || transition_id_value >= transition_id)
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Transition references unknown numeric id", native_id);
            }
          }
          else
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Transition references unknown id", native_id);
          }

          ParquetFile::appendOrThrow(ft_feature_id_builder.Append(feature_id), "feature_id");
          ParquetFile::appendOrThrow(ft_run_id_builder.Append(run_id_clean), "run_id");
          ParquetFile::appendOrThrow(ft_transition_id_builder.Append(transition_id_value), "transition_id");
          appendOptionalFloat_(ft_area_builder, true, sub_it.getIntensity(), "area_intensity");
          appendOptionalFloat_(ft_total_area_builder, extractMetaDouble_(sub_it, "total_xic", value), value, "total_area_intensity");
          appendOptionalFloat_(ft_apex_int_builder, extractMetaDouble_(sub_it, "peak_apex_int", value), value, "apex_intensity");
          appendOptionalFloat_(ft_apex_rt_builder, extractMetaDouble_(sub_it, "peak_apex_position", value), value, "apex_rt");
          appendOptionalFloat_(ft_rt_fwhm_builder, extractMetaDouble_(sub_it, "width_at_50", value), value, "rt_fwhm");
          appendOptionalFloatFromList_(ft_masserror_builder, masserror_ppm, i, "masserror_ppm");
          appendOptionalFloat_(ft_total_mi_builder, extractMetaDouble_(sub_it, "total_mi", value), value, "total_mi");
          ParquetFile::appendOrThrow(ft_var_intensity_builder.AppendNull(), "var_intensity_score");
          ParquetFile::appendOrThrow(ft_var_intensity_ratio_builder.AppendNull(), "var_intensity_ratio_score");
          ParquetFile::appendOrThrow(ft_var_log_intensity_builder.AppendNull(), "var_log_intensity");
          ParquetFile::appendOrThrow(ft_var_xcorr_coelution_builder.AppendNull(), "var_xcorr_coelution");
          ParquetFile::appendOrThrow(ft_var_xcorr_shape_builder.AppendNull(), "var_xcorr_shape");
          ParquetFile::appendOrThrow(ft_var_log_sn_builder.AppendNull(), "var_log_sn_score");
          ParquetFile::appendOrThrow(ft_var_massdev_builder.AppendNull(), "var_massdev_score");
          ParquetFile::appendOrThrow(ft_var_mi_builder.AppendNull(), "var_mi_score");
          ParquetFile::appendOrThrow(ft_var_mi_ratio_builder.AppendNull(), "var_mi_ratio_score");
          ParquetFile::appendOrThrow(ft_var_isotope_corr_builder.AppendNull(), "var_isotope_correlation_score");
          ParquetFile::appendOrThrow(ft_var_isotope_overlap_builder.AppendNull(), "var_isotope_overlap_score");
          ParquetFile::appendOrThrow(ft_exp_im_builder.AppendNull(), "exp_im");
          ParquetFile::appendOrThrow(ft_exp_im_left_builder.AppendNull(), "exp_im_leftwidth");
          ParquetFile::appendOrThrow(ft_exp_im_right_builder.AppendNull(), "exp_im_rightwidth");
          ParquetFile::appendOrThrow(ft_delta_im_builder.AppendNull(), "delta_im");
          ParquetFile::appendOrThrow(ft_var_im_delta_builder.AppendNull(), "var_im_delta_score");
          ParquetFile::appendOrThrow(ft_var_im_log_intensity_builder.AppendNull(), "var_im_log_intensity");
          ParquetFile::appendOrThrow(ft_var_im_xcorr_coelution_contrast_builder.AppendNull(), "var_im_xcorr_coelution_contrast");
          ParquetFile::appendOrThrow(ft_var_im_xcorr_shape_contrast_builder.AppendNull(), "var_im_xcorr_shape_contrast");
          ParquetFile::appendOrThrow(ft_var_im_xcorr_coelution_combined_builder.AppendNull(), "var_im_xcorr_coelution_combined");
          ParquetFile::appendOrThrow(ft_var_im_xcorr_shape_combined_builder.AppendNull(), "var_im_xcorr_shape_combined");

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
          ParquetFile::appendOrThrow(fp_feature_id_builder.Append(feature_id), "feature_id");
          ParquetFile::appendOrThrow(fp_run_id_builder.Append(run_id_clean), "run_id");
          ParquetFile::appendOrThrow(fp_isotope_builder.Append(static_cast<int32_t>(isotope_value)), "precursor_isotope");
          appendOptionalFloat_(fp_area_builder, true, sub_it.getIntensity(), "precursor_area_intensity");
          appendOptionalFloat_(fp_apex_builder, extractMetaDouble_(sub_it, "peak_apex_int", value), value, "precursor_apex_intensity");
        }
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

          ParquetFile::appendOrThrow(ft_feature_id_builder.Append(feature_id), "feature_id");
          ParquetFile::appendOrThrow(ft_run_id_builder.Append(run_id_clean), "run_id");
          ParquetFile::appendOrThrow(ft_transition_id_builder.Append(it->second), "transition_id");
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
            ParquetFile::appendOrThrow(ft_start_position_at_5_builder.AppendNull(), "start_position_at_5");
            ParquetFile::appendOrThrow(ft_end_position_at_5_builder.AppendNull(), "end_position_at_5");
            ParquetFile::appendOrThrow(ft_start_position_at_10_builder.AppendNull(), "start_position_at_10");
            ParquetFile::appendOrThrow(ft_end_position_at_10_builder.AppendNull(), "end_position_at_10");
            ParquetFile::appendOrThrow(ft_start_position_at_50_builder.AppendNull(), "start_position_at_50");
            ParquetFile::appendOrThrow(ft_end_position_at_50_builder.AppendNull(), "end_position_at_50");
            ParquetFile::appendOrThrow(ft_total_width_builder.AppendNull(), "total_width");
            ParquetFile::appendOrThrow(ft_tailing_factor_builder.AppendNull(), "tailing_factor");
            ParquetFile::appendOrThrow(ft_asymmetry_factor_builder.AppendNull(), "asymmetry_factor");
            ParquetFile::appendOrThrow(ft_slope_of_baseline_builder.AppendNull(), "slope_of_baseline");
            ParquetFile::appendOrThrow(ft_baseline_delta_2_height_builder.AppendNull(), "baseline_delta_2_height");
            ParquetFile::appendOrThrow(ft_points_across_baseline_builder.AppendNull(), "points_across_baseline");
            ParquetFile::appendOrThrow(ft_points_across_half_height_builder.AppendNull(), "points_across_half_height");
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

          ParquetFile::appendOrThrow(ft_feature_id_builder.Append(feature_id), "feature_id");
          ParquetFile::appendOrThrow(ft_run_id_builder.Append(run_id_clean), "run_id");
          ParquetFile::appendOrThrow(ft_transition_id_builder.Append(it->second), "transition_id");
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
            ParquetFile::appendOrThrow(ft_start_position_at_5_builder.AppendNull(), "start_position_at_5");
            ParquetFile::appendOrThrow(ft_end_position_at_5_builder.AppendNull(), "end_position_at_5");
            ParquetFile::appendOrThrow(ft_start_position_at_10_builder.AppendNull(), "start_position_at_10");
            ParquetFile::appendOrThrow(ft_end_position_at_10_builder.AppendNull(), "end_position_at_10");
            ParquetFile::appendOrThrow(ft_start_position_at_50_builder.AppendNull(), "start_position_at_50");
            ParquetFile::appendOrThrow(ft_end_position_at_50_builder.AppendNull(), "end_position_at_50");
            ParquetFile::appendOrThrow(ft_total_width_builder.AppendNull(), "total_width");
            ParquetFile::appendOrThrow(ft_tailing_factor_builder.AppendNull(), "tailing_factor");
            ParquetFile::appendOrThrow(ft_asymmetry_factor_builder.AppendNull(), "asymmetry_factor");
            ParquetFile::appendOrThrow(ft_slope_of_baseline_builder.AppendNull(), "slope_of_baseline");
            ParquetFile::appendOrThrow(ft_baseline_delta_2_height_builder.AppendNull(), "baseline_delta_2_height");
            ParquetFile::appendOrThrow(ft_points_across_baseline_builder.AppendNull(), "points_across_baseline");
            ParquetFile::appendOrThrow(ft_points_across_half_height_builder.AppendNull(), "points_across_half_height");
          }
        }
      }
    }

    std::vector<std::future<void>> write_tasks;
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
        arrow::field("exp_im_rightwidth", arrow::float64()),
        arrow::field("ms1_area_intensity", arrow::float64()),
        arrow::field("ms1_apex_intensity", arrow::float64()),
        arrow::field("ms1_exp_im", arrow::float64()),
        arrow::field("ms1_delta_im", arrow::float64()),
        arrow::field("var_ms1_massdev_score", arrow::float64()),
        arrow::field("var_ms1_im_ms1_delta_score", arrow::float64()),
        arrow::field("var_ms1_mi_score", arrow::float64()),
        arrow::field("var_ms1_mi_contrast_score", arrow::float64()),
        arrow::field("var_ms1_mi_combined_score", arrow::float64()),
        arrow::field("var_ms1_isotope_correlation_score", arrow::float64()),
        arrow::field("var_ms1_isotope_overlap_score", arrow::float64()),
        arrow::field("var_ms1_xcorr_coelution", arrow::float64()),
        arrow::field("var_ms1_xcorr_coelution_contrast", arrow::float64()),
        arrow::field("var_ms1_xcorr_coelution_combined", arrow::float64()),
        arrow::field("var_ms1_xcorr_shape", arrow::float64()),
        arrow::field("var_ms1_xcorr_shape_contrast", arrow::float64()),
        arrow::field("var_ms1_xcorr_shape_combined", arrow::float64()),
        arrow::field("ms2_area_intensity", arrow::float64()),
        arrow::field("ms2_total_area_intensity", arrow::float64()),
        arrow::field("ms2_apex_intensity", arrow::float64()),
        arrow::field("ms2_exp_im", arrow::float64()),
        arrow::field("ms2_exp_im_leftwidth", arrow::float64()),
        arrow::field("ms2_exp_im_rightwidth", arrow::float64()),
        arrow::field("ms2_delta_im", arrow::float64()),
        arrow::field("ms2_total_mi", arrow::float64()),
        arrow::field("var_ms2_bseries_score", arrow::float64()),
        arrow::field("var_ms2_dotprod_score", arrow::float64()),
        arrow::field("var_ms2_intensity_score", arrow::float64()),
        arrow::field("var_ms2_isotope_correlation_score", arrow::float64()),
        arrow::field("var_ms2_isotope_overlap_score", arrow::float64()),
        arrow::field("var_ms2_library_corr", arrow::float64()),
        arrow::field("var_ms2_library_dotprod", arrow::float64()),
        arrow::field("var_ms2_library_manhattan", arrow::float64()),
        arrow::field("var_ms2_library_rmsd", arrow::float64()),
        arrow::field("var_ms2_library_rootmeansquare", arrow::float64()),
        arrow::field("var_ms2_library_sangle", arrow::float64()),
        arrow::field("var_ms2_log_sn_score", arrow::float64()),
        arrow::field("var_ms2_manhattan_score", arrow::float64()),
        arrow::field("var_ms2_massdev_score", arrow::float64()),
        arrow::field("var_ms2_massdev_score_weighted", arrow::float64()),
        arrow::field("var_ms2_mi_score", arrow::float64()),
        arrow::field("var_ms2_mi_weighted_score", arrow::float64()),
        arrow::field("var_ms2_mi_ratio_score", arrow::float64()),
        arrow::field("var_ms2_norm_rt_score", arrow::float64()),
        arrow::field("var_ms2_xcorr_coelution", arrow::float64()),
        arrow::field("var_ms2_xcorr_coelution_weighted", arrow::float64()),
        arrow::field("var_ms2_xcorr_shape", arrow::float64()),
        arrow::field("var_ms2_xcorr_shape_weighted", arrow::float64()),
        arrow::field("var_ms2_yseries_score", arrow::float64()),
        arrow::field("var_ms2_elution_model_fit_score", arrow::float64()),
        arrow::field("var_ms2_im_xcorr_shape", arrow::float64()),
        arrow::field("var_ms2_im_xcorr_coelution", arrow::float64()),
        arrow::field("var_ms2_im_delta_score", arrow::float64()),
        arrow::field("var_ms2_im_log_intensity", arrow::float64())
      });
      auto features_table = arrow::Table::Make(features_schema, {
        ParquetFile::finishArray(feature_id_builder, "feature_id"),
        ParquetFile::finishArray(feature_run_id_builder, "run_id"),
        ParquetFile::finishArray(precursor_id_builder, "precursor_id"),
        ParquetFile::finishArray(exp_rt_builder, "exp_rt"),
        ParquetFile::finishArray(exp_im_builder, "exp_im"),
        ParquetFile::finishArray(norm_rt_builder, "norm_rt"),
        ParquetFile::finishArray(delta_rt_builder, "delta_rt"),
        ParquetFile::finishArray(left_width_builder, "left_width"),
        ParquetFile::finishArray(right_width_builder, "right_width"),
        ParquetFile::finishArray(exp_im_left_builder, "exp_im_leftwidth"),
        ParquetFile::finishArray(exp_im_right_builder, "exp_im_rightwidth"),
        ParquetFile::finishArray(ms1_area_builder, "ms1_area_intensity"),
        ParquetFile::finishArray(ms1_apex_builder, "ms1_apex_intensity"),
        ParquetFile::finishArray(ms1_exp_im_builder, "ms1_exp_im"),
        ParquetFile::finishArray(ms1_delta_im_builder, "ms1_delta_im"),
        ParquetFile::finishArray(var_ms1_massdev_builder, "var_ms1_massdev_score"),
        ParquetFile::finishArray(var_ms1_im_ms1_delta_builder, "var_ms1_im_ms1_delta_score"),
        ParquetFile::finishArray(var_ms1_mi_builder, "var_ms1_mi_score"),
        ParquetFile::finishArray(var_ms1_mi_contrast_builder, "var_ms1_mi_contrast_score"),
        ParquetFile::finishArray(var_ms1_mi_combined_builder, "var_ms1_mi_combined_score"),
        ParquetFile::finishArray(var_ms1_iso_corr_builder, "var_ms1_isotope_correlation_score"),
        ParquetFile::finishArray(var_ms1_iso_overlap_builder, "var_ms1_isotope_overlap_score"),
        ParquetFile::finishArray(var_ms1_xcorr_coelution_builder, "var_ms1_xcorr_coelution"),
        ParquetFile::finishArray(var_ms1_xcorr_coelution_contrast_builder, "var_ms1_xcorr_coelution_contrast"),
        ParquetFile::finishArray(var_ms1_xcorr_coelution_combined_builder, "var_ms1_xcorr_coelution_combined"),
        ParquetFile::finishArray(var_ms1_xcorr_shape_builder, "var_ms1_xcorr_shape"),
        ParquetFile::finishArray(var_ms1_xcorr_shape_contrast_builder, "var_ms1_xcorr_shape_contrast"),
        ParquetFile::finishArray(var_ms1_xcorr_shape_combined_builder, "var_ms1_xcorr_shape_combined"),
        ParquetFile::finishArray(ms2_area_builder, "ms2_area_intensity"),
        ParquetFile::finishArray(ms2_total_area_builder, "ms2_total_area_intensity"),
        ParquetFile::finishArray(ms2_apex_builder, "ms2_apex_intensity"),
        ParquetFile::finishArray(ms2_exp_im_builder, "ms2_exp_im"),
        ParquetFile::finishArray(ms2_exp_im_left_builder, "ms2_exp_im_leftwidth"),
        ParquetFile::finishArray(ms2_exp_im_right_builder, "ms2_exp_im_rightwidth"),
        ParquetFile::finishArray(ms2_delta_im_builder, "ms2_delta_im"),
        ParquetFile::finishArray(ms2_total_mi_builder, "ms2_total_mi"),
        ParquetFile::finishArray(var_ms2_bseries_builder, "var_ms2_bseries_score"),
        ParquetFile::finishArray(var_ms2_dotprod_builder, "var_ms2_dotprod_score"),
        ParquetFile::finishArray(var_ms2_intensity_builder, "var_ms2_intensity_score"),
        ParquetFile::finishArray(var_ms2_iso_corr_builder, "var_ms2_isotope_correlation_score"),
        ParquetFile::finishArray(var_ms2_iso_overlap_builder, "var_ms2_isotope_overlap_score"),
        ParquetFile::finishArray(var_ms2_library_corr_builder, "var_ms2_library_corr"),
        ParquetFile::finishArray(var_ms2_library_dotprod_builder, "var_ms2_library_dotprod"),
        ParquetFile::finishArray(var_ms2_library_manhattan_builder, "var_ms2_library_manhattan"),
        ParquetFile::finishArray(var_ms2_library_rmsd_builder, "var_ms2_library_rmsd"),
        ParquetFile::finishArray(var_ms2_library_rootmeansquare_builder, "var_ms2_library_rootmeansquare"),
        ParquetFile::finishArray(var_ms2_library_sangle_builder, "var_ms2_library_sangle"),
        ParquetFile::finishArray(var_ms2_log_sn_builder, "var_ms2_log_sn_score"),
        ParquetFile::finishArray(var_ms2_manhattan_builder, "var_ms2_manhattan_score"),
        ParquetFile::finishArray(var_ms2_massdev_builder, "var_ms2_massdev_score"),
        ParquetFile::finishArray(var_ms2_massdev_weighted_builder, "var_ms2_massdev_score_weighted"),
        ParquetFile::finishArray(var_ms2_mi_builder, "var_ms2_mi_score"),
        ParquetFile::finishArray(var_ms2_mi_weighted_builder, "var_ms2_mi_weighted_score"),
        ParquetFile::finishArray(var_ms2_mi_ratio_builder, "var_ms2_mi_ratio_score"),
        ParquetFile::finishArray(var_ms2_norm_rt_builder, "var_ms2_norm_rt_score"),
        ParquetFile::finishArray(var_ms2_xcorr_coelution_builder, "var_ms2_xcorr_coelution"),
        ParquetFile::finishArray(var_ms2_xcorr_coelution_weighted_builder, "var_ms2_xcorr_coelution_weighted"),
        ParquetFile::finishArray(var_ms2_xcorr_shape_builder, "var_ms2_xcorr_shape"),
        ParquetFile::finishArray(var_ms2_xcorr_shape_weighted_builder, "var_ms2_xcorr_shape_weighted"),
        ParquetFile::finishArray(var_ms2_yseries_builder, "var_ms2_yseries_score"),
        ParquetFile::finishArray(var_ms2_elution_model_fit_builder, "var_ms2_elution_model_fit_score"),
        ParquetFile::finishArray(var_ms2_im_xcorr_shape_builder, "var_ms2_im_xcorr_shape"),
        ParquetFile::finishArray(var_ms2_im_xcorr_coelution_builder, "var_ms2_im_xcorr_coelution"),
        ParquetFile::finishArray(var_ms2_im_delta_builder, "var_ms2_im_delta_score"),
        ParquetFile::finishArray(var_ms2_im_log_intensity_builder, "var_ms2_im_log_intensity")
      });
      write_tasks.emplace_back(std::async(std::launch::async, [features_table, run_path]() {
        ParquetFile::writeTable(features_table, run_path + "/features.parquet");
      }));
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
        ParquetFile::finishArray(ft_feature_id_builder, "feature_id"),
        ParquetFile::finishArray(ft_run_id_builder, "run_id"),
        ParquetFile::finishArray(ft_transition_id_builder, "transition_id"),
        ParquetFile::finishArray(ft_area_builder, "area_intensity"),
        ParquetFile::finishArray(ft_total_area_builder, "total_area_intensity"),
        ParquetFile::finishArray(ft_apex_int_builder, "apex_intensity"),
        ParquetFile::finishArray(ft_apex_rt_builder, "apex_rt"),
        ParquetFile::finishArray(ft_rt_fwhm_builder, "rt_fwhm"),
        ParquetFile::finishArray(ft_masserror_builder, "masserror_ppm"),
        ParquetFile::finishArray(ft_total_mi_builder, "total_mi"),
        ParquetFile::finishArray(ft_var_intensity_builder, "var_intensity_score"),
        ParquetFile::finishArray(ft_var_intensity_ratio_builder, "var_intensity_ratio_score"),
        ParquetFile::finishArray(ft_var_log_intensity_builder, "var_log_intensity"),
        ParquetFile::finishArray(ft_var_xcorr_coelution_builder, "var_xcorr_coelution"),
        ParquetFile::finishArray(ft_var_xcorr_shape_builder, "var_xcorr_shape"),
        ParquetFile::finishArray(ft_var_log_sn_builder, "var_log_sn_score"),
        ParquetFile::finishArray(ft_var_massdev_builder, "var_massdev_score"),
        ParquetFile::finishArray(ft_var_mi_builder, "var_mi_score"),
        ParquetFile::finishArray(ft_var_mi_ratio_builder, "var_mi_ratio_score"),
        ParquetFile::finishArray(ft_var_isotope_corr_builder, "var_isotope_correlation_score"),
        ParquetFile::finishArray(ft_var_isotope_overlap_builder, "var_isotope_overlap_score"),
        ParquetFile::finishArray(ft_exp_im_builder, "exp_im"),
        ParquetFile::finishArray(ft_exp_im_left_builder, "exp_im_leftwidth"),
        ParquetFile::finishArray(ft_exp_im_right_builder, "exp_im_rightwidth"),
        ParquetFile::finishArray(ft_delta_im_builder, "delta_im"),
        ParquetFile::finishArray(ft_var_im_delta_builder, "var_im_delta_score"),
        ParquetFile::finishArray(ft_var_im_log_intensity_builder, "var_im_log_intensity"),
        ParquetFile::finishArray(ft_var_im_xcorr_coelution_contrast_builder, "var_im_xcorr_coelution_contrast"),
        ParquetFile::finishArray(ft_var_im_xcorr_shape_contrast_builder, "var_im_xcorr_shape_contrast"),
        ParquetFile::finishArray(ft_var_im_xcorr_coelution_combined_builder, "var_im_xcorr_coelution_combined"),
        ParquetFile::finishArray(ft_var_im_xcorr_shape_combined_builder, "var_im_xcorr_shape_combined"),
        ParquetFile::finishArray(ft_start_position_at_5_builder, "start_position_at_5"),
        ParquetFile::finishArray(ft_end_position_at_5_builder, "end_position_at_5"),
        ParquetFile::finishArray(ft_start_position_at_10_builder, "start_position_at_10"),
        ParquetFile::finishArray(ft_end_position_at_10_builder, "end_position_at_10"),
        ParquetFile::finishArray(ft_start_position_at_50_builder, "start_position_at_50"),
        ParquetFile::finishArray(ft_end_position_at_50_builder, "end_position_at_50"),
        ParquetFile::finishArray(ft_total_width_builder, "total_width"),
        ParquetFile::finishArray(ft_tailing_factor_builder, "tailing_factor"),
        ParquetFile::finishArray(ft_asymmetry_factor_builder, "asymmetry_factor"),
        ParquetFile::finishArray(ft_slope_of_baseline_builder, "slope_of_baseline"),
        ParquetFile::finishArray(ft_baseline_delta_2_height_builder, "baseline_delta_2_height"),
        ParquetFile::finishArray(ft_points_across_baseline_builder, "points_across_baseline"),
        ParquetFile::finishArray(ft_points_across_half_height_builder, "points_across_half_height")
      });
      write_tasks.emplace_back(std::async(std::launch::async, [ft_table, run_path]() {
        ParquetFile::writeTable(ft_table, run_path + "/feature_transition.parquet");
      }));
    }

    // Always write a feature_precursor.parquet file for this run. Consumers
    // expect a consistent per-run layout; writing an empty table when no
    // precursor rows are present prevents downstream breakage.
    {
      auto fp_schema = arrow::schema({
        arrow::field("feature_id", arrow::int64()),
        arrow::field("run_id", arrow::int64()),
        arrow::field("precursor_isotope", arrow::int32()),
        arrow::field("precursor_area_intensity", arrow::float64()),
        arrow::field("precursor_apex_intensity", arrow::float64())
      });
      auto fp_table = arrow::Table::Make(fp_schema, {
        ParquetFile::finishArray(fp_feature_id_builder, "feature_id"),
        ParquetFile::finishArray(fp_run_id_builder, "run_id"),
        ParquetFile::finishArray(fp_isotope_builder, "precursor_isotope"),
        ParquetFile::finishArray(fp_area_builder, "precursor_area_intensity"),
        ParquetFile::finishArray(fp_apex_builder, "precursor_apex_intensity")
      });
      write_tasks.emplace_back(std::async(std::launch::async, [fp_table, run_path]() {
        ParquetFile::writeTable(fp_table, run_path + "/feature_precursor.parquet");
      }));
    }

    try
    {
      waitForWriteTasks_(write_tasks);
    }
    catch (...)
    {
      // If any asynchronous write failed, remove the partially-created run
      // directory to avoid leaving a stale run_id=<id> that blocks retries.
      if (File::exists(run_path))
      {
        File::removeDirRecursively(run_path);
      }
      throw;
    }

      // Now that all per-run files have been written (or their async tasks
      // have thrown), update the global runs.parquet file. Doing this after
      // successful writes prevents adding a runs entry that points to a
      // partially-written run directory.
      {
        arrow::Int64Builder run_id_builder;
        arrow::StringBuilder filename_builder;
        ParquetFile::appendOrThrow(run_id_builder.Append(run_id_clean), "run_id");
        ParquetFile::appendOrThrow(filename_builder.Append(std::string(input_filename)), "filename");

        auto runs_schema = arrow::schema({
          arrow::field("run_id", arrow::int64()),
          arrow::field("filename", arrow::utf8())
        });
        auto runs_table = arrow::Table::Make(runs_schema, {
          ParquetFile::finishArray(run_id_builder, "run_id"),
          ParquetFile::finishArray(filename_builder, "filename")
        });
        const String runs_parquet = runs_dir + "/runs.parquet";
        if (File::exists(runs_parquet))
        {
          auto existing_table = ParquetFile::readTable(runs_parquet);
          auto existing_run_ids = std::static_pointer_cast<arrow::Int64Array>(ParquetFile::getColumn(existing_table, "run_id"));
          const int64_t existing_rows = existing_table->num_rows();
          for (int64_t row = 0; row < existing_rows; ++row)
          {
            if (existing_run_ids->Value(row) == static_cast<int64_t>(run_id_clean))
            {
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Run id already present in runs.parquet", String(run_id_clean));
            }
          }
          auto combined_result = arrow::ConcatenateTables({existing_table, runs_table});
          if (!combined_result.ok())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Failed to append runs table", combined_result.status().ToString());
          }
          ParquetFile::writeTable(combined_result.ValueOrDie(), runs_parquet);
        }
        else
        {
          ParquetFile::writeTable(runs_table, runs_parquet);
        }
      }
      const String runs_parquet = runs_dir + "/runs.parquet";
    auto runs = readRuns_(runs_parquet);
    std::vector<RunCounts> run_counts;
    run_counts.reserve(runs.size());

    RunCounts total_counts;
    for (const auto& run : runs)
    {
      const String current_run_path = runs_dir + "/run_id=" + String(run.run_id);
      RunCounts counts;
      counts.features = getParquetRowCount_(current_run_path + "/features.parquet");
      counts.feature_precursor = getParquetRowCount_(current_run_path + "/feature_precursor.parquet");
      counts.feature_transition = getParquetRowCount_(current_run_path + "/feature_transition.parquet");
      run_counts.push_back(counts);

      total_counts.features += counts.features;
      total_counts.feature_precursor += counts.feature_precursor;
      total_counts.feature_transition += counts.feature_transition;
    }

    writeMetadata_(base_dir, runs, run_counts, total_counts);
#endif
  }

} // namespace OpenMS
