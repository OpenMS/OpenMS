// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h>

#include "OpenSwathCanonicalLibraryMappingHelper.h"
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionParquetFile.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/FORMAT/SqliteConnector_impl.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

#include <fstream>
#include <cmath>
#include <future>
#include <limits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <filesystem>

#ifdef signals
#undef signals
#endif
#include <arrow/api.h>
#include <parquet/file_reader.h>

namespace OpenMS
{
  namespace
  {
    using OpenMS::Size;

    bool parseOptionalInt64_(const std::string& text, int64_t& value)
    {
      try
      {
        value = StringUtils::toInt64(text);
        return true;
      }
      catch (Exception::ConversionError&)
      {
        return false;
      }
    }

    bool parseOptionalDouble_(const std::string& text, double& value)
    {
      try
      {
        value = StringUtils::toDouble(text);
        return true;
      }
      catch (Exception::ConversionError&)
      {
        return false;
      }
    }

    std::string oswValueToDebugString_(const OpenSwathOSWWriter::OSWValue& value)
    {
      switch (value.type)
      {
        case OpenSwathOSWWriter::OSWValue::Type::Null:
          return "NULL";
        case OpenSwathOSWWriter::OSWValue::Type::Int64:
          return StringUtils::toStr(value.asInt());
        case OpenSwathOSWWriter::OSWValue::Type::Double:
          return StringUtils::toStr(value.asDouble());
        case OpenSwathOSWWriter::OSWValue::Type::Text:
          return value.asText();
      }
      return "NULL";
    }

    bool oswValueToInt64_(const OpenSwathOSWWriter::OSWValue& value, int64_t& result)
    {
      if (value.isNull())
      {
        return false;
      }

      switch (value.type)
      {
        case OpenSwathOSWWriter::OSWValue::Type::Int64:
          result = value.asInt();
          return true;
        case OpenSwathOSWWriter::OSWValue::Type::Double:
        {
          const double double_value = value.asDouble();
          if (!std::isfinite(double_value))
          {
            return false;
          }
          const double rounded = std::round(double_value);
          if (std::fabs(double_value - rounded) > 1e-9)
          {
            return false;
          }
          if (rounded < static_cast<double>(std::numeric_limits<int64_t>::min()) ||
              rounded > static_cast<double>(std::numeric_limits<int64_t>::max()))
          {
            return false;
          }
          result = static_cast<int64_t>(rounded);
          return true;
        }
        case OpenSwathOSWWriter::OSWValue::Type::Text:
          return parseOptionalInt64_(value.asText(), result);
        case OpenSwathOSWWriter::OSWValue::Type::Null:
          return false;
      }

      return false;
    }

    bool oswValueToDouble_(const OpenSwathOSWWriter::OSWValue& value, double& result)
    {
      if (value.isNull())
      {
        return false;
      }

      switch (value.type)
      {
        case OpenSwathOSWWriter::OSWValue::Type::Int64:
          result = static_cast<double>(value.asInt());
          return true;
        case OpenSwathOSWWriter::OSWValue::Type::Double:
          result = value.asDouble();
          return std::isfinite(result);
        case OpenSwathOSWWriter::OSWValue::Type::Text:
          return parseOptionalDouble_(value.asText(), result) && std::isfinite(result);
        case OpenSwathOSWWriter::OSWValue::Type::Null:
          return false;
      }

      return false;
    }

    int64_t requireOSWInt64_(const OpenSwathOSWWriter::OSWValue& value, const char* column)
    {
      int64_t result = 0;
      if (oswValueToInt64_(value, result))
      {
        return result;
      }

      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Expected integer OSW value for column '" + std::string(column) + "'",
                                    oswValueToDebugString_(value));
    }

    void appendOSWInt64_(arrow::Int64Builder& builder,
                         const OpenSwathOSWWriter::OSWValue& value,
                         const char* column)
    {
      if (value.isNull())
      {
        ParquetFile::appendOrThrow(builder.AppendNull(), column);
        return;
      }

      int64_t int_value = 0;
      if (!oswValueToInt64_(value, int_value))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Expected integer OSW value for parquet column '" + std::string(column) + "'",
                                      oswValueToDebugString_(value));
      }
      ParquetFile::appendOrThrow(builder.Append(int_value), column);
    }

    void appendOSWInt32_(arrow::Int32Builder& builder,
                         const OpenSwathOSWWriter::OSWValue& value,
                         const char* column)
    {
      if (value.isNull())
      {
        ParquetFile::appendOrThrow(builder.AppendNull(), column);
        return;
      }

      int64_t int_value = 0;
      if (!oswValueToInt64_(value, int_value) ||
          int_value < std::numeric_limits<int32_t>::min() ||
          int_value > std::numeric_limits<int32_t>::max())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Expected 32-bit integer OSW value for parquet column '" + std::string(column) + "'",
                                      oswValueToDebugString_(value));
      }
      ParquetFile::appendOrThrow(builder.Append(static_cast<int32_t>(int_value)), column);
    }

    void appendOSWDouble_(arrow::DoubleBuilder& builder,
                          const OpenSwathOSWWriter::OSWValue& value,
                          const char* column)
    {
      if (value.isNull())
      {
        ParquetFile::appendOrThrow(builder.AppendNull(), column);
        return;
      }

      double double_value = 0.0;
      if (!oswValueToDouble_(value, double_value))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Expected floating-point OSW value for parquet column '" + std::string(column) + "'",
                                      oswValueToDebugString_(value));
      }
      ParquetFile::appendOrThrow(builder.Append(double_value), column);
    }

    int64_t resolvePrecursorId_(const OpenSwathOSWWriter::OSWValue& value,
                                const std::unordered_map<std::string, int64_t>& compound_to_precursor,
                                const std::unordered_set<int64_t>& valid_precursor_ids)
    {
      const std::string precursor_ref = oswValueToDebugString_(value);
      auto precursor_it = compound_to_precursor.find(precursor_ref);
      if (precursor_it != compound_to_precursor.end())
      {
        return precursor_it->second;
      }

      int64_t precursor_id = 0;
      if (oswValueToInt64_(value, precursor_id) && valid_precursor_ids.contains(precursor_id))
      {
        return precursor_id;
      }

      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "OSW feature row references unknown precursor",
                                    precursor_ref);
    }

    int64_t resolveTransitionId_(const OpenSwathOSWWriter::OSWValue& value,
                                 const std::unordered_map<std::string, int64_t>& transition_to_id,
                                 int64_t next_transition_id)
    {
      const std::string transition_ref = oswValueToDebugString_(value);
      auto transition_it = transition_to_id.find(transition_ref);
      if (transition_it != transition_to_id.end())
      {
        return transition_it->second;
      }

      int64_t transition_id = 0;
      if (oswValueToInt64_(value, transition_id))
      {
        if (transition_id >= 0 && transition_id < next_transition_id)
        {
          return transition_id;
        }
      }

      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "OSW transition row references unknown transition",
                                    transition_ref);
    }

    OpenSwathOSWWriter::OSWData canonicalRowsFromFeatureMap_(const FeatureMap& feature_map,
                                                             UInt64 run_id,
                                                             bool enable_uis_scoring)
    {
      OpenSwathOSWWriter row_writer("", enable_uis_scoring);
      row_writer.setRunId(run_id);

      OpenSwathOSWWriter::OSWData osw_output;
      Size transition_row_estimate = 0;
      for (const auto& feature : feature_map)
      {
        transition_row_estimate += feature.getSubordinates().size();
      }
      osw_output.reserve(feature_map.size(), transition_row_estimate);

      std::map<std::string, FeatureMap> grouped_features;
      for (const auto& feature : feature_map)
      {
        if (!feature.metaValueExists("PeptideRef"))
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Feature missing PeptideRef meta value required for canonical OSW row generation.");
        }
        const std::string peptide_ref = StringUtils::toStr(feature.getMetaValue("PeptideRef"));
        grouped_features[peptide_ref].push_back(feature);
      }

      for (const auto& grouped_entry : grouped_features)
      {
        row_writer.prepareRowsInto(osw_output, OpenSwath::LightCompound(), nullptr, grouped_entry.second, grouped_entry.first);
      }

      return osw_output;
    }

    struct RunEntry
    {
      int64_t run_id = 0;
      std::string filename;
    };

    struct RunCounts
    {
      int64_t features = 0;
      int64_t feature_precursor = 0;
      int64_t feature_transition = 0;
    };

    int64_t getParquetRowCount_(const std::string& filename)
    {
      return OpenMS::ParquetFile::rowCount(filename);
    }

    std::vector<RunEntry> readRuns_(const std::string& runs_parquet)
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

    std::string jsonEscape_(const std::string& input)
    {
      return OpenMS::ParquetFile::jsonEscape(input);
    }

    void writeMetadata_(const std::string& base_dir,
                        const std::vector<RunEntry>& runs,
                        const std::vector<RunCounts>& run_counts,
                        const RunCounts& total_counts)
    {
      const std::string metadata_path = base_dir + "/metadata.json";
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
                                      "Run count metadata mismatch",StringUtils::toStr(run_counts.size()));
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

  } // namespace

  void OpenSwathOSWParquetWriter::write(const std::string& output_path,
                                        const OpenSwath::LightTargetedExperiment& assay_library,
                                        const FeatureMap& feature_map,
                                        UInt64 run_id,
                                        const std::string& input_filename,
                                        bool enable_uis_scoring) const
  {
    const auto osw_output = canonicalRowsFromFeatureMap_(feature_map, run_id, enable_uis_scoring);
    write(output_path, assay_library, osw_output, run_id, input_filename, enable_uis_scoring);
  }

  void OpenSwathOSWParquetWriter::write(const std::string& output_path,
                                        const OpenSwath::LightTargetedExperiment& assay_library,
                                        const OpenSwathOSWWriter::OSWData& osw_output,
                                        UInt64 run_id,
                                        const std::string& input_filename,
                                        bool /*enable_uis_scoring*/) const
  {
    const UInt64 run_id_clean = Internal::SqliteHelper::clearSignBit(run_id);
    const bool output_is_dir = File::isDirectory(output_path);
    std::unique_ptr<File::TempDir> temp_dir;
    std::string base_dir = output_path;
    if (!output_is_dir)
    {
      if (File::exists(output_path) && preserve_existing_)
      {
        // Preserve prior runs/library by unpacking existing archive first.
          base_dir = ZipArchiveFile::unzipDirectory(output_path, temp_dir);
        // Inform the user that existing runs were found and will be preserved
        OPENMS_LOG_WARN << "Existing archive '" << File::absolutePath(output_path)
                        << "' was detected — its contents will be preserved and new runs will be appended.\n"
                        << "To overwrite the archive remove it before running the writer." << std::endl;
      }
      else
      {
        temp_dir = std::make_unique<File::TempDir>();
        base_dir = temp_dir->getPath() + "/oswpq_output";
      }
    }

    const std::string library_dir = base_dir + "/library";
    const std::string runs_dir = base_dir + "/runs";
    if (!File::exists(base_dir))
    {
      File::makeDir(base_dir);
    }
    if (!File::exists(runs_dir))
    {
      File::makeDir(runs_dir);
    }

    const std::string precursors_path = library_dir + "/precursors.parquet";
    const std::string transitions_path = library_dir + "/transitions.parquet";
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
      // Use deduplicated counts from the provided assay_library (unique compound ids
      // and unique transition names) to compare against existing parquet tables. The
      // writer deduplicates compounds/transitions when generating the library, so a
      // direct comparison to raw vector sizes can give false incompatibility errors.
      std::unordered_set<std::string> unique_compounds;
      unique_compounds.reserve(assay_library.compounds.size());
      for (const auto& c : assay_library.compounds) unique_compounds.insert(c.id);
      const int64_t expected_precursors = static_cast<int64_t>(unique_compounds.size());

      std::unordered_set<std::string> unique_transitions;
      unique_transitions.reserve(assay_library.transitions.size());
      for (const auto& t : assay_library.transitions) unique_transitions.insert(t.transition_name);
      const int64_t expected_transitions = static_cast<int64_t>(unique_transitions.size());

      if (existing_precursors != expected_precursors || existing_transitions != expected_transitions)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Existing library at '" + library_dir + "' appears incompatible with provided assay_library. Please rebuild the library or use a different output path.",
                                      "existing_precursors=" + StringUtils::toStr(existing_precursors) + ", expected_precursors=" + StringUtils::toStr(expected_precursors));
      }
    }
    if (!library_ready)
    {
      const std::string library_tmp_dir = base_dir + "/library_tmp";
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
      File::copyDirRecursively(library_tmp_dir + "/library", library_dir);
      File::removeDirRecursively(library_tmp_dir);
    }

    const std::string run_path = runs_dir + "/run_id=" + StringUtils::toStr(run_id_clean);
    if (File::exists(run_path))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Run output directory already exists", run_path);
    }
    File::makeDir(run_path);

    const auto canonical_mapping = Internal::buildOpenSwathCanonicalLibraryMapping(assay_library);
    const auto& compound_to_precursor = canonical_mapping.compound_to_precursor;
    const auto& transition_to_id = canonical_mapping.transition_to_id;
    std::unordered_set<int64_t> valid_precursor_ids;
    valid_precursor_ids.reserve(compound_to_precursor.size());
    for (const auto& [_, precursor_id] : compound_to_precursor)
    {
      valid_precursor_ids.insert(precursor_id);
    }
    const int64_t transition_id_count = static_cast<int64_t>(assay_library.transitions.size());

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

    // Reserve Arrow builders where possible to avoid repeated reallocations.
    const int64_t n_features = static_cast<int64_t>(osw_output.feature_rows.size());
    ParquetFile::appendOrThrow(feature_id_builder.Reserve(n_features), "feature_id (reserve)");
    ParquetFile::appendOrThrow(feature_run_id_builder.Reserve(n_features), "run_id (reserve)");
    ParquetFile::appendOrThrow(precursor_id_builder.Reserve(n_features), "precursor_id (reserve)");
    ParquetFile::appendOrThrow(exp_rt_builder.Reserve(n_features), "exp_rt (reserve)");

    const int64_t approx_ft = std::max<int64_t>(1, static_cast<int64_t>(osw_output.feature_transition_rows.size()));
    ParquetFile::appendOrThrow(ft_feature_id_builder.Reserve(approx_ft), "ft.feature_id (reserve)");
    ParquetFile::appendOrThrow(ft_run_id_builder.Reserve(approx_ft), "ft.run_id (reserve)");
    ParquetFile::appendOrThrow(ft_transition_id_builder.Reserve(approx_ft), "ft.transition_id (reserve)");

    const int64_t approx_fp = std::max<int64_t>(1, static_cast<int64_t>(osw_output.feature_precursor_rows.size()));
    ParquetFile::appendOrThrow(fp_feature_id_builder.Reserve(approx_fp), "fp.feature_id (reserve)");
    ParquetFile::appendOrThrow(fp_run_id_builder.Reserve(approx_fp), "fp.run_id (reserve)");

    std::unordered_map<int64_t, const OpenSwathOSWWriter::FeatureMS1Row*> feature_to_ms1;
    feature_to_ms1.reserve(osw_output.feature_ms1_rows.size());
    for (const auto& ms1_row : osw_output.feature_ms1_rows)
    {
      feature_to_ms1[requireOSWInt64_(ms1_row[0], "FEATURE_MS1.FEATURE_ID")] = &ms1_row;
    }

    std::unordered_map<int64_t, const OpenSwathOSWWriter::FeatureMS2Row*> feature_to_ms2;
    feature_to_ms2.reserve(osw_output.feature_ms2_rows.size());
    for (const auto& ms2_row : osw_output.feature_ms2_rows)
    {
      feature_to_ms2[requireOSWInt64_(ms2_row[0], "FEATURE_MS2.FEATURE_ID")] = &ms2_row;
    }

    const std::vector<arrow::DoubleBuilder*> feature_numeric_builders = {
      &exp_rt_builder, &exp_im_builder, &norm_rt_builder, &delta_rt_builder,
      &left_width_builder, &right_width_builder, &exp_im_left_builder, &exp_im_right_builder
    };
    const std::vector<const char*> feature_numeric_names = {
      "exp_rt", "exp_im", "norm_rt", "delta_rt",
      "left_width", "right_width", "exp_im_leftwidth", "exp_im_rightwidth"
    };

    const std::vector<arrow::DoubleBuilder*> ms1_builders = {
      &ms1_area_builder, &ms1_apex_builder, &ms1_exp_im_builder, &ms1_delta_im_builder,
      &var_ms1_massdev_builder, &var_ms1_im_ms1_delta_builder, &var_ms1_mi_builder,
      &var_ms1_mi_contrast_builder, &var_ms1_mi_combined_builder, &var_ms1_iso_corr_builder,
      &var_ms1_iso_overlap_builder, &var_ms1_xcorr_coelution_builder,
      &var_ms1_xcorr_coelution_contrast_builder, &var_ms1_xcorr_coelution_combined_builder,
      &var_ms1_xcorr_shape_builder, &var_ms1_xcorr_shape_contrast_builder,
      &var_ms1_xcorr_shape_combined_builder
    };
    const std::vector<const char*> ms1_names = {
      "ms1_area_intensity", "ms1_apex_intensity", "ms1_exp_im", "ms1_delta_im",
      "var_ms1_massdev_score", "var_ms1_im_ms1_delta_score", "var_ms1_mi_score",
      "var_ms1_mi_contrast_score", "var_ms1_mi_combined_score", "var_ms1_isotope_correlation_score",
      "var_ms1_isotope_overlap_score", "var_ms1_xcorr_coelution",
      "var_ms1_xcorr_coelution_contrast", "var_ms1_xcorr_coelution_combined",
      "var_ms1_xcorr_shape", "var_ms1_xcorr_shape_contrast", "var_ms1_xcorr_shape_combined"
    };

    const std::vector<arrow::DoubleBuilder*> ms2_builders = {
      &ms2_area_builder, &ms2_total_area_builder, &ms2_apex_builder, &ms2_exp_im_builder,
      &ms2_exp_im_left_builder, &ms2_exp_im_right_builder, &ms2_delta_im_builder,
      &ms2_total_mi_builder, &var_ms2_bseries_builder, &var_ms2_dotprod_builder,
      &var_ms2_intensity_builder, &var_ms2_iso_corr_builder, &var_ms2_iso_overlap_builder,
      &var_ms2_library_corr_builder, &var_ms2_library_dotprod_builder, &var_ms2_library_manhattan_builder,
      &var_ms2_library_rmsd_builder, &var_ms2_library_rootmeansquare_builder, &var_ms2_library_sangle_builder,
      &var_ms2_log_sn_builder, &var_ms2_manhattan_builder, &var_ms2_massdev_builder,
      &var_ms2_massdev_weighted_builder, &var_ms2_mi_builder, &var_ms2_mi_weighted_builder,
      &var_ms2_mi_ratio_builder, &var_ms2_norm_rt_builder, &var_ms2_xcorr_coelution_builder,
      &var_ms2_xcorr_coelution_weighted_builder, &var_ms2_xcorr_shape_builder,
      &var_ms2_xcorr_shape_weighted_builder, &var_ms2_yseries_builder,
      &var_ms2_elution_model_fit_builder, &var_ms2_im_xcorr_shape_builder,
      &var_ms2_im_xcorr_coelution_builder, &var_ms2_im_delta_builder,
      &var_ms2_im_log_intensity_builder
    };
    const std::vector<const char*> ms2_names = {
      "ms2_area_intensity", "ms2_total_area_intensity", "ms2_apex_intensity", "ms2_exp_im",
      "ms2_exp_im_leftwidth", "ms2_exp_im_rightwidth", "ms2_delta_im", "ms2_total_mi",
      "var_ms2_bseries_score", "var_ms2_dotprod_score", "var_ms2_intensity_score",
      "var_ms2_isotope_correlation_score", "var_ms2_isotope_overlap_score", "var_ms2_library_corr",
      "var_ms2_library_dotprod", "var_ms2_library_manhattan", "var_ms2_library_rmsd",
      "var_ms2_library_rootmeansquare", "var_ms2_library_sangle", "var_ms2_log_sn_score",
      "var_ms2_manhattan_score", "var_ms2_massdev_score", "var_ms2_massdev_score_weighted",
      "var_ms2_mi_score", "var_ms2_mi_weighted_score", "var_ms2_mi_ratio_score",
      "var_ms2_norm_rt_score", "var_ms2_xcorr_coelution", "var_ms2_xcorr_coelution_weighted",
      "var_ms2_xcorr_shape", "var_ms2_xcorr_shape_weighted", "var_ms2_yseries_score",
      "var_ms2_elution_model_fit_score", "var_ms2_im_xcorr_shape", "var_ms2_im_xcorr_coelution",
      "var_ms2_im_delta_score", "var_ms2_im_log_intensity"
    };

    const std::vector<arrow::DoubleBuilder*> transition_numeric_builders = {
      &ft_area_builder, &ft_total_area_builder, &ft_apex_int_builder, &ft_apex_rt_builder,
      &ft_rt_fwhm_builder, &ft_masserror_builder, &ft_total_mi_builder, &ft_var_intensity_builder,
      &ft_var_intensity_ratio_builder, &ft_var_log_intensity_builder, &ft_var_xcorr_coelution_builder,
      &ft_var_xcorr_shape_builder, &ft_var_log_sn_builder, &ft_var_massdev_builder,
      &ft_var_mi_builder, &ft_var_mi_ratio_builder, &ft_var_isotope_corr_builder,
      &ft_var_isotope_overlap_builder, &ft_exp_im_builder, &ft_exp_im_left_builder,
      &ft_exp_im_right_builder, &ft_delta_im_builder, &ft_var_im_delta_builder,
      &ft_var_im_log_intensity_builder, &ft_var_im_xcorr_coelution_contrast_builder,
      &ft_var_im_xcorr_shape_contrast_builder, &ft_var_im_xcorr_coelution_combined_builder,
      &ft_var_im_xcorr_shape_combined_builder, &ft_start_position_at_5_builder,
      &ft_end_position_at_5_builder, &ft_start_position_at_10_builder, &ft_end_position_at_10_builder,
      &ft_start_position_at_50_builder, &ft_end_position_at_50_builder, &ft_total_width_builder,
      &ft_tailing_factor_builder, &ft_asymmetry_factor_builder, &ft_slope_of_baseline_builder,
      &ft_baseline_delta_2_height_builder, &ft_points_across_baseline_builder,
      &ft_points_across_half_height_builder
    };
    const std::vector<const char*> transition_numeric_names = {
      "area_intensity", "total_area_intensity", "apex_intensity", "apex_rt",
      "rt_fwhm", "masserror_ppm", "total_mi", "var_intensity_score",
      "var_intensity_ratio_score", "var_log_intensity", "var_xcorr_coelution",
      "var_xcorr_shape", "var_log_sn_score", "var_massdev_score",
      "var_mi_score", "var_mi_ratio_score", "var_isotope_correlation_score",
      "var_isotope_overlap_score", "exp_im", "exp_im_leftwidth",
      "exp_im_rightwidth", "delta_im", "var_im_delta_score",
      "var_im_log_intensity", "var_im_xcorr_coelution_contrast",
      "var_im_xcorr_shape_contrast", "var_im_xcorr_coelution_combined",
      "var_im_xcorr_shape_combined", "start_position_at_5",
      "end_position_at_5", "start_position_at_10", "end_position_at_10",
      "start_position_at_50", "end_position_at_50", "total_width",
      "tailing_factor", "asymmetry_factor", "slope_of_baseline",
      "baseline_delta_2_height", "points_across_baseline",
      "points_across_half_height"
    };

    for (const auto& feature_row : osw_output.feature_rows)
    {
      const int64_t feature_id = requireOSWInt64_(feature_row[0], "FEATURE.ID");
      const int64_t feature_run_id = requireOSWInt64_(feature_row[1], "FEATURE.RUN_ID");
      if (feature_run_id != static_cast<int64_t>(run_id_clean))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Canonical OSW feature row belongs to a different run than the requested parquet partition.",
                                      StringUtils::toStr(feature_run_id));
      }
      const int64_t precursor_id = resolvePrecursorId_(feature_row[2], compound_to_precursor, valid_precursor_ids);

      ParquetFile::appendOrThrow(feature_id_builder.Append(feature_id), "feature_id");
      ParquetFile::appendOrThrow(feature_run_id_builder.Append(feature_run_id), "run_id");
      ParquetFile::appendOrThrow(precursor_id_builder.Append(precursor_id), "precursor_id");

      for (Size i = 0; i < feature_numeric_builders.size(); ++i)
      {
        appendOSWDouble_(*feature_numeric_builders[i], feature_row[i + 3], feature_numeric_names[i]);
      }

      auto ms1_it = feature_to_ms1.find(feature_id);
      if (ms1_it != feature_to_ms1.end())
      {
        const auto& ms1_row = *ms1_it->second;
        for (Size i = 0; i < ms1_builders.size(); ++i)
        {
          appendOSWDouble_(*ms1_builders[i], ms1_row[i + 1], ms1_names[i]);
        }
      }
      else
      {
        for (Size i = 0; i < ms1_builders.size(); ++i)
        {
          ParquetFile::appendOrThrow(ms1_builders[i]->AppendNull(), ms1_names[i]);
        }
      }

      auto ms2_it = feature_to_ms2.find(feature_id);
      if (ms2_it != feature_to_ms2.end())
      {
        const auto& ms2_row = *ms2_it->second;
        for (Size i = 0; i < ms2_builders.size(); ++i)
        {
          appendOSWDouble_(*ms2_builders[i], ms2_row[i + 1], ms2_names[i]);
        }
      }
      else
      {
        for (Size i = 0; i < ms2_builders.size(); ++i)
        {
          ParquetFile::appendOrThrow(ms2_builders[i]->AppendNull(), ms2_names[i]);
        }
      }
    }

    for (const auto& transition_row : osw_output.feature_transition_rows)
    {
      const int64_t feature_id = requireOSWInt64_(transition_row[0], "FEATURE_TRANSITION.FEATURE_ID");
      const int64_t transition_id_value = resolveTransitionId_(transition_row[1], transition_to_id, transition_id_count);

      ParquetFile::appendOrThrow(ft_feature_id_builder.Append(feature_id), "feature_id");
      ParquetFile::appendOrThrow(ft_run_id_builder.Append(static_cast<int64_t>(run_id_clean)), "run_id");
      ParquetFile::appendOrThrow(ft_transition_id_builder.Append(transition_id_value), "transition_id");

      for (Size i = 0; i < transition_numeric_builders.size(); ++i)
      {
        appendOSWDouble_(*transition_numeric_builders[i], transition_row[i + 2], transition_numeric_names[i]);
      }
    }

    for (const auto& precursor_row : osw_output.feature_precursor_rows)
    {
      const int64_t feature_id = requireOSWInt64_(precursor_row[0], "FEATURE_PRECURSOR.FEATURE_ID");
      ParquetFile::appendOrThrow(fp_feature_id_builder.Append(feature_id), "feature_id");
      ParquetFile::appendOrThrow(fp_run_id_builder.Append(static_cast<int64_t>(run_id_clean)), "run_id");
      appendOSWInt32_(fp_isotope_builder, precursor_row[1], "precursor_isotope");
      appendOSWDouble_(fp_area_builder, precursor_row[2], "precursor_area_intensity");
      appendOSWDouble_(fp_apex_builder, precursor_row[3], "precursor_apex_intensity");
    }

    std::vector<std::future<void>> write_tasks;
    try
    {
    {
      auto features_schema = OSWFeatureSchema::schema();
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
      // Validate features table against registry schema
      auto feat_validation = ArrowSchemaValidation::validate(features_table, OSWFeatureSchema::schema());
      if (!feat_validation.valid)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Features table schema validation failed: " + feat_validation.toString(), "");
      }
      write_tasks.emplace_back(std::async(std::launch::async, [features_table, run_path]() {
        ParquetFile::writeTable(features_table, run_path + "/features.parquet");
      }));
    }

    {
      auto ft_schema = OSWFeatureTransitionSchema::schema();
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
      // Validate feature_transition table against registry schema
      auto ft_validation = ArrowSchemaValidation::validate(ft_table, OSWFeatureTransitionSchema::schema());
      if (!ft_validation.valid)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Feature transition table schema validation failed: " + ft_validation.toString(), "");
      }
      write_tasks.emplace_back(std::async(std::launch::async, [ft_table, run_path]() {
        ParquetFile::writeTable(ft_table, run_path + "/feature_transition.parquet");
      }));
    }

    // Always write a feature_precursor.parquet file for this run. Consumers
    // expect a consistent per-run layout; writing an empty table when no
    // precursor rows are present prevents downstream breakage.
    {
      auto fp_schema = OSWFeaturePrecursorSchema::schema();
      auto fp_table = arrow::Table::Make(fp_schema, {
        ParquetFile::finishArray(fp_feature_id_builder, "feature_id"),
        ParquetFile::finishArray(fp_run_id_builder, "run_id"),
        ParquetFile::finishArray(fp_isotope_builder, "precursor_isotope"),
        ParquetFile::finishArray(fp_area_builder, "precursor_area_intensity"),
        ParquetFile::finishArray(fp_apex_builder, "precursor_apex_intensity")
      });
      // Validate feature_precursor table against registry schema
      auto fp_validation = ArrowSchemaValidation::validate(fp_table, OSWFeaturePrecursorSchema::schema());
      if (!fp_validation.valid)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Feature precursor table schema validation failed: " + fp_validation.toString(), "");
      }
      write_tasks.emplace_back(std::async(std::launch::async, [fp_table, run_path]() {
        ParquetFile::writeTable(fp_table, run_path + "/feature_precursor.parquet");
      }));
    }

      waitForWriteTasks_(write_tasks);
    }
    catch (...)
    {
      // If any validation or asynchronous write failed, remove the partially-created run
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

        auto runs_schema = OSWRunSchema::schema();
        auto runs_table = arrow::Table::Make(runs_schema, {
          ParquetFile::finishArray(run_id_builder, "run_id"),
          ParquetFile::finishArray(filename_builder, "filename")
        });
        // Validate runs table against registry schema
        auto runs_validation = ArrowSchemaValidation::validate(runs_table, OSWRunSchema::schema());
        if (!runs_validation.valid)
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Runs table schema validation failed: " + runs_validation.toString(), "");
        }
        const std::string runs_parquet = runs_dir + "/runs.parquet";
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
                                            "Run id already present in runs.parquet",StringUtils::toStr(run_id_clean));
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
      const std::string runs_parquet = runs_dir + "/runs.parquet";
    auto runs = readRuns_(runs_parquet);
    std::vector<RunCounts> run_counts;
    run_counts.reserve(runs.size());

    RunCounts total_counts;
    for (const auto& run : runs)
    {
      const std::string current_run_path = runs_dir + "/run_id=" + StringUtils::toStr(run.run_id);
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
    // If the requested output was a single archive file (not a directory),
    // stream files into the archive and write a sidecar index to enable
    // robust random-access reads later.
    if (!output_is_dir)
    {
      const std::filesystem::path dirpath = std::filesystem::u8path(std::string(base_dir));
      const std::string output_zip_abs = File::absolutePath(output_path);
      // If we're preserving an existing archive (we unpacked it above), don't
      // remove it here. Otherwise remove any existing file to start fresh.
      if (File::exists(output_zip_abs) && !preserve_existing_)
      {
        File::remove(output_zip_abs);
      }
      for (auto it = std::filesystem::recursive_directory_iterator(dirpath); it != std::filesystem::recursive_directory_iterator(); ++it)
      {
        if (it->is_directory()) continue;
        const auto full = it->path();
        std::string rel = std::filesystem::relative(full, dirpath).generic_string();
        ZipArchiveFile::addOrReplaceFromFile(output_path,std::string(rel),std::string(full.string()));
      }
      ZipArchiveFile::writeSidecarIndex(output_zip_abs);
    }
  }

void OpenSwathOSWParquetWriter::setPreserveExisting(bool preserve)
{
  preserve_existing_ = preserve;
}

} // namespace OpenMS
