// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetReader.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Exception.h>

#ifdef WITH_PARQUET
#include <arrow/api.h>
#endif

#include <fstream>
#include <unordered_map>

using namespace OpenMS;

OpenSwathOSWParquetReader::OpenSwathOSWParquetReader(const String& oswpq_dir)
{
  load(oswpq_dir);
}

void OpenSwathOSWParquetReader::load(const String& oswpq_dir)
{
#ifdef WITH_PARQUET
  rows_.clear();
  // remember the provided path for later fetch calls
  oswpq_dir_ = oswpq_dir;
  std::unique_ptr<File::TempDir> temp_dir;
  const String base_dir = ParquetFile::unzipDirectory(oswpq_dir, temp_dir);

  const String library_dir = base_dir + "/library";
  const String precursors_path = library_dir + "/precursors.parquet";
  const String transitions_path = library_dir + "/transitions.parquet";

  if (!File::exists(precursors_path) || !File::exists(transitions_path))
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Missing library parquet files in '" + oswpq_dir + "'");
  }

  // Read precursors -> charge / decoy
  auto precursors_table = ParquetFile::readTable(precursors_path);
  auto precursor_id_col = ParquetFile::getColumn(precursors_table, "precursor_id");
  auto charge_col = ParquetFile::getColumn(precursors_table, "charge");
  auto decoy_col = ParquetFile::getOptionalColumn(precursors_table, "decoy");

  std::unordered_map<int64_t, std::pair<int, bool>> precursor_info;
  precursor_info.reserve(precursors_table->num_rows());
  for (int64_t row = 0; row < precursors_table->num_rows(); ++row)
  {
    const int64_t pid = ParquetFile::getInt64(precursor_id_col, row, 0, false);
    const int charge = static_cast<int>(ParquetFile::getInt64(charge_col, row, 0, false));
    const bool decoy = ParquetFile::getBool(decoy_col, row, false, true);
    precursor_info.emplace(pid, std::make_pair(charge, decoy));
  }

  // Read transitions to compute detecting transition counts per precursor
  auto transitions_table = ParquetFile::readTable(transitions_path);
  auto transition_precursor_col = ParquetFile::getColumn(transitions_table, "precursor_id");
  auto detecting_col = ParquetFile::getColumn(transitions_table, "detecting");

  std::unordered_map<int64_t, int64_t> transition_counts;
  transition_counts.reserve(transitions_table->num_rows());
  for (int64_t row = 0; row < transitions_table->num_rows(); ++row)
  {
    const int64_t pid = ParquetFile::getInt64(transition_precursor_col, row, 0, false);
    const bool detecting = ParquetFile::getBool(detecting_col, row, true, true);
    if (detecting)
    {
      transition_counts[pid]++;
    }
  }

  // Read runs list
  const String runs_parquet = base_dir + "/runs/runs.parquet";
  if (!File::exists(runs_parquet))
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Missing runs.parquet in '" + oswpq_dir + "'");
  }
  auto runs_table = ParquetFile::readTable(runs_parquet);
  auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");

  const int64_t num_runs = runs_table->num_rows();
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String run_dir = base_dir + "/runs/run_id=" + String(run_id);
    const String features_path = run_dir + "/features.parquet";
    if (!File::exists(features_path)) continue;

    auto features_table = ParquetFile::readTable(features_path);
    auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
    auto feature_run_id_col = ParquetFile::getColumn(features_table, "run_id");
    auto precursor_id_col = ParquetFile::getColumn(features_table, "precursor_id");
    auto exp_rt_col = ParquetFile::getColumn(features_table, "exp_rt");
    auto ms2_area_col = ParquetFile::getOptionalColumn(features_table, "ms2_area_intensity");
    auto ms2_total_area_col = ParquetFile::getOptionalColumn(features_table, "ms2_total_area_intensity");
    auto ms2_apex_col = ParquetFile::getOptionalColumn(features_table, "ms2_apex_intensity");

    const int64_t num_rows = features_table->num_rows();
    rows_.reserve(rows_.size() + static_cast<size_t>(num_rows));
    for (int64_t row = 0; row < num_rows; ++row)
    {
      Row out;
      out.feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
      out.run_id = ParquetFile::getInt64(feature_run_id_col, row, 0, false);
      out.precursor_id = ParquetFile::getInt64(precursor_id_col, row, 0, false);
      out.exp_rt = ParquetFile::getDouble(exp_rt_col, row, 0.0, true);
      out.ms2_area_intensity = ParquetFile::getDouble(ms2_area_col, row, 0.0, true);
      out.ms2_total_area_intensity = ParquetFile::getDouble(ms2_total_area_col, row, 0.0, true);
      out.ms2_apex_intensity = ParquetFile::getDouble(ms2_apex_col, row, 0.0, true);

      auto it = precursor_info.find(out.precursor_id);
      if (it != precursor_info.end())
      {
        out.precursor_charge = it->second.first;
        out.decoy = it->second.second;
      }

      auto tcit = transition_counts.find(out.precursor_id);
      if (tcit != transition_counts.end()) out.transition_count = tcit->second;

      rows_.push_back(std::move(out));
    }
  }
#else
  (void)oswpq_dir;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}

OpenSwathOSWParquetReader::MS2FeaturesResult OpenSwathOSWParquetReader::fetchMS2Features(const String& oswpq_dir, const String& level, const String& main_score) const
{
#ifdef WITH_PARQUET
  MS2FeaturesResult result;
  std::unique_ptr<File::TempDir> temp_dir;
  const String base_dir = ParquetFile::unzipDirectory(oswpq_dir, temp_dir);

  // Read precursors -> charge / decoy
  const String precursors_path = base_dir + "/library/precursors.parquet";
  const String transitions_path = base_dir + "/library/transitions.parquet";
  if (!File::exists(precursors_path) || !File::exists(transitions_path))
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Missing library parquet files in '" + base_dir + "'");
  }

  auto precursors_table = ParquetFile::readTable(precursors_path);
  auto precursor_id_col = ParquetFile::getColumn(precursors_table, "precursor_id");
  auto charge_col = ParquetFile::getColumn(precursors_table, "charge");
  auto decoy_col = ParquetFile::getOptionalColumn(precursors_table, "decoy");

  std::unordered_map<int64_t, std::pair<int, bool>> precursor_info;
  precursor_info.reserve(precursors_table->num_rows());
  for (int64_t row = 0; row < precursors_table->num_rows(); ++row)
  {
    const int64_t pid = ParquetFile::getInt64(precursor_id_col, row, 0, false);
    const int charge = static_cast<int>(ParquetFile::getInt64(charge_col, row, 0, false));
    const bool decoy = ParquetFile::getBool(decoy_col, row, false, true);
    precursor_info.emplace(pid, std::make_pair(charge, decoy));
  }

  // Count detecting transitions per precursor
  auto transitions_table = ParquetFile::readTable(transitions_path);
  auto transition_precursor_col = ParquetFile::getColumn(transitions_table, "precursor_id");
  auto detecting_col = ParquetFile::getColumn(transitions_table, "detecting");
  std::unordered_map<int64_t, int64_t> transition_counts;
  transition_counts.reserve(transitions_table->num_rows());
  for (int64_t row = 0; row < transitions_table->num_rows(); ++row)
  {
    const int64_t pid = ParquetFile::getInt64(transition_precursor_col, row, 0, false);
    const bool detecting = ParquetFile::getBool(detecting_col, row, true, true);
    if (detecting) transition_counts[pid]++;
  }

  // Read runs and per-run features
  const String runs_parquet = base_dir + "/runs/runs.parquet";
  if (!File::exists(runs_parquet))
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Missing runs.parquet in '" + base_dir + "'");
  }
  auto runs_table = ParquetFile::readTable(runs_parquet);
  auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
  const int64_t num_runs = runs_table->num_rows();

  std::unordered_set<std::string> ms2_set;
  std::unordered_set<std::string> ms1_set;

  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String run_dir = base_dir + "/runs/run_id=" + String(run_id);
    const String features_path = run_dir + "/features.parquet";
    if (!File::exists(features_path)) continue;

    auto features_table = ParquetFile::readTable(features_path);

    // core columns
    auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
    auto precursor_id_col = ParquetFile::getColumn(features_table, "precursor_id");
    auto exp_rt_col = ParquetFile::getColumn(features_table, "exp_rt");

    // collect column names and arrays for score extraction
    const auto& schema = features_table->schema();
    std::vector<std::string> col_names;
    for (const auto& f : schema->fields()) col_names.push_back(f->name());

    // precompute indices of ms2 and ms1 score columns
    std::vector<std::string> ms2_cols;
    std::vector<std::string> ms1_cols;
    for (const auto& name : col_names)
    {
      const std::string lname = name;
      if (lname.rfind("var_ms2_", 0) == 0 || lname.rfind("ms2_", 0) == 0)
      {
        ms2_cols.push_back(lname);
        ms2_set.insert(lname);
      }
      else if (lname.rfind("var_ms1_", 0) == 0 || lname.rfind("ms1_", 0) == 0)
      {
        ms1_cols.push_back(lname);
        ms1_set.insert(lname);
      }
    }

    const int64_t num_rows = features_table->num_rows();
    result.rows.reserve(result.rows.size() + static_cast<size_t>(num_rows));
    for (int64_t row = 0; row < num_rows; ++row)
    {
      OpenSwathOSWParquetReaderRowMS2 out;
      out.feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
      out.run_id = run_id;
      out.precursor_id = ParquetFile::getInt64(precursor_id_col, row, 0, false);
      out.exp_rt = ParquetFile::getDouble(exp_rt_col, row, 0.0, true);
      out.group_id = String(out.run_id) + "_" + String(out.precursor_id);

      auto pit = precursor_info.find(out.precursor_id);
      if (pit != precursor_info.end())
      {
        out.precursor_charge = pit->second.first;
        out.decoy = pit->second.second;
      }
      auto tcit = transition_counts.find(out.precursor_id);
      if (tcit != transition_counts.end()) out.transition_count = tcit->second;

      // extract ms2 scores
      for (const auto& cname : ms2_cols)
      {
        double v = ParquetFile::getDouble(ParquetFile::getColumn(features_table, cname), row, std::numeric_limits<double>::quiet_NaN(), true);
        out.ms2_scores[String(cname)] = v;
      }

      if (level == "ms1ms2")
      {
        for (const auto& cname : ms1_cols)
        {
          double v = ParquetFile::getDouble(ParquetFile::getColumn(features_table, cname), row, std::numeric_limits<double>::quiet_NaN(), true);
          out.ms1_scores[String(cname)] = v;
        }
      }

      result.rows.push_back(std::move(out));
    }
  }

  // sort by run_id, precursor_id, exp_rt
  std::sort(result.rows.begin(), result.rows.end(), [](const OpenSwathOSWParquetReaderRowMS2& a, const OpenSwathOSWParquetReaderRowMS2& b) {
    if (a.run_id != b.run_id) return a.run_id < b.run_id;
    if (a.precursor_id != b.precursor_id) return a.precursor_id < b.precursor_id;
    return a.exp_rt < b.exp_rt;
  });

  // populate column lists
  result.ms2_columns.reserve(ms2_set.size());
  for (const auto &s : ms2_set) result.ms2_columns.emplace_back(s);
  result.ms1_columns.reserve(ms1_set.size());
  for (const auto &s : ms1_set) result.ms1_columns.emplace_back(s);

  // main_score handling is left to caller; we include the maps for flexible selection
  return result;
#else
  (void)level; (void)main_score; (void)oswpq_dir;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}


