// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetReader.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/FORMAT/ZipRandomAccessFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <arrow/api.h>

#include <fstream>
#include <unordered_map>

using namespace OpenMS;

OpenSwathOSWParquetReader::OpenSwathOSWParquetReader(const String& oswpq_dir)
{
  load(oswpq_dir);
}

void OpenSwathOSWParquetReader::load(const String& oswpq_dir)
{
  rows_.clear();
  // remember the provided path for later fetch calls
  oswpq_dir_ = oswpq_dir;
  std::unique_ptr<File::TempDir> temp_dir;
  // Extract small library and runs index files only; per-run parquet files will be
  // read on-demand from the archive when possible to avoid unpacking the whole archive.

  auto open_table_from_entry = [&](const String& entry) -> std::shared_ptr<arrow::Table>
  {
    // Try to open an Arrow RandomAccessFile backed by the archive entry.
    auto ra_res = ZipRandomAccessFile::Open(oswpq_dir, entry, temp_dir);
    if (ra_res.ok())
    {
      const auto& raf = ra_res.ValueOrDie();
      return ParquetFile::readTable(std::static_pointer_cast<arrow::io::RandomAccessFile>(raf));
    }
    // Fallback: extract to temp file and read
    const String path = ZipArchiveFile::extractEntryToTempFile(oswpq_dir, entry, temp_dir);
    if (!File::exists(path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing required parquet file '" + entry + "' in '" + oswpq_dir + "'");
    }
    return ParquetFile::readTable(path);
  };

  // Read precursors -> charge / decoy
  auto precursors_table = open_table_from_entry("library/precursors.parquet");
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
  auto transitions_table = open_table_from_entry("library/transitions.parquet");
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
  auto runs_table = open_table_from_entry("runs/runs.parquet");
  auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");

  const int64_t num_runs = runs_table->num_rows();
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
  const String features_entry = String("runs/run_id=") + String(run_id) + "/features.parquet";
  auto features_table = open_table_from_entry(features_entry);
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
      if (it == precursor_info.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Unknown precursor_id=" + String(out.precursor_id) +
                                            " for run_id=" + String(out.run_id));
      }
      out.precursor_charge = it->second.first;
      out.decoy = it->second.second;

      auto tcit = transition_counts.find(out.precursor_id);
      if (tcit != transition_counts.end()) out.transition_count = tcit->second;

      // provide a stable group id for downstream grouping/join workflows
      // include feature_id to avoid collisions for multiple features with same run/precursor
      out.group_id = String(out.run_id) + "_" + String(out.precursor_id) + "_" + String(out.feature_id);

      rows_.push_back(std::move(out));
    }
  }
}

OpenSwathOSWParquetReader::PeakGroupFeatureScoresResult OpenSwathOSWParquetReader::fetchPeakGroupFeatures(const String& oswpq_dir, const String& level, const String& main_score) const
{
  PeakGroupFeatureScoresResult result;
  std::unique_ptr<File::TempDir> temp_dir;
  // Use RandomAccessFile-backed reads when possible; fall back to extracting
  // entries to temp files when not available.
  auto open_table_from_entry = [&](const String& entry) -> std::shared_ptr<arrow::Table>
  {
    auto ra_res = ZipRandomAccessFile::Open(oswpq_dir, entry, temp_dir);
    if (ra_res.ok())
    {
      const auto& raf = ra_res.ValueOrDie();
      return ParquetFile::readTable(std::static_pointer_cast<arrow::io::RandomAccessFile>(raf));
    }
    const String path = ZipArchiveFile::extractEntryToTempFile(oswpq_dir, entry, temp_dir);
    if (!File::exists(path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing library parquet file '" + entry + "' in '" + oswpq_dir + "'");
    }
    return ParquetFile::readTable(path);
  };

  auto precursors_table = open_table_from_entry("library/precursors.parquet");
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
  auto transitions_table = open_table_from_entry("library/transitions.parquet");
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
  auto runs_table = open_table_from_entry("runs/runs.parquet");
  auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
  const int64_t num_runs = runs_table->num_rows();

  // First pass: discover all score columns across runs (ms2 + optionally ms1)
  std::vector<std::string> all_ms2_cols;
  std::vector<std::string> all_ms1_cols;
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
  const String features_entry = String("runs/run_id=") + String(run_id) + "/features.parquet";
  auto features_table = open_table_from_entry(features_entry);
    const auto& schema = features_table->schema();
    for (const auto& f : schema->fields())
    {
      const std::string name = f->name();
      if (name.rfind("var_ms2_", 0) == 0 || name.rfind("ms2_", 0) == 0)
      {
        if (std::find(all_ms2_cols.begin(), all_ms2_cols.end(), name) == all_ms2_cols.end())
          all_ms2_cols.push_back(name);
      }
      else if (level == "ms1ms2" && (name.rfind("var_ms1_", 0) == 0 || name.rfind("ms1_", 0) == 0))
      {
        if (std::find(all_ms1_cols.begin(), all_ms1_cols.end(), name) == all_ms1_cols.end())
          all_ms1_cols.push_back(name);
      }
    }
  }

  // Sort columns deterministically
  std::sort(all_ms2_cols.begin(), all_ms2_cols.end());
  std::sort(all_ms1_cols.begin(), all_ms1_cols.end());

  // Prepare output vectors
  std::vector<int64_t> feature_id_v;
  std::vector<int64_t> run_id_v;
  std::vector<int64_t> precursor_id_v;
  std::vector<double> exp_rt_v;
  std::vector<int> precursor_charge_v;
  std::vector<bool> decoy_v;
  std::vector<int64_t> transition_count_v;
  std::vector<String> group_id_v;

  std::vector<std::vector<double>> ms2_values(all_ms2_cols.size());
  std::vector<std::vector<double>> ms1_values(all_ms1_cols.size());

  // Second pass: populate vectors
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String features_entry = String("runs/run_id=") + String(run_id) + "/features.parquet";
    auto features_table = open_table_from_entry(features_entry);

    auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
    auto precursor_id_col = ParquetFile::getColumn(features_table, "precursor_id");
    auto exp_rt_col = ParquetFile::getColumn(features_table, "exp_rt");

    // Pre-fetch optional columns for ms2/ms1 (Arrow arrays)
    std::vector<std::shared_ptr<arrow::Array>> ms2_cols_readers;
    ms2_cols_readers.reserve(all_ms2_cols.size());
    for (const auto &cname : all_ms2_cols)
    {
      ms2_cols_readers.push_back(ParquetFile::getOptionalColumn(features_table, cname));
    }
    std::vector<std::shared_ptr<arrow::Array>> ms1_cols_readers;
    ms1_cols_readers.reserve(all_ms1_cols.size());
    for (const auto &cname : all_ms1_cols)
    {
      ms1_cols_readers.push_back(ParquetFile::getOptionalColumn(features_table, cname));
    }

    const int64_t num_rows = features_table->num_rows();
    feature_id_v.reserve(feature_id_v.size() + static_cast<size_t>(num_rows));
    run_id_v.reserve(run_id_v.size() + static_cast<size_t>(num_rows));
    precursor_id_v.reserve(precursor_id_v.size() + static_cast<size_t>(num_rows));
    exp_rt_v.reserve(exp_rt_v.size() + static_cast<size_t>(num_rows));
    precursor_charge_v.reserve(precursor_charge_v.size() + static_cast<size_t>(num_rows));
    decoy_v.reserve(decoy_v.size() + static_cast<size_t>(num_rows));
    transition_count_v.reserve(transition_count_v.size() + static_cast<size_t>(num_rows));
    group_id_v.reserve(group_id_v.size() + static_cast<size_t>(num_rows));

    for (int64_t row = 0; row < num_rows; ++row)
    {
      const int64_t fid = ParquetFile::getInt64(feature_id_col, row, 0, false);
      const int64_t pid = ParquetFile::getInt64(precursor_id_col, row, 0, false);
      const double rt = ParquetFile::getDouble(exp_rt_col, row, 0.0, true);

      feature_id_v.push_back(fid);
      run_id_v.push_back(run_id);
      precursor_id_v.push_back(pid);
      exp_rt_v.push_back(rt);

      auto pit = precursor_info.find(pid);
      if (pit == precursor_info.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Unknown precursor_id=" + String(pid) +
                                            " for run_id=" + String(run_id));
      }
      precursor_charge_v.push_back(pit->second.first);
      decoy_v.push_back(pit->second.second);

      auto tcit = transition_counts.find(pid);
      if (tcit != transition_counts.end()) transition_count_v.push_back(tcit->second);
      else transition_count_v.push_back(0);

  // include feature id in group id to avoid collisions when multiple features exist
  const String gid = String(run_id) + "_" + String(pid) + "_" + String(fid);
      group_id_v.push_back(gid);

      // ms2 columns
      for (size_t ci = 0; ci < all_ms2_cols.size(); ++ci)
      {
        double v = ParquetFile::getDouble(ms2_cols_readers[ci], row, std::numeric_limits<double>::quiet_NaN(), true);
        ms2_values[ci].push_back(v);
      }

      // ms1 columns
      for (size_t ci = 0; ci < all_ms1_cols.size(); ++ci)
      {
        double v = ParquetFile::getDouble(ms1_cols_readers[ci], row, std::numeric_limits<double>::quiet_NaN(), true);
        ms1_values[ci].push_back(v);
      }
    }
  }

  // sort by run_id, precursor_id, exp_rt: produce an index permutation and apply to all vectors
  const size_t N = feature_id_v.size();
  std::vector<size_t> idx(N);
  for (size_t i = 0; i < N; ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){
    if (run_id_v[a] != run_id_v[b]) return run_id_v[a] < run_id_v[b];
    if (precursor_id_v[a] != precursor_id_v[b]) return precursor_id_v[a] < precursor_id_v[b];
    return exp_rt_v[a] < exp_rt_v[b];
  });

  auto permute = [&](auto &vec)
  {
    using T = typename std::decay<decltype(vec)>::type::value_type;
    std::vector<T> tmp;
    tmp.reserve(vec.size());
    for (size_t i = 0; i < idx.size(); ++i) tmp.push_back(vec[idx[i]]);
    vec.swap(tmp);
  };

  permute(feature_id_v);
  permute(run_id_v);
  permute(precursor_id_v);
  permute(exp_rt_v);
  permute(precursor_charge_v);
  permute(decoy_v);
  permute(transition_count_v);
  permute(group_id_v);
  for (auto &col : ms2_values) permute(col);
  for (auto &col : ms1_values) permute(col);

  // fill result
  result.feature_id = std::move(feature_id_v);
  result.run_id = std::move(run_id_v);
  result.precursor_id = std::move(precursor_id_v);
  result.exp_rt = std::move(exp_rt_v);
  result.precursor_charge = std::move(precursor_charge_v);
  result.decoy = std::move(decoy_v);
  result.transition_count = std::move(transition_count_v);
  result.group_id = std::move(group_id_v);

  result.ms2_columns.reserve(all_ms2_cols.size());
  for (const auto &s : all_ms2_cols) result.ms2_columns.emplace_back(s);
  result.ms2_values = std::move(ms2_values);

  result.ms1_columns.reserve(all_ms1_cols.size());
  for (const auto &s : all_ms1_cols) result.ms1_columns.emplace_back(s);
  result.ms1_values = std::move(ms1_values);

  // If a main_score was requested, place it first among ms2 columns/values
  if (!main_score.empty() && !result.ms2_columns.empty())
  {
    SignedSize found = -1;
    for (size_t i = 0; i < result.ms2_columns.size(); ++i)
    {
      if (result.ms2_columns[i] == main_score)
      {
        found = static_cast<SignedSize>(i);
        break;
      }
    }
    if (found > 0)
    {
      // move chosen column to front
      auto col_name = result.ms2_columns[found];
      auto col_values = std::move(result.ms2_values[found]);
      for (size_t i = static_cast<size_t>(found); i > 0; --i)
      {
        result.ms2_columns[i] = result.ms2_columns[i-1];
        result.ms2_values[i] = std::move(result.ms2_values[i-1]);
      }
      result.ms2_columns[0] = col_name;
      result.ms2_values[0] = std::move(col_values);
    }
    else if (found == -1)
    {
      OPENMS_LOG_WARN << "Requested main_score '" << main_score << "' not found among discovered MS2 columns\n";
    }
  }

  return result;
}

OpenSwathOSWParquetReader::TransitionFeaturesResult OpenSwathOSWParquetReader::fetchTransitionFeatures(const String& oswpq_dir) const
{
  TransitionFeaturesResult result;
  std::unique_ptr<File::TempDir> temp_dir;
  // Use RandomAccessFile-backed reads when possible; fall back to extraction when not.
  auto open_table_from_entry = [&](const String& entry) -> std::shared_ptr<arrow::Table>
  {
    auto ra_res = ZipRandomAccessFile::Open(oswpq_dir, entry, temp_dir);
    if (ra_res.ok())
    {
      const auto& raf = ra_res.ValueOrDie();
      return ParquetFile::readTable(std::static_pointer_cast<arrow::io::RandomAccessFile>(raf));
    }
    const String path = ZipArchiveFile::extractEntryToTempFile(oswpq_dir, entry, temp_dir);
    if (!File::exists(path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing parquet file '" + entry + "' in '" + oswpq_dir + "'");
    }
    return ParquetFile::readTable(path);
  };

  auto precursors_table = open_table_from_entry("library/precursors.parquet");
  auto precursor_id_col = ParquetFile::getColumn(precursors_table, "precursor_id");
  auto charge_col = ParquetFile::getColumn(precursors_table, "charge");

  std::unordered_map<int64_t, int> precursor_charge;
  precursor_charge.reserve(precursors_table->num_rows());
  for (int64_t row = 0; row < precursors_table->num_rows(); ++row)
  {
    const int64_t pid = ParquetFile::getInt64(precursor_id_col, row, 0, false);
    const int charge = static_cast<int>(ParquetFile::getInt64(charge_col, row, 0, false));
    precursor_charge.emplace(pid, charge);
  }

  // Read library transitions for product charge and decoy
  auto transitions_table = open_table_from_entry("library/transitions.parquet");
  auto transition_id_col = ParquetFile::getColumn(transitions_table, "transition_id");
  auto transition_charge_col = ParquetFile::getColumn(transitions_table, "charge");
  auto transition_decoy_col = ParquetFile::getOptionalColumn(transitions_table, "decoy");

  std::unordered_map<int64_t, std::pair<int, bool>> transition_info;
  transition_info.reserve(transitions_table->num_rows());
  for (int64_t row = 0; row < transitions_table->num_rows(); ++row)
  {
    const int64_t tid = ParquetFile::getInt64(transition_id_col, row, 0, false);
    const int charge = static_cast<int>(ParquetFile::getInt64(transition_charge_col, row, 0, false));
    const bool dec = ParquetFile::getBool(transition_decoy_col, row, false, true);
    transition_info.emplace(tid, std::make_pair(charge, dec));
  }

  // Read runs list
  auto runs_table = open_table_from_entry("runs/runs.parquet");
  auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
  const int64_t num_runs = runs_table->num_rows();

  // First pass: discover transition-level var columns across runs
  std::vector<std::string> all_var_cols;
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String ft_entry = String("runs/run_id=") + String(run_id) + "/feature_transition.parquet";
    auto ft_table = open_table_from_entry(ft_entry);
    const auto& schema = ft_table->schema();
    for (const auto& f : schema->fields())
    {
      const std::string name = f->name();
      if (name.rfind("var_", 0) == 0)
      {
        if (std::find(all_var_cols.begin(), all_var_cols.end(), name) == all_var_cols.end())
          all_var_cols.push_back(name);
      }
    }
  }

  std::sort(all_var_cols.begin(), all_var_cols.end());
  result.transition_var_columns.reserve(all_var_cols.size());
  for (const auto &s : all_var_cols) result.transition_var_columns.emplace_back(s);
  result.transition_var_values.assign(all_var_cols.size(), std::vector<double>());

  // Second pass: populate vectors
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String ft_entry2 = String("runs/run_id=") + String(run_id) + "/feature_transition.parquet";
    const String features_entry2 = String("runs/run_id=") + String(run_id) + "/features.parquet";

    // Build mapping feature_id -> (precursor_id, exp_rt) from features.parquet
    // Fail fast if the features.parquet file is missing to avoid silently
    // emitting transition rows with defaulted precursor/RT values.
    std::unordered_map<int64_t, std::pair<int64_t, double>> feature_map;
    {
      auto features_table = open_table_from_entry(features_entry2);
      auto feat_id_col = ParquetFile::getColumn(features_table, "feature_id");
      auto precursor_id_col = ParquetFile::getColumn(features_table, "precursor_id");
      auto exp_rt_col = ParquetFile::getColumn(features_table, "exp_rt");
      const int64_t num_feat_rows = features_table->num_rows();
      for (int64_t i = 0; i < num_feat_rows; ++i)
      {
        const int64_t fid = ParquetFile::getInt64(feat_id_col, i, 0, false);
        const int64_t pid = ParquetFile::getInt64(precursor_id_col, i, 0, false);
        const double rt = ParquetFile::getDouble(exp_rt_col, i, 0.0, true);
        feature_map.emplace(fid, std::make_pair(pid, rt));
      }
    }

    auto ft_table = open_table_from_entry(ft_entry2);
    auto ft_feature_id_col = ParquetFile::getColumn(ft_table, "feature_id");
    auto ft_transition_id_col = ParquetFile::getColumn(ft_table, "transition_id");
    auto ft_area_col = ParquetFile::getOptionalColumn(ft_table, "area_intensity");
    auto ft_total_area_col = ParquetFile::getOptionalColumn(ft_table, "total_area_intensity");
    auto ft_apex_col = ParquetFile::getOptionalColumn(ft_table, "apex_intensity");
    auto ft_apex_rt_col = ParquetFile::getOptionalColumn(ft_table, "apex_rt");
    auto ft_rt_fwhm_col = ParquetFile::getOptionalColumn(ft_table, "rt_fwhm");
    auto ft_masserror_col = ParquetFile::getOptionalColumn(ft_table, "masserror_ppm");
    auto ft_total_mi_col = ParquetFile::getOptionalColumn(ft_table, "total_mi");

    // Pre-fetch var columns
    std::vector<std::shared_ptr<arrow::Array>> var_cols_readers;
    var_cols_readers.reserve(all_var_cols.size());
    for (const auto &cname : all_var_cols)
    {
      var_cols_readers.push_back(ParquetFile::getOptionalColumn(ft_table, cname));
    }

    const int64_t num_rows = ft_table->num_rows();
    result.feature_id.reserve(result.feature_id.size() + static_cast<size_t>(num_rows));
    result.run_id.reserve(result.run_id.size() + static_cast<size_t>(num_rows));
    result.precursor_id.reserve(result.precursor_id.size() + static_cast<size_t>(num_rows));
    result.exp_rt.reserve(result.exp_rt.size() + static_cast<size_t>(num_rows));
    result.precursor_charge.reserve(result.precursor_charge.size() + static_cast<size_t>(num_rows));
    result.transition_id.reserve(result.transition_id.size() + static_cast<size_t>(num_rows));
    result.product_charge.reserve(result.product_charge.size() + static_cast<size_t>(num_rows));
    result.decoy.reserve(result.decoy.size() + static_cast<size_t>(num_rows));
    result.area_intensity.reserve(result.area_intensity.size() + static_cast<size_t>(num_rows));
    result.total_area_intensity.reserve(result.total_area_intensity.size() + static_cast<size_t>(num_rows));
    result.apex_intensity.reserve(result.apex_intensity.size() + static_cast<size_t>(num_rows));
    result.apex_rt.reserve(result.apex_rt.size() + static_cast<size_t>(num_rows));
    result.rt_fwhm.reserve(result.rt_fwhm.size() + static_cast<size_t>(num_rows));
    result.masserror_ppm.reserve(result.masserror_ppm.size() + static_cast<size_t>(num_rows));
    result.total_mi.reserve(result.total_mi.size() + static_cast<size_t>(num_rows));
    result.group_id.reserve(result.group_id.size() + static_cast<size_t>(num_rows));
    for (auto &col : result.transition_var_values) col.reserve(col.size() + static_cast<size_t>(num_rows));

    for (int64_t row = 0; row < num_rows; ++row)
    {
      const int64_t fid = ParquetFile::getInt64(ft_feature_id_col, row, 0, false);
      const int64_t tid = ParquetFile::getInt64(ft_transition_id_col, row, 0, false);

      auto fit = feature_map.find(fid);
      if (fit == feature_map.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "feature_transition.parquet references unknown feature_id=" + String(fid) +
                                            " for run_id=" + String(run_id));
      }
      const int64_t pid = fit->second.first;
      const double rt = fit->second.second;

      result.feature_id.push_back(fid);
      result.run_id.push_back(run_id);
      result.precursor_id.push_back(pid);
      result.exp_rt.push_back(rt);

      auto pcit = precursor_charge.find(pid);
      if (pcit == precursor_charge.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Unknown precursor_id=" + String(pid) + " for run_id=" + String(run_id));
      }
      result.precursor_charge.push_back(pcit->second);

      result.transition_id.push_back(tid);
      auto ticit = transition_info.find(tid);
      if (ticit == transition_info.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Unknown transition_id=" + String(tid) + " for run_id=" + String(run_id));
      }
      result.product_charge.push_back(ticit->second.first);
      result.decoy.push_back(ticit->second.second);

      const double area = ParquetFile::getDouble(ft_area_col, row, std::numeric_limits<double>::quiet_NaN(), true);
      const double total_area = ParquetFile::getDouble(ft_total_area_col, row, std::numeric_limits<double>::quiet_NaN(), true);
      const double apex = ParquetFile::getDouble(ft_apex_col, row, std::numeric_limits<double>::quiet_NaN(), true);
      const double apexrt = ParquetFile::getDouble(ft_apex_rt_col, row, std::numeric_limits<double>::quiet_NaN(), true);
      const double fwhm = ParquetFile::getDouble(ft_rt_fwhm_col, row, std::numeric_limits<double>::quiet_NaN(), true);
      const double masserr = ParquetFile::getDouble(ft_masserror_col, row, std::numeric_limits<double>::quiet_NaN(), true);
      const double tmi = ParquetFile::getDouble(ft_total_mi_col, row, std::numeric_limits<double>::quiet_NaN(), true);

      result.area_intensity.push_back(area);
      result.total_area_intensity.push_back(total_area);
      result.apex_intensity.push_back(apex);
      result.apex_rt.push_back(apexrt);
      result.rt_fwhm.push_back(fwhm);
      result.masserror_ppm.push_back(masserr);
      result.total_mi.push_back(tmi);

      const String gid = String(run_id) + "_" + String(fid) + "_" + String(pid) + "_" + String(tid);
      result.group_id.emplace_back(gid);

      // var columns
      for (size_t ci = 0; ci < var_cols_readers.size(); ++ci)
      {
        double v = ParquetFile::getDouble(var_cols_readers[ci], row, std::numeric_limits<double>::quiet_NaN(), true);
        result.transition_var_values[ci].push_back(v);
      }
    }
  }

  // sort by run_id, precursor_id, exp_rt, transition_id
  const size_t N = result.feature_id.size();
  std::vector<size_t> idx(N);
  for (size_t i = 0; i < N; ++i) idx[i] = i;
  std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){
    if (result.run_id[a] != result.run_id[b]) return result.run_id[a] < result.run_id[b];
    if (result.precursor_id[a] != result.precursor_id[b]) return result.precursor_id[a] < result.precursor_id[b];
    if (result.exp_rt[a] != result.exp_rt[b]) return result.exp_rt[a] < result.exp_rt[b];
    return result.transition_id[a] < result.transition_id[b];
  });

  auto permute = [&](auto &vec)
  {
    using T = typename std::decay<decltype(vec)>::type::value_type;
    std::vector<T> tmp;
    tmp.reserve(vec.size());
    for (size_t i = 0; i < idx.size(); ++i) tmp.push_back(vec[idx[i]]);
    vec.swap(tmp);
  };

  permute(result.feature_id);
  permute(result.run_id);
  permute(result.precursor_id);
  permute(result.exp_rt);
  permute(result.precursor_charge);
  permute(result.transition_id);
  permute(result.product_charge);
  permute(result.decoy);
  permute(result.area_intensity);
  permute(result.total_area_intensity);
  permute(result.apex_intensity);
  permute(result.apex_rt);
  permute(result.rt_fwhm);
  permute(result.masserror_ppm);
  permute(result.total_mi);
  permute(result.group_id);
  for (auto &col : result.transition_var_values) permute(col);

  return result;
}

OpenSwathOSWParquetReader::UnscoredResult OpenSwathOSWParquetReader::fetchUnscoredData(const String& oswpq_dir) const
{
  UnscoredResult result;
  std::unique_ptr<File::TempDir> temp_dir;
  // Use RandomAccessFile-backed reads when possible; fall back to extraction when not.
  auto open_table_from_entry = [&](const String& entry) -> std::shared_ptr<arrow::Table>
  {
    auto ra_res = ZipRandomAccessFile::Open(oswpq_dir, entry, temp_dir);
    if (ra_res.ok())
    {
      const auto& raf = ra_res.ValueOrDie();
      return ParquetFile::readTable(std::static_pointer_cast<arrow::io::RandomAccessFile>(raf));
    }
    const String path = ZipArchiveFile::extractEntryToTempFile(oswpq_dir, entry, temp_dir);
    if (!File::exists(path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing parquet file '" + entry + "' in '" + oswpq_dir + "'");
    }
    return ParquetFile::readTable(path);
  };

  auto precursors_table = open_table_from_entry("library/precursors.parquet");
  auto precursor_id_col = ParquetFile::getColumn(precursors_table, "precursor_id");
  auto charge_col = ParquetFile::getOptionalColumn(precursors_table, "charge");
  auto decoy_col = ParquetFile::getOptionalColumn(precursors_table, "decoy");
  auto precursor_mz_col = ParquetFile::getOptionalColumn(precursors_table, "precursor_mz");
  auto library_rt_col = ParquetFile::getOptionalColumn(precursors_table, "library_rt");

  // build precursor info maps
  std::unordered_map<int64_t, int> precursor_charge;
  std::unordered_map<int64_t, bool> precursor_decoy;
  std::unordered_map<int64_t, double> precursor_mz;
  std::unordered_map<int64_t, double> precursor_library_rt;

  for (int64_t row = 0; row < precursors_table->num_rows(); ++row)
  {
    const int64_t pid = ParquetFile::getInt64(precursor_id_col, row, 0, false);
    if (charge_col) precursor_charge[pid] = static_cast<int>(ParquetFile::getInt64(charge_col, row, 0, false));
    if (decoy_col) precursor_decoy[pid] = ParquetFile::getBool(decoy_col, row, false, true);
    if (precursor_mz_col) precursor_mz[pid] = ParquetFile::getDouble(precursor_mz_col, row, std::numeric_limits<double>::quiet_NaN(), true);
    if (library_rt_col) precursor_library_rt[pid] = ParquetFile::getDouble(library_rt_col, row, std::numeric_limits<double>::quiet_NaN(), true);
  }

  // try to read peptide mapping if present: precursor_peptide_mapping.parquet and peptides.parquet
  std::unordered_map<int64_t, int64_t> precursor_to_peptide;
  bool has_ppm = false;
  bool has_pep = false;
  if (File::isDirectory(oswpq_dir))
  {
    // Directory input: check files directly
    has_ppm = File::exists(oswpq_dir + "/library/precursor_peptide_mapping.parquet");
    has_pep = File::exists(oswpq_dir + "/library/peptides.parquet");
  }
  else
  {
    auto entries = ZipArchiveFile::listEntries(oswpq_dir);
    has_ppm = std::find(entries.begin(), entries.end(), std::string("library/precursor_peptide_mapping.parquet")) != entries.end();
    has_pep = std::find(entries.begin(), entries.end(), std::string("library/peptides.parquet")) != entries.end();
  }
  if (has_ppm && has_pep)
  {
    auto ppm_table = open_table_from_entry("library/precursor_peptide_mapping.parquet");
    auto pep_table = open_table_from_entry("library/peptides.parquet");
    (void)pep_table; // loaded to validate presence; fields used via ppm_table
    auto ppm_prec_col = ParquetFile::getColumn(ppm_table, "precursor_id");
    auto ppm_pep_col = ParquetFile::getColumn(ppm_table, "peptide_id");
    for (int64_t r = 0; r < ppm_table->num_rows(); ++r)
    {
      const int64_t pid = ParquetFile::getInt64(ppm_prec_col, r, 0, false);
      const int64_t pep = ParquetFile::getInt64(ppm_pep_col, r, 0, false);
      precursor_to_peptide[pid] = pep;
    }
  }

  // Read runs
  auto runs_table = open_table_from_entry("runs/runs.parquet");
  auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
  auto run_filename_col = ParquetFile::getOptionalColumn(runs_table, "filename");
  const int64_t num_runs = runs_table->num_rows();

  // First pass: discover ms2 and ms1 score columns and IM presence
  std::vector<std::string> all_ms2_cols;
  std::vector<std::string> all_ms1_cols;

  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String features_entry = String("runs/run_id=") + String(run_id) + "/features.parquet";
    auto features_table = open_table_from_entry(features_entry);
    const auto& schema = features_table->schema();
    for (const auto &f : schema->fields())
    {
      const std::string name = f->name();
      if (name.rfind("var_ms2_", 0) == 0 || name.rfind("ms2_", 0) == 0)
      {
        if (std::find(all_ms2_cols.begin(), all_ms2_cols.end(), name) == all_ms2_cols.end())
          all_ms2_cols.push_back(name);
      }
      else if (name.rfind("var_ms1_", 0) == 0 || name.rfind("ms1_", 0) == 0)
      {
        if (std::find(all_ms1_cols.begin(), all_ms1_cols.end(), name) == all_ms1_cols.end())
          all_ms1_cols.push_back(name);
      }
    }
  }

  std::sort(all_ms2_cols.begin(), all_ms2_cols.end());
  std::sort(all_ms1_cols.begin(), all_ms1_cols.end());

  // Prepare ms score storage
  result.ms2_columns.reserve(all_ms2_cols.size());
  for (const auto &s : all_ms2_cols) result.ms2_columns.emplace_back(s);
  result.ms2_values.assign(all_ms2_cols.size(), std::vector<double>());
  result.ms1_columns.reserve(all_ms1_cols.size());
  for (const auto &s : all_ms1_cols) result.ms1_columns.emplace_back(s);
  result.ms1_values.assign(all_ms1_cols.size(), std::vector<double>());

  // Second pass: populate rows
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String features_entry = String("runs/run_id=") + String(run_id) + "/features.parquet";
    auto features_table = open_table_from_entry(features_entry);

    auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
    auto feature_prec_col = ParquetFile::getColumn(features_table, "precursor_id");
    auto exp_rt_col = ParquetFile::getOptionalColumn(features_table, "exp_rt");
    auto delta_rt_col = ParquetFile::getOptionalColumn(features_table, "delta_rt");
    auto norm_rt_col = ParquetFile::getOptionalColumn(features_table, "norm_rt");
    auto ms2_area_col = ParquetFile::getOptionalColumn(features_table, "ms2_area_intensity");
    auto ms1_area_col = ParquetFile::getOptionalColumn(features_table, "ms1_area_intensity");
    auto ms1_apex_col = ParquetFile::getOptionalColumn(features_table, "ms1_apex_intensity");
    auto left_col = ParquetFile::getOptionalColumn(features_table, "left_width");
    auto right_col = ParquetFile::getOptionalColumn(features_table, "right_width");
    auto exp_im_col = ParquetFile::getOptionalColumn(features_table, "exp_im");
    auto exp_im_l_col = ParquetFile::getOptionalColumn(features_table, "exp_im_leftwidth");
    auto exp_im_r_col = ParquetFile::getOptionalColumn(features_table, "exp_im_rightwidth");

    // prefetch ms1/ms2 optional columns
    std::vector<std::shared_ptr<arrow::Array>> ms2_cols_readers;
    for (const auto &c : all_ms2_cols) ms2_cols_readers.push_back(ParquetFile::getOptionalColumn(features_table, c));
    std::vector<std::shared_ptr<arrow::Array>> ms1_cols_readers;
    for (const auto &c : all_ms1_cols) ms1_cols_readers.push_back(ParquetFile::getOptionalColumn(features_table, c));

    const int64_t num_rows = features_table->num_rows();
    for (int64_t i = 0; i < num_rows; ++i)
    {
      const int64_t fid = ParquetFile::getInt64(feature_id_col, i, 0, false);
      const int64_t pid = ParquetFile::getInt64(feature_prec_col, i, 0, false);
      const double rt = ParquetFile::getDouble(exp_rt_col, i, std::numeric_limits<double>::quiet_NaN(), true);
      const double drt = ParquetFile::getDouble(delta_rt_col, i, std::numeric_limits<double>::quiet_NaN(), true);
      const double nrt = ParquetFile::getDouble(norm_rt_col, i, std::numeric_limits<double>::quiet_NaN(), true);

      result.id_run.push_back(run_id);
      result.id_peptide.push_back(precursor_to_peptide.count(pid) ? precursor_to_peptide[pid] : 0);
      result.transition_group_id.push_back(pid);
      result.decoy.push_back(precursor_decoy.count(pid) ? precursor_decoy[pid] : false);
      result.run_id.push_back(run_id);
  result.filename.push_back(run_filename_col ? String(ParquetFile::getString(run_filename_col, r)) : String(""));
      result.RT.push_back(rt);
      result.assay_rt.push_back(rt - drt);
      result.delta_rt.push_back(drt);
      result.assay_RT.push_back(precursor_library_rt.count(pid) ? precursor_library_rt[pid] : std::numeric_limits<double>::quiet_NaN());
      result.delta_RT.push_back(nrt - (precursor_library_rt.count(pid) ? precursor_library_rt[pid] : std::numeric_limits<double>::quiet_NaN()));
      result.id.push_back(fid);
      result.Charge.push_back(precursor_charge.count(pid) ? precursor_charge[pid] : 0);
      result.mz.push_back(precursor_mz.count(pid) ? precursor_mz[pid] : std::numeric_limits<double>::quiet_NaN());
      result.Intensity.push_back(ParquetFile::getDouble(ms2_area_col, i, std::numeric_limits<double>::quiet_NaN(), true));
      result.aggr_prec_Peak_Area.push_back(ParquetFile::getDouble(ms1_area_col, i, std::numeric_limits<double>::quiet_NaN(), true));
      result.aggr_prec_Peak_Apex.push_back(ParquetFile::getDouble(ms1_apex_col, i, std::numeric_limits<double>::quiet_NaN(), true));
      result.leftWidth.push_back(ParquetFile::getDouble(left_col, i, std::numeric_limits<double>::quiet_NaN(), true));
      result.rightWidth.push_back(ParquetFile::getDouble(right_col, i, std::numeric_limits<double>::quiet_NaN(), true));

      result.EXP_IM.push_back(ParquetFile::getDouble(exp_im_col, i, std::numeric_limits<double>::quiet_NaN(), true));
      result.IM_leftWidth.push_back(ParquetFile::getDouble(exp_im_l_col, i, std::numeric_limits<double>::quiet_NaN(), true));
      result.IM_rightWidth.push_back(ParquetFile::getDouble(exp_im_r_col, i, std::numeric_limits<double>::quiet_NaN(), true));

      // ms2 scores
      for (size_t ci = 0; ci < ms2_cols_readers.size(); ++ci)
      {
        double v = ParquetFile::getDouble(ms2_cols_readers[ci], i, std::numeric_limits<double>::quiet_NaN(), true);
        result.ms2_values[ci].push_back(v);
      }
      for (size_t ci = 0; ci < ms1_cols_readers.size(); ++ci)
      {
        double v = ParquetFile::getDouble(ms1_cols_readers[ci], i, std::numeric_limits<double>::quiet_NaN(), true);
        result.ms1_values[ci].push_back(v);
      }
    }
  }

  return result;
}



