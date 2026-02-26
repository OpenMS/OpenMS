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
    if (!File::exists(features_path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing features.parquet for run_id=" + String(run_id) + " in '" + oswpq_dir + "'");
    }

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

OpenSwathOSWParquetReader::PeakGroupFeatureScoresResult OpenSwathOSWParquetReader::fetchPeakGroupFeatures(const String& oswpq_dir, const String& level, const String& main_score) const
{
#ifdef WITH_PARQUET
  PeakGroupFeatureScoresResult result;
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

  // First pass: discover all score columns across runs (ms2 + optionally ms1)
  std::vector<std::string> all_ms2_cols;
  std::vector<std::string> all_ms1_cols;
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String run_dir = base_dir + "/runs/run_id=" + String(run_id);
    const String features_path = run_dir + "/features.parquet";
    if (!File::exists(features_path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing features.parquet for run_id=" + String(run_id) + " in '" + base_dir + "'");
    }
    auto features_table = ParquetFile::readTable(features_path);
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
    const String run_dir = base_dir + "/runs/run_id=" + String(run_id);
    const String features_path = run_dir + "/features.parquet";
    if (!File::exists(features_path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing features.parquet for run_id=" + String(run_id) + " in '" + base_dir + "'");
    }
    auto features_table = ParquetFile::readTable(features_path);

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
      if (pit != precursor_info.end())
      {
        precursor_charge_v.push_back(pit->second.first);
        decoy_v.push_back(pit->second.second);
      }
      else
      {
        precursor_charge_v.push_back(0);
        decoy_v.push_back(false);
      }

      auto tcit = transition_counts.find(pid);
      if (tcit != transition_counts.end()) transition_count_v.push_back(tcit->second);
      else transition_count_v.push_back(0);

      const String gid = String(run_id) + "_" + String(pid);
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
  for (const auto &s : all_ms2_cols) result.ms2_columns.emplace_back(String(s));
  result.ms2_values = std::move(ms2_values);

  result.ms1_columns.reserve(all_ms1_cols.size());
  for (const auto &s : all_ms1_cols) result.ms1_columns.emplace_back(String(s));
  result.ms1_values = std::move(ms1_values);

  return result;
#else
  (void)level; (void)main_score; (void)oswpq_dir;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}

OpenSwathOSWParquetReader::TransitionFeaturesResult OpenSwathOSWParquetReader::fetchTransitionFeatures(const String& oswpq_dir) const
{
#ifdef WITH_PARQUET
  TransitionFeaturesResult result;
  std::unique_ptr<File::TempDir> temp_dir;
  const String base_dir = ParquetFile::unzipDirectory(oswpq_dir, temp_dir);

  // Read precursors -> charge
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

  std::unordered_map<int64_t, int> precursor_charge;
  precursor_charge.reserve(precursors_table->num_rows());
  for (int64_t row = 0; row < precursors_table->num_rows(); ++row)
  {
    const int64_t pid = ParquetFile::getInt64(precursor_id_col, row, 0, false);
    const int charge = static_cast<int>(ParquetFile::getInt64(charge_col, row, 0, false));
    precursor_charge.emplace(pid, charge);
  }

  // Read library transitions for product charge and decoy
  auto transitions_table = ParquetFile::readTable(transitions_path);
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
  const String runs_parquet = base_dir + "/runs/runs.parquet";
  if (!File::exists(runs_parquet))
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Missing runs.parquet in '" + base_dir + "'");
  }
  auto runs_table = ParquetFile::readTable(runs_parquet);
  auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
  const int64_t num_runs = runs_table->num_rows();

  // First pass: discover transition-level var columns across runs
  std::vector<std::string> all_var_cols;
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String run_dir = base_dir + "/runs/run_id=" + String(run_id);
    const String ft_path = run_dir + "/feature_transition.parquet";
    if (!File::exists(ft_path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing feature_transition.parquet for run_id=" + String(run_id) + " in '" + base_dir + "'");
    }
    auto ft_table = ParquetFile::readTable(ft_path);
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
  for (const auto &s : all_var_cols) result.transition_var_columns.emplace_back(String(s));
  result.transition_var_values.assign(all_var_cols.size(), std::vector<double>());

  // Second pass: populate vectors
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String run_dir = base_dir + "/runs/run_id=" + String(run_id);
    const String features_path = run_dir + "/features.parquet";
    const String ft_path = run_dir + "/feature_transition.parquet";
    if (!File::exists(ft_path))
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Missing feature_transition.parquet for run_id=" + String(run_id) + " in '" + base_dir + "'");
    }

    // Build mapping feature_id -> (precursor_id, exp_rt) from features.parquet
    std::unordered_map<int64_t, std::pair<int64_t, double>> feature_map;
    if (File::exists(features_path))
    {
      auto features_table = ParquetFile::readTable(features_path);
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

    auto ft_table = ParquetFile::readTable(ft_path);
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

      int64_t pid = 0;
      double rt = 0.0;
      auto fit = feature_map.find(fid);
      if (fit != feature_map.end())
      {
        pid = fit->second.first;
        rt = fit->second.second;
      }

      result.feature_id.push_back(fid);
      result.run_id.push_back(run_id);
      result.precursor_id.push_back(pid);
      result.exp_rt.push_back(rt);

      auto pcit = precursor_charge.find(pid);
      if (pcit != precursor_charge.end()) result.precursor_charge.push_back(pcit->second);
      else result.precursor_charge.push_back(0);

      result.transition_id.push_back(tid);
      auto ticit = transition_info.find(tid);
      if (ticit != transition_info.end())
      {
        result.product_charge.push_back(ticit->second.first);
        result.decoy.push_back(ticit->second.second);
      }
      else
      {
        result.product_charge.push_back(0);
        result.decoy.push_back(false);
      }

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
#else
  (void)oswpq_dir;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}



