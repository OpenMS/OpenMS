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

void OpenSwathOSWParquetReader::writeCSV(const String& filename) const
{
  std::ofstream out(filename.c_str(), std::ios::out | std::ios::trunc);
  if (!out.is_open())
  {
    throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
  }
  out << "feature_id,run_id,precursor_id,exp_rt,ms2_area,ms2_total_area,ms2_apex,precursor_charge,decoy,transition_count\n";
  for (const auto& r : rows_)
  {
    out << r.feature_id << "," << r.run_id << "," << r.precursor_id << "," << r.exp_rt << ","
        << r.ms2_area_intensity << "," << r.ms2_total_area_intensity << "," << r.ms2_apex_intensity << ","
        << r.precursor_charge << "," << (r.decoy ? 1 : 0) << "," << r.transition_count << "\n";
  }
}

void OpenSwathOSWParquetReader::fetchMS2Features(const String& oswpq_dir, std::vector<OpenSwathOSWParquetReaderRowMS2>& out_rows, const String& level, const String& main_score) const
{
#ifdef WITH_PARQUET
  out_rows.clear();
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
      }
      else if (lname.rfind("var_ms1_", 0) == 0 || lname.rfind("ms1_", 0) == 0)
      {
        ms1_cols.push_back(lname);
      }
    }

    const int64_t num_rows = features_table->num_rows();
    out_rows.reserve(out_rows.size() + static_cast<size_t>(num_rows));
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

      out_rows.push_back(std::move(out));
    }
  }

  // sort by run_id, precursor_id, exp_rt
  std::sort(out_rows.begin(), out_rows.end(), [](const OpenSwathOSWParquetReaderRowMS2& a, const OpenSwathOSWParquetReaderRowMS2& b) {
    if (a.run_id != b.run_id) return a.run_id < b.run_id;
    if (a.precursor_id != b.precursor_id) return a.precursor_id < b.precursor_id;
    return a.exp_rt < b.exp_rt;
  });

  // main_score handling is left to caller; we include the maps for flexible selection
#else
  (void)out_rows; (void)level; (void)main_score;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}

void OpenSwathOSWParquetReader::writeMS2FeaturesParquet(const String& oswpq_dir, const String& output_path, const String& level) const
{
#ifdef WITH_PARQUET
  std::unique_ptr<File::TempDir> temp_dir;
  const String base_dir = ParquetFile::unzipDirectory(oswpq_dir, temp_dir);

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

  const String runs_parquet = base_dir + "/runs/runs.parquet";
  if (!File::exists(runs_parquet))
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Missing runs.parquet in '" + base_dir + "'");
  }
  auto runs_table = ParquetFile::readTable(runs_parquet);
  auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");

  // First pass: discover all score columns across runs (ms2 + optionally ms1)
  std::vector<std::string> all_ms2_cols;
  std::vector<std::string> all_ms1_cols;
  const int64_t num_runs = runs_table->num_rows();
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String run_dir = base_dir + "/runs/run_id=" + String(run_id);
    const String features_path = run_dir + "/features.parquet";
    if (!File::exists(features_path)) continue;
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

  // Prepare builders
  arrow::Int64Builder feature_id_b;
  arrow::Int64Builder run_id_b;
  arrow::Int64Builder precursor_id_b;
  arrow::DoubleBuilder exp_rt_b;
  arrow::Int64Builder precursor_charge_b;
  arrow::BooleanBuilder decoy_b;
  arrow::Int64Builder transition_count_b;
  arrow::StringBuilder group_id_b;

  std::vector<std::unique_ptr<arrow::DoubleBuilder>> ms2_builders;
  for (size_t i = 0; i < all_ms2_cols.size(); ++i)
  {
    ms2_builders.emplace_back(std::make_unique<arrow::DoubleBuilder>());
  }
  std::vector<std::unique_ptr<arrow::DoubleBuilder>> ms1_builders;
  for (size_t i = 0; i < all_ms1_cols.size(); ++i)
  {
    ms1_builders.emplace_back(std::make_unique<arrow::DoubleBuilder>());
  }

  // Second pass: populate builders
  for (int64_t r = 0; r < num_runs; ++r)
  {
    const int64_t run_id = ParquetFile::getInt64(run_id_col, r, 0, false);
    const String run_dir = base_dir + "/runs/run_id=" + String(run_id);
    const String features_path = run_dir + "/features.parquet";
    if (!File::exists(features_path)) continue;
    auto features_table = ParquetFile::readTable(features_path);

    auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
    auto precursor_id_col = ParquetFile::getColumn(features_table, "precursor_id");
    auto exp_rt_col = ParquetFile::getColumn(features_table, "exp_rt");

    const int64_t num_rows = features_table->num_rows();
    for (int64_t row = 0; row < num_rows; ++row)
    {
      const int64_t fid = ParquetFile::getInt64(feature_id_col, row, 0, false);
      const int64_t pid = ParquetFile::getInt64(precursor_id_col, row, 0, false);
      const double rt = ParquetFile::getDouble(exp_rt_col, row, 0.0, true);

      ParquetFile::appendOrThrow(feature_id_b.Append(fid), "feature_id");
      ParquetFile::appendOrThrow(run_id_b.Append(run_id), "run_id");
      ParquetFile::appendOrThrow(precursor_id_b.Append(pid), "precursor_id");
      ParquetFile::appendOrThrow(exp_rt_b.Append(rt), "exp_rt");

      auto pit = precursor_info.find(pid);
      if (pit != precursor_info.end())
      {
        ParquetFile::appendOrThrow(precursor_charge_b.Append(pit->second.first), "precursor_charge");
        ParquetFile::appendOrThrow(decoy_b.Append(pit->second.second), "decoy");
      }
      else
      {
        ParquetFile::appendOrThrow(precursor_charge_b.AppendNull(), "precursor_charge");
        ParquetFile::appendOrThrow(decoy_b.AppendNull(), "decoy");
      }

      auto tcit = transition_counts.find(pid);
      if (tcit != transition_counts.end())
      {
        ParquetFile::appendOrThrow(transition_count_b.Append(tcit->second), "transition_count");
      }
      else
      {
        ParquetFile::appendOrThrow(transition_count_b.AppendNull(), "transition_count");
      }

      const String group_id = String(run_id) + "_" + String(pid);
      ParquetFile::appendOrThrow(group_id_b.Append(std::string(group_id)), "group_id");

      // ms2 columns
      for (size_t ci = 0; ci < all_ms2_cols.size(); ++ci)
      {
        const std::string& cname = all_ms2_cols[ci];
        double v = ParquetFile::getDouble(ParquetFile::getOptionalColumn(features_table, cname), row, std::numeric_limits<double>::quiet_NaN(), true);
        if (std::isnan(v)) ParquetFile::appendOrThrow(ms2_builders[ci]->AppendNull(), cname);
        else ParquetFile::appendOrThrow(ms2_builders[ci]->Append(v), cname);
      }

      // ms1 columns
      for (size_t ci = 0; ci < all_ms1_cols.size(); ++ci)
      {
        const std::string& cname = all_ms1_cols[ci];
        double v = ParquetFile::getDouble(ParquetFile::getOptionalColumn(features_table, cname), row, std::numeric_limits<double>::quiet_NaN(), true);
        if (std::isnan(v)) ParquetFile::appendOrThrow(ms1_builders[ci]->AppendNull(), cname);
        else ParquetFile::appendOrThrow(ms1_builders[ci]->Append(v), cname);
      }
    }
  }

  // Finish arrays
  std::shared_ptr<arrow::Array> feature_id_arr;
  std::shared_ptr<arrow::Array> run_id_arr;
  std::shared_ptr<arrow::Array> precursor_id_arr;
  std::shared_ptr<arrow::Array> exp_rt_arr;
  std::shared_ptr<arrow::Array> precursor_charge_arr;
  std::shared_ptr<arrow::Array> decoy_arr;
  std::shared_ptr<arrow::Array> transition_count_arr;
  std::shared_ptr<arrow::Array> group_id_arr;

  ParquetFile::appendOrThrow(feature_id_b.Finish(&feature_id_arr), "feature_id");
  ParquetFile::appendOrThrow(run_id_b.Finish(&run_id_arr), "run_id");
  ParquetFile::appendOrThrow(precursor_id_b.Finish(&precursor_id_arr), "precursor_id");
  ParquetFile::appendOrThrow(exp_rt_b.Finish(&exp_rt_arr), "exp_rt");
  ParquetFile::appendOrThrow(precursor_charge_b.Finish(&precursor_charge_arr), "precursor_charge");
  ParquetFile::appendOrThrow(decoy_b.Finish(&decoy_arr), "decoy");
  ParquetFile::appendOrThrow(transition_count_b.Finish(&transition_count_arr), "transition_count");
  ParquetFile::appendOrThrow(group_id_b.Finish(&group_id_arr), "group_id");

  std::vector<std::shared_ptr<arrow::Array>> arrays;
  std::vector<std::shared_ptr<arrow::Field>> fields;
  arrays.push_back(feature_id_arr);
  fields.push_back(arrow::field("feature_id", arrow::int64()));
  arrays.push_back(run_id_arr);
  fields.push_back(arrow::field("run_id", arrow::int64()));
  arrays.push_back(precursor_id_arr);
  fields.push_back(arrow::field("precursor_id", arrow::int64()));
  arrays.push_back(exp_rt_arr);
  fields.push_back(arrow::field("exp_rt", arrow::float64()));
  arrays.push_back(precursor_charge_arr);
  fields.push_back(arrow::field("precursor_charge", arrow::int64()));
  arrays.push_back(decoy_arr);
  fields.push_back(arrow::field("decoy", arrow::boolean()));
  arrays.push_back(transition_count_arr);
  fields.push_back(arrow::field("transition_count", arrow::int64()));
  arrays.push_back(group_id_arr);
  fields.push_back(arrow::field("group_id", arrow::utf8()));

  for (size_t i = 0; i < all_ms2_cols.size(); ++i)
  {
    std::shared_ptr<arrow::Array> a;
    ParquetFile::appendOrThrow(ms2_builders[i]->Finish(&a), all_ms2_cols[i]);
    arrays.push_back(a);
    fields.push_back(arrow::field(all_ms2_cols[i], arrow::float64()));
  }
  for (size_t i = 0; i < all_ms1_cols.size(); ++i)
  {
    std::shared_ptr<arrow::Array> a;
    ParquetFile::appendOrThrow(ms1_builders[i]->Finish(&a), all_ms1_cols[i]);
    arrays.push_back(a);
    fields.push_back(arrow::field(all_ms1_cols[i], arrow::float64()));
  }

  auto schema = std::make_shared<arrow::Schema>(fields);
  auto table = arrow::Table::Make(schema, arrays);
  ParquetFile::writeTable(table, output_path);
#else
  (void)oswpq_dir; (void)output_path; (void)level;
  throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
#endif
}
