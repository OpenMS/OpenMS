// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetReader.h>
#include <OpenMS/SYSTEM/File.h>
#include <iostream>
#include <cmath>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

using namespace OpenMS;

namespace
{
  template <typename BuilderT, typename ValueT>
  void appendOk_(BuilderT& builder, const ValueT& value)
  {
    TEST_EQUAL(builder.Append(value).ok(), true)
  }

  template <typename BuilderT>
  std::shared_ptr<arrow::Array> finishArray_(BuilderT& builder)
  {
    auto result = builder.Finish();
    TEST_EQUAL(result.ok(), true)
    if (result.ok())
    {
      return result.ValueOrDie();
    }
    return nullptr;
  }

  ::arrow::Status writeParquetTable_(const std::shared_ptr<arrow::Table>& table, const String& filename)
  {
    auto outfile_result = arrow::io::FileOutputStream::Open(std::string(filename));
    if (!outfile_result.ok())
    {
      return outfile_result.status();
    }
    auto outfile = outfile_result.ValueOrDie();
    return parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), outfile, 1024);
  }
}

START_TEST(OpenSwathOSWParquetReader, "$Id$")

START_SECTION(OpenSwathOSWParquetReader())
{
  OpenSwathOSWParquetReader* ptr = nullptr;
  ptr = new OpenSwathOSWParquetReader();
  TEST_NOT_EQUAL(ptr, nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(void load(const String& oswpq_dir))
{
  // Create a minimal .oswpq directory with library, runs and features
  File::TempDir tmp_dir;
  const String base_dir = tmp_dir.getPath() + "/test.oswpq";
  const String library_dir = base_dir + "/library";
  const String runs_dir = base_dir + "/runs";
  const String run_subdir = runs_dir + "/run_id=1";
  File::makeDir(base_dir);
  File::makeDir(library_dir);
  File::makeDir(runs_dir);
  File::makeDir(run_subdir);

  // precursors.parquet (one precursor)
  arrow::Int64Builder precursor_id_builder;
  arrow::DoubleBuilder precursor_mz_builder;
  arrow::Int32Builder precursor_charge_builder;
  arrow::BooleanBuilder decoy_builder;

  appendOk_(precursor_id_builder, (int64_t)1);
  appendOk_(precursor_mz_builder, 500.0);
  appendOk_(precursor_charge_builder, (int32_t)2);
  appendOk_(decoy_builder, false);

  auto precursor_id_array = finishArray_(precursor_id_builder);
  auto precursor_mz_array = finishArray_(precursor_mz_builder);
  auto precursor_charge_array = finishArray_(precursor_charge_builder);
  auto decoy_array = finishArray_(decoy_builder);

  auto precursor_schema = arrow::schema({
    arrow::field("precursor_id", arrow::int64()),
    arrow::field("precursor_mz", arrow::float64()),
    arrow::field("charge", arrow::int32()),
    arrow::field("decoy", arrow::boolean())
  });

  auto precursors_table = arrow::Table::Make(precursor_schema,
    {precursor_id_array, precursor_mz_array, precursor_charge_array, decoy_array});

  TEST_EQUAL(writeParquetTable_(precursors_table, library_dir + "/precursors.parquet").ok(), true)

  // transitions.parquet (one transition detecting true)
  arrow::Int64Builder transition_id_builder;
  arrow::Int64Builder transition_precursor_id_builder;
  arrow::Int32Builder transition_charge_builder;
  arrow::BooleanBuilder transition_decoy_builder;
  arrow::BooleanBuilder detecting_builder;

  appendOk_(transition_id_builder, (int64_t)1);
  appendOk_(transition_precursor_id_builder, (int64_t)1);
  appendOk_(transition_charge_builder, (int32_t)1);
  appendOk_(transition_decoy_builder, false);
  appendOk_(detecting_builder, true);

  auto transition_id_array = finishArray_(transition_id_builder);
  auto transition_precursor_id_array = finishArray_(transition_precursor_id_builder);
  auto transition_charge_array = finishArray_(transition_charge_builder);
  auto transition_decoy_array = finishArray_(transition_decoy_builder);
  auto detecting_array = finishArray_(detecting_builder);

  auto transition_schema = arrow::schema({
    arrow::field("transition_id", arrow::int64()),
    arrow::field("precursor_id", arrow::int64()),
    arrow::field("charge", arrow::int32()),
    arrow::field("decoy", arrow::boolean()),
    arrow::field("detecting", arrow::boolean())
  });

  auto transitions_table = arrow::Table::Make(transition_schema,
    {transition_id_array, transition_precursor_id_array, transition_charge_array, transition_decoy_array, detecting_array});

  TEST_EQUAL(writeParquetTable_(transitions_table, library_dir + "/transitions.parquet").ok(), true)

  // runs/runs.parquet (one run with id 1)
  arrow::Int64Builder run_id_builder;
  appendOk_(run_id_builder, (int64_t)1);
  auto run_id_array = finishArray_(run_id_builder);
  auto runs_schema = arrow::schema({ arrow::field("run_id", arrow::int64()) });
  auto runs_table = arrow::Table::Make(runs_schema, {run_id_array});
  TEST_EQUAL(writeParquetTable_(runs_table, runs_dir + "/runs.parquet").ok(), true)

  // features.parquet for run_id=1 (two features)
  arrow::Int64Builder feature_id_builder;
  arrow::Int64Builder feature_run_id_builder;
  arrow::Int64Builder feature_precursor_id_builder;
  arrow::DoubleBuilder exp_rt_builder;

  appendOk_(feature_id_builder, (int64_t)1);
  appendOk_(feature_run_id_builder, (int64_t)1);
  appendOk_(feature_precursor_id_builder, (int64_t)1);
  appendOk_(exp_rt_builder, 10.0);

  appendOk_(feature_id_builder, (int64_t)2);
  appendOk_(feature_run_id_builder, (int64_t)1);
  appendOk_(feature_precursor_id_builder, (int64_t)1);
  appendOk_(exp_rt_builder, 11.0);

  auto feature_id_array = finishArray_(feature_id_builder);
  auto feature_run_id_array = finishArray_(feature_run_id_builder);
  auto feature_precursor_id_array = finishArray_(feature_precursor_id_builder);
  auto exp_rt_array = finishArray_(exp_rt_builder);

  auto features_schema = arrow::schema({
    arrow::field("feature_id", arrow::int64()),
    arrow::field("run_id", arrow::int64()),
    arrow::field("precursor_id", arrow::int64()),
    arrow::field("exp_rt", arrow::float64())
  });

  auto features_table = arrow::Table::Make(features_schema,
    {feature_id_array, feature_run_id_array, feature_precursor_id_array, exp_rt_array});

  TEST_EQUAL(writeParquetTable_(features_table, run_subdir + "/features.parquet").ok(), true)

  // feature_transition.parquet for run_id=1 (each feature linked to transition 1)
  arrow::Int64Builder ft_feature_id_builder;
  arrow::Int64Builder ft_transition_id_builder;
  arrow::DoubleBuilder ft_area_builder;

  appendOk_(ft_feature_id_builder, (int64_t)1);
  appendOk_(ft_transition_id_builder, (int64_t)1);
  appendOk_(ft_area_builder, 100.0);

  appendOk_(ft_feature_id_builder, (int64_t)2);
  appendOk_(ft_transition_id_builder, (int64_t)1);
  appendOk_(ft_area_builder, 110.0);

  auto ft_feature_id_array = finishArray_(ft_feature_id_builder);
  auto ft_transition_id_array = finishArray_(ft_transition_id_builder);
  auto ft_area_array = finishArray_(ft_area_builder);

  auto ft_schema = arrow::schema({
    arrow::field("feature_id", arrow::int64()),
    arrow::field("transition_id", arrow::int64()),
    arrow::field("area_intensity", arrow::float64())
  });

  auto ft_table = arrow::Table::Make(ft_schema, {ft_feature_id_array, ft_transition_id_array, ft_area_array});
  TEST_EQUAL(writeParquetTable_(ft_table, run_subdir + "/feature_transition.parquet").ok(), true)

  // Now read with the reader
  OpenSwathOSWParquetReader reader(base_dir);
  // rows() should be populated by load()
  const auto& rows = reader.rows();
  TEST_EQUAL(rows.size(), 2)
  // Check some values
  TEST_EQUAL(rows[0].feature_id == 1 || rows[1].feature_id == 1, true)

  // fetchPeakGroupFeatures should return two entries as well
  auto result = reader.fetchPeakGroupFeatures(base_dir);
  TEST_EQUAL(result.feature_id.size(), 2)
  TEST_EQUAL(result.precursor_charge.size(), 2)
  TEST_EQUAL(result.precursor_charge[0], 2)

  // fetchTransitionFeatures should return two transition rows as well
  auto tresult = reader.fetchTransitionFeatures(base_dir);
  TEST_EQUAL(tresult.feature_id.size(), 2)
  TEST_EQUAL(tresult.transition_id.size(), 2)
  TEST_EQUAL(tresult.product_charge.size(), 2)
  TEST_EQUAL(tresult.product_charge[0], 1)
  TEST_EQUAL(tresult.precursor_charge[0], 2)
  TEST_EQUAL(tresult.group_id[0], String("1_1_1_1"))

  // Test fetchUnscoredData on the synthetic temp .oswpq we just created
  auto uresult = reader.fetchUnscoredData(base_dir);
  // We created two features, so many vectors should be length 2
  TEST_EQUAL(uresult.id.size(), 2)
  TEST_EQUAL(uresult.RT.size(), 2)
  // Print a short summary for inspection
  std::cout << "Unscored synthetic rows: " << uresult.id.size() << "\n";
  for (size_t i = 0; i < uresult.id.size(); ++i)
  {
    std::cout << "urow " << i << ": id=" << uresult.id[i]
              << " transition_group_id=" << uresult.transition_group_id[i]
              << " RT=" << uresult.RT[i]
              << " Intensity=" << uresult.Intensity[i]
              << " Charge=" << uresult.Charge[i]
              << "\n";
  }

    
  // Test using test file
  const String real_path = OPENMS_GET_TEST_DATA_PATH("OpenSwathWorkflow_tworuns_1_17.output.oswpq");

  // Read once 
  OpenSwathOSWParquetReader real_reader(real_path);

  // Peak-group (peak-level) checks
  auto pf = real_reader.fetchPeakGroupFeatures(real_path);
  TEST_NOT_EQUAL(pf.feature_id.size(), 0)
  TEST_EQUAL(pf.feature_id.size(), pf.exp_rt.size())
  // spot-check first few rows for sanity
  for (size_t i = 0; i < std::min<size_t>(pf.feature_id.size(), 5); ++i)
  {
      std::cout << "Peak group feature " << i << ": feature_id=" << pf.feature_id[i]
              << " exp_rt=" << pf.exp_rt[i]
              << " precursor_charge=" << pf.precursor_charge[i]
              << "\n";
      TEST_EQUAL(pf.feature_id[i] > 0, true)
      // exp_rt should be a finite number
      TEST_EQUAL(std::isfinite(pf.exp_rt[i]), true)
  }
  // Transition-level checks
  auto tf = real_reader.fetchTransitionFeatures(real_path);
  TEST_NOT_EQUAL(tf.feature_id.size(), 0)
  TEST_EQUAL(tf.feature_id.size(), tf.area_intensity.size())
  for (size_t i = 0; i < std::min<size_t>(tf.feature_id.size(), 5); ++i)
  {
      std::cout << "Transition feature " << i << ": feature_id=" << tf.feature_id[i]
              << " transition_id=" << tf.transition_id[i]
              << " precursor_charge=" << tf.precursor_charge[i]
              << " product_charge=" << tf.product_charge[i]
              << " area_intensity=" << tf.area_intensity[i]
              << "\n";
      // area_intensity should be finite and non-negative
      TEST_EQUAL(std::isfinite(tf.area_intensity[i]), true)
      TEST_EQUAL(tf.area_intensity[i] >= 0.0, true)
  }
  // Unscored data (pyprophet-style) checks
  auto uf = real_reader.fetchUnscoredData(real_path);
  TEST_NOT_EQUAL(uf.id.size(), 0)
  TEST_EQUAL(uf.id.size(), uf.RT.size())

  if (uf.leftWidth.size() == uf.id.size())
  {
      for (size_t i = 0; i < std::min<size_t>(uf.id.size(), 5); ++i)
      {
        std::cout << "Unscored feature " << i << ": id=" << uf.id[i]
                << " transition_group_id=" << uf.transition_group_id[i]
                << " RT=" << uf.RT[i]
                << " Intensity=" << uf.Intensity[i]
                << " Charge=" << uf.Charge[i]
                << " leftWidth=" << uf.leftWidth[i]
                << " rightWidth=" << uf.rightWidth[i]
                << "\n";
       TEST_EQUAL(std::isfinite(uf.leftWidth[i]), true)
       TEST_EQUAL(std::isfinite(uf.rightWidth[i]), true)
      }
  }

  // Basic consistency: at least one of the result tables should contain more than one row
  TEST_EQUAL((pf.feature_id.size() > 1) || (tf.feature_id.size() > 1) || (uf.id.size() > 1), true)
}
END_SECTION

END_TEST
