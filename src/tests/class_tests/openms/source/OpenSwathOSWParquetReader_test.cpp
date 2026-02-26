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

#ifdef WITH_PARQUET
#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#endif

using namespace OpenMS;

#ifdef WITH_PARQUET
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
#endif

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
#ifdef WITH_PARQUET
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
  arrow::BooleanBuilder detecting_builder;

  appendOk_(transition_id_builder, (int64_t)1);
  appendOk_(transition_precursor_id_builder, (int64_t)1);
  appendOk_(detecting_builder, true);

  auto transition_id_array = finishArray_(transition_id_builder);
  auto transition_precursor_id_array = finishArray_(transition_precursor_id_builder);
  auto detecting_array = finishArray_(detecting_builder);

  auto transition_schema = arrow::schema({
    arrow::field("transition_id", arrow::int64()),
    arrow::field("precursor_id", arrow::int64()),
    arrow::field("detecting", arrow::boolean())
  });

  auto transitions_table = arrow::Table::Make(transition_schema,
    {transition_id_array, transition_precursor_id_array, detecting_array});

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

#else
  NOT_TESTABLE
#endif
}
END_SECTION

END_TEST
