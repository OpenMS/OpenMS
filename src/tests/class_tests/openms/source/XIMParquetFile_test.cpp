// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/XIMParquetFile.h>
#include <OpenMS/FORMAT/ParquetFilter.h>

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>

using namespace OpenMS;
using namespace std;

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

  void appendBinaryOk_(arrow::BinaryBuilder& builder, const std::string& value)
  {
    auto status = builder.Append(reinterpret_cast<const uint8_t*>(value.data()),
                                 static_cast<int32_t>(value.size()));
    TEST_EQUAL(status.ok(), true)
  }

  std::string encodeDoubles_(std::initializer_list<double> values)
  {
    std::vector<double> buffer(values);
    return std::string(reinterpret_cast<const char*>(buffer.data()),
                       buffer.size() * sizeof(double));
  }

  ::arrow::Status writeParquetTable_(const std::shared_ptr<arrow::Table>& table, const std::string& filename)
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

START_TEST(XIMParquetFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(void load(std::vector<XIMMobilogram>& output) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));
  std::vector<XIMMobilogram> mobilograms;
  xim.load(mobilograms);

  TEST_EQUAL(mobilograms.empty(), false)
  TEST_EQUAL(mobilograms[0].mobility.empty(), false)
  TEST_EQUAL(mobilograms[0].intensity.empty(), false)

  for (const auto& m : mobilograms)
  {
    TEST_EQUAL(m.mobility.size(), m.intensity.size())
  }
}
END_SECTION

START_SECTION(void getMobilograms(std::vector<XIMMobilogram>&, Int64, Int64, const std::string&, Int64, Int64, Int64, Int64, const std::string&, Int64, double, const std::string&) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));

  std::vector<XIMMobilogram> all;
  xim.getMobilograms(all);
  TEST_EQUAL(all.empty(), false)

  std::vector<XIMMobilogram> none;
  xim.getMobilograms(none, -1, -1, "", -1, -1, -1, -1, "", -1, -1.0, "feature_id=99999");
  TEST_EQUAL(none.size(), 0)

  std::vector<XIMMobilogram> filtered_string;
  xim.getMobilograms(filtered_string, -1, -1, "", -1, -1, -1, all[0].run_id, "", -1, -1.0, "");
  TEST_EQUAL(filtered_string.empty(), false)
  TEST_EQUAL(filtered_string.size() <= all.size(), true)

  std::vector<XIMMobilogram> filtered_args;
  xim.getMobilograms(filtered_args, -1, -1, "", -1, -1, -1, all[0].run_id, "", -1, -1.0, "");
  TEST_EQUAL(filtered_args.empty(), false)
  TEST_EQUAL(filtered_args.size() <= all.size(), true)

  ParquetFilter typed_filter;
  typed_filter.eq("RUN_ID", all[0].run_id);
  std::vector<XIMMobilogram> typed;
  xim.getMobilograms(typed, typed_filter);
  TEST_EQUAL(typed.empty(), false)
  TEST_EQUAL(typed.size() <= all.size(), true)

  std::vector<XIMMobilogram> invalid_filter;
  TEST_EXCEPTION(Exception::InvalidValue,
                 xim.getMobilograms(invalid_filter, -1, -1, "", -1, -1, -1, -1, "", -1, -1.0, "feature_id=="))
}
END_SECTION

START_SECTION(void getMobilograms_multi_file)
{
  const std::string file = OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim");
  XIMParquetFile xim_single(file);
  std::vector<XIMMobilogram> mobilograms_single;
  xim_single.getMobilograms(mobilograms_single);

  std::vector<std::string> files;
  files.push_back(file);
  files.push_back(file);

  XIMParquetFile xim_multi(files);
  std::vector<XIMMobilogram> mobilograms_multi;
  xim_multi.getMobilograms(mobilograms_multi);
  TEST_EQUAL(mobilograms_multi.size(), 2 * mobilograms_single.size())
}
END_SECTION

START_SECTION(void getMobilograms_large_string_and_double_feature_rt)
{
  std::string filename;
  NEW_TMP_FILE(filename);

  arrow::Int64Builder run_id_builder;
  arrow::Int64Builder precursor_id_builder;
  arrow::Int64Builder feature_id_builder;
  arrow::DoubleBuilder feature_rt_builder;
  arrow::LargeStringBuilder source_file_builder;
  arrow::LargeStringBuilder mobilogram_type_builder;
  arrow::LargeStringBuilder modified_sequence_builder;
  arrow::BinaryBuilder mobility_data_builder;
  arrow::BinaryBuilder intensity_data_builder;
  arrow::Int64Builder mobility_compression_builder;
  arrow::Int64Builder intensity_compression_builder;

  appendOk_(run_id_builder, static_cast<int64_t>(11));
  appendOk_(precursor_id_builder, static_cast<int64_t>(101));
  appendOk_(feature_id_builder, static_cast<int64_t>(202));
  appendOk_(feature_rt_builder, 12.5);
  appendOk_(source_file_builder, std::string("large_mobi.raw"));
  appendOk_(mobilogram_type_builder, std::string("precursor"));
  appendOk_(modified_sequence_builder, std::string("PEPTIDE"));
  appendBinaryOk_(mobility_data_builder, encodeDoubles_({0.8, 0.9}));
  appendBinaryOk_(intensity_data_builder, encodeDoubles_({10.0, 20.0}));
  appendOk_(mobility_compression_builder, static_cast<int64_t>(0));
  appendOk_(intensity_compression_builder, static_cast<int64_t>(0));

  auto table = arrow::Table::Make(
    arrow::schema({
      arrow::field(XIMSchema::RUN_ID, arrow::int64()),
      arrow::field(XIMSchema::PRECURSOR_ID, arrow::int64()),
      arrow::field(XIMSchema::FEATURE_ID, arrow::int64()),
      arrow::field(XIMSchema::FEATURE_RT, arrow::float64()),
      arrow::field(XIMSchema::SOURCE_FILE, arrow::large_utf8()),
      arrow::field(XIMSchema::MOBILOGRAM_TYPE, arrow::large_utf8()),
      arrow::field(XIMSchema::MODIFIED_SEQUENCE, arrow::large_utf8()),
      arrow::field(XIMSchema::MOBILITY_DATA, arrow::binary()),
      arrow::field(XIMSchema::INTENSITY_DATA, arrow::binary()),
      arrow::field(XIMSchema::MOBILITY_COMPRESSION, arrow::int64()),
      arrow::field(XIMSchema::INTENSITY_COMPRESSION, arrow::int64())
    }),
    {
      finishArray_(run_id_builder),
      finishArray_(precursor_id_builder),
      finishArray_(feature_id_builder),
      finishArray_(feature_rt_builder),
      finishArray_(source_file_builder),
      finishArray_(mobilogram_type_builder),
      finishArray_(modified_sequence_builder),
      finishArray_(mobility_data_builder),
      finishArray_(intensity_data_builder),
      finishArray_(mobility_compression_builder),
      finishArray_(intensity_compression_builder)
    });

  TEST_EQUAL(writeParquetTable_(table, filename).ok(), true)

  XIMParquetFile xim(filename);
  std::vector<XIMMobilogram> mobilograms;
  xim.getMobilograms(mobilograms, -1, -1, "", -1, -1, -1, -1, "", -1, 12.5, "MOBILOGRAM_TYPE=precursor");

  TEST_EQUAL(mobilograms.size(), 1)
  TEST_STRING_EQUAL(mobilograms[0].source_file, "large_mobi.raw")
  TEST_STRING_EQUAL(mobilograms[0].mobilogram_type, "precursor")
  TEST_REAL_SIMILAR(mobilograms[0].feature_rt, 12.5)
  TEST_EQUAL(mobilograms[0].mobility.size(), 2)
  TEST_EQUAL(mobilograms[0].intensity.size(), 2)
}
END_SECTION

START_SECTION(void getMobilograms_float_feature_rt)
{
  std::string filename;
  NEW_TMP_FILE(filename);

  arrow::Int64Builder run_id_builder;
  arrow::Int64Builder precursor_id_builder;
  arrow::FloatBuilder feature_rt_builder;
  arrow::BinaryBuilder mobility_data_builder;
  arrow::BinaryBuilder intensity_data_builder;
  arrow::Int64Builder mobility_compression_builder;
  arrow::Int64Builder intensity_compression_builder;

  appendOk_(run_id_builder, static_cast<int64_t>(12));
  appendOk_(precursor_id_builder, static_cast<int64_t>(303));
  appendOk_(feature_rt_builder, 42.5f);
  appendBinaryOk_(mobility_data_builder, encodeDoubles_({1.1, 1.2}));
  appendBinaryOk_(intensity_data_builder, encodeDoubles_({30.0, 40.0}));
  appendOk_(mobility_compression_builder, static_cast<int64_t>(0));
  appendOk_(intensity_compression_builder, static_cast<int64_t>(0));

  auto table = arrow::Table::Make(
    arrow::schema({
      arrow::field(XIMSchema::RUN_ID, arrow::int64()),
      arrow::field(XIMSchema::PRECURSOR_ID, arrow::int64()),
      arrow::field(XIMSchema::FEATURE_RT, arrow::float32()),
      arrow::field(XIMSchema::MOBILITY_DATA, arrow::binary()),
      arrow::field(XIMSchema::INTENSITY_DATA, arrow::binary()),
      arrow::field(XIMSchema::MOBILITY_COMPRESSION, arrow::int64()),
      arrow::field(XIMSchema::INTENSITY_COMPRESSION, arrow::int64())
    }),
    {
      finishArray_(run_id_builder),
      finishArray_(precursor_id_builder),
      finishArray_(feature_rt_builder),
      finishArray_(mobility_data_builder),
      finishArray_(intensity_data_builder),
      finishArray_(mobility_compression_builder),
      finishArray_(intensity_compression_builder)
    });

  TEST_EQUAL(writeParquetTable_(table, filename).ok(), true)

  XIMParquetFile xim(filename);
  std::vector<XIMMobilogram> mobilograms;
  xim.getMobilograms(mobilograms, -1, -1, "", -1, -1, -1, -1, "", -1, 42.5, "");

  TEST_EQUAL(mobilograms.size(), 1)
  TEST_EQUAL(mobilograms[0].has_feature_rt, true)
  TEST_REAL_SIMILAR(mobilograms[0].feature_rt, 42.5)
}
END_SECTION

START_SECTION(void getRuns(std::vector<XIMRunInfo>& output) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));
  std::vector<XIMRunInfo> runs;
  xim.getRuns(runs);
  TEST_EQUAL(runs.empty(), false)
  TEST_NOT_EQUAL(runs[0].run_id, 0)
}
END_SECTION

START_SECTION(void load_invalid_path)
{
  TEST_EXCEPTION(Exception::FileNotFound, XIMParquetFile("no_such_file.xim"))
}
END_SECTION

START_SECTION(void getColumns(std::vector<std::string>&) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));
  std::vector<std::string> columns;
  xim.getColumns(columns);

  TEST_EQUAL(columns.empty(), false)
  // Schema must contain at least the mandatory data columns
  bool has_run_id = false;
  bool has_mobility = false;
  bool has_intensity = false;
  for (const auto& col : columns)
  {
    if (col == "RUN_ID") has_run_id = true;
    if (col == "MOBILITY_DATA") has_mobility = true;
    if (col == "INTENSITY_DATA") has_intensity = true;
  }
  TEST_EQUAL(has_run_id, true)
  TEST_EQUAL(has_mobility, true)
  TEST_EQUAL(has_intensity, true)
}
END_SECTION

START_SECTION(void getAnalytes(std::vector<XIMAnalyte>&, const std::vector<std::string>&, bool) const)
{
  XIMParquetFile xim(OPENMS_GET_TEST_DATA_PATH("XIMParquetFile_23_input.xim"));

  // Default: all analytes with nested transitions
  std::vector<XIMAnalyte> analytes;
  xim.getAnalytes(analytes);
  TEST_EQUAL(analytes.empty(), false)

  // Non-nested: each row has scalar transition fields
  std::vector<XIMAnalyte> flat;
  xim.getAnalytes(flat, {}, false);
  TEST_EQUAL(flat.empty(), false)

  // With a valid analyte column filter
  std::vector<XIMAnalyte> filtered;
  xim.getAnalytes(filtered, {"PRECURSOR_ID", "MODIFIED_SEQUENCE", "PRECURSOR_CHARGE"}, false);
  TEST_EQUAL(filtered.empty(), false)

  // Passing a non-analyte column (e.g. a signal column) must throw
  std::vector<XIMAnalyte> bad;
  TEST_EXCEPTION(Exception::InvalidValue, xim.getAnalytes(bad, {"RUN_ID"}, false))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
