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
#include <OpenMS/FORMAT/XICParquetFile.h>
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

START_TEST(XICParquetFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(void load(std::vector<XICChromatogram>& output) const)
{
  XICParquetFile xic(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_1_input.xic"));
  std::vector<XICChromatogram> chroms;
  xic.load(chroms);
  TEST_EQUAL(chroms.size(), 18)
  TEST_EQUAL(chroms[0].rt.empty(), false)
  TEST_EQUAL(chroms[0].intensity.empty(), false)

  // RT/intensity lengths should match for all rows
  for (const auto& c : chroms)
  {
    TEST_EQUAL(c.rt.size(), c.intensity.size())
  }
}
END_SECTION

START_SECTION(void getChromatograms(std::vector<XICChromatogram>&, Int64, Int64, const std::string&, Int64, Int64, Int64, Int64, const std::string&) const)
{
  XICParquetFile xic(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_1_input.xic"));
  std::vector<XICChromatogram> chroms_all;
  xic.getChromatograms(chroms_all);
  TEST_EQUAL(chroms_all.size(), 18)

  std::vector<XICChromatogram> chroms_none;
  xic.getChromatograms(chroms_none, -1, -1, "", -1, -1, -1, -1, "precursor_id=99999");
  TEST_EQUAL(chroms_none.size(), 0)

  std::vector<XICChromatogram> chroms_filtered;
  xic.getChromatograms(chroms_filtered, -1, -1, "", -1, -1, -1, -1, "precursor_id=2");
  TEST_EQUAL(chroms_filtered.size() > 0, true)
  TEST_EQUAL(chroms_filtered.size() <= chroms_all.size(), true)

  std::vector<XICChromatogram> chroms_filtered_eqeq;
  xic.getChromatograms(chroms_filtered_eqeq, -1, -1, "", -1, -1, -1, -1, "precursor_id==2");
  TEST_EQUAL(chroms_filtered_eqeq.size(), chroms_filtered.size())

  ParquetFilter typed_filter;
  typed_filter.eq("PRECURSOR_ID", 2);
  std::vector<XICChromatogram> chroms_typed;
  xic.getChromatograms(chroms_typed, typed_filter);
  TEST_EQUAL(chroms_typed.size() > 0, true)
  TEST_EQUAL(chroms_typed.size() <= chroms_all.size(), true)

  // Invalid filter syntax should throw a parse/bind error
  std::vector<XICChromatogram> chroms_invalid_filter;
  TEST_EXCEPTION(Exception::InvalidValue,
                 xic.getChromatograms(chroms_invalid_filter, -1, -1, "", -1, -1, -1, -1, "precursor_id=="))

  std::vector<XICChromatogram> chroms_invalid_in;
  TEST_EXCEPTION(Exception::InvalidValue,
                 xic.getChromatograms(chroms_invalid_in, -1, -1, "", -1, -1, -1, -1, "precursor_id IN []"))
}
END_SECTION

START_SECTION(void getChromatograms_multi_file)
{
  std::vector<std::string> files;
  files.emplace_back(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_1_input.xic"));
  files.emplace_back(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_2_input.xic"));

  XICParquetFile xic_multi(files);
  std::vector<XICChromatogram> chroms_multi;
  xic_multi.getChromatograms(chroms_multi);
  TEST_EQUAL(chroms_multi.size(), 36)
}
END_SECTION

START_SECTION(void getChromatograms_large_string_columns)
{
  std::string filename;
  NEW_TMP_FILE(filename);

  arrow::Int64Builder run_id_builder;
  arrow::Int64Builder precursor_id_builder;
  arrow::LargeStringBuilder source_file_builder;
  arrow::LargeStringBuilder modified_sequence_builder;
  arrow::BinaryBuilder rt_data_builder;
  arrow::BinaryBuilder intensity_data_builder;
  arrow::Int64Builder rt_compression_builder;
  arrow::Int64Builder intensity_compression_builder;

  appendOk_(run_id_builder, static_cast<int64_t>(7));
  appendOk_(precursor_id_builder, static_cast<int64_t>(42));
  appendOk_(source_file_builder, std::string("large_source.raw"));
  appendOk_(modified_sequence_builder, std::string("PEPTIDE"));
  appendBinaryOk_(rt_data_builder, encodeDoubles_({100.0, 101.0}));
  appendBinaryOk_(intensity_data_builder, encodeDoubles_({1000.0, 1001.0}));
  appendOk_(rt_compression_builder, static_cast<int64_t>(0));
  appendOk_(intensity_compression_builder, static_cast<int64_t>(0));

  auto table = arrow::Table::Make(
    arrow::schema({
      arrow::field(XICSchema::RUN_ID, arrow::int64()),
      arrow::field(XICSchema::PRECURSOR_ID, arrow::int64()),
      arrow::field(XICSchema::SOURCE_FILE, arrow::large_utf8()),
      arrow::field(XICSchema::MODIFIED_SEQUENCE, arrow::large_utf8()),
      arrow::field(XICSchema::RT_DATA, arrow::binary()),
      arrow::field(XICSchema::INTENSITY_DATA, arrow::binary()),
      arrow::field(XICSchema::RT_COMPRESSION, arrow::int64()),
      arrow::field(XICSchema::INTENSITY_COMPRESSION, arrow::int64())
    }),
    {
      finishArray_(run_id_builder),
      finishArray_(precursor_id_builder),
      finishArray_(source_file_builder),
      finishArray_(modified_sequence_builder),
      finishArray_(rt_data_builder),
      finishArray_(intensity_data_builder),
      finishArray_(rt_compression_builder),
      finishArray_(intensity_compression_builder)
    });

  TEST_EQUAL(writeParquetTable_(table, filename).ok(), true)

  XICParquetFile xic(filename);
  std::vector<XICChromatogram> chroms;
  xic.getChromatograms(chroms, -1, -1, "", -1, -1, -1, -1, "MODIFIED_SEQUENCE=PEPTIDE");

  TEST_EQUAL(chroms.size(), 1)
  TEST_STRING_EQUAL(chroms[0].source_file, "large_source.raw")
  TEST_STRING_EQUAL(chroms[0].modified_sequence, "PEPTIDE")
  TEST_EQUAL(chroms[0].rt.size(), 2)
  TEST_EQUAL(chroms[0].intensity.size(), 2)
  TEST_REAL_SIMILAR(chroms[0].rt[0], 100.0)
  TEST_REAL_SIMILAR(chroms[0].intensity[1], 1001.0)
}
END_SECTION

START_SECTION(void getRuns(std::vector<XICRunInfo>& output) const)
{
  XICParquetFile xic(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_1_input.xic"));
  std::vector<XICRunInfo> runs;
  xic.getRuns(runs);
  TEST_EQUAL(runs.size(), 1)
  TEST_NOT_EQUAL(runs[0].run_id, 0)
}
END_SECTION

START_SECTION(void getAnalytes(std::vector<XICAnalyte>& output, bool) const)
{
  XICParquetFile xic(OPENMS_GET_TEST_DATA_PATH("XICParquetFile_1_input.xic"));

  std::vector<XICAnalyte> analytes_exploded;
  std::vector<std::string> columns;
  xic.getAnalytes(analytes_exploded, columns, false);
  TEST_EQUAL(analytes_exploded.size(), 18)

  std::vector<XICAnalyte> analytes_nested;
  xic.getAnalytes(analytes_nested, columns, true);
  TEST_EQUAL(analytes_nested.size(), 6)

  Size transition_count = 0;
  for (const auto& a : analytes_nested)
  {
    transition_count += a.transition_ids.size();
  }
  TEST_EQUAL(transition_count, analytes_exploded.size())
}
END_SECTION

START_SECTION(void load_invalid_path)
{
  TEST_EXCEPTION(Exception::FileNotFound, XICParquetFile("no_such_file.xic"))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
