// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathMatrixExporter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathParquetExporter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathResultsExporter.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/SYSTEM/File.h>

#include <arrow/api.h>
#include <arrow/io/file.h>
#include <parquet/arrow/reader.h>

#include <fstream>

using namespace OpenMS;

namespace
{
  void copySharedExportFixture_(const String& filename)
  {
    File::remove(filename);
    TEST_EQUAL(File::copy(OPENMS_GET_TEST_DATA_PATH("PyProphet_inference_test.osw"), filename), true)
  }

  std::shared_ptr<arrow::Table> readParquetTable_(const String& filename)
  {
    auto infile_result = arrow::io::ReadableFile::Open(std::string(filename));
    TEST_EQUAL(infile_result.ok(), true)
    if (!infile_result.ok())
    {
      return {};
    }
    const auto& infile = *infile_result;

    auto reader_result = parquet::arrow::OpenFile(infile, arrow::default_memory_pool());
    TEST_EQUAL(reader_result.ok(), true)
    if (!reader_result.ok())
    {
      return {};
    }

    std::unique_ptr<parquet::arrow::FileReader> reader = std::move(reader_result.ValueOrDie());
    std::shared_ptr<arrow::Table> table;
    auto read_status = reader->ReadTable(&table);
    TEST_EQUAL(read_status.ok(), true)
    if (!read_status.ok())
    {
      return {};
    }
    return table;
  }

  void dropTable_(const String& filename, const String& table_name)
  {
    SqliteConnector conn(filename, SqliteConnector::SqlOpenMode::READWRITE);
    conn.executeStatement("DROP TABLE " + table_name + ";");
  }
} // namespace

START_TEST(OpenSwathExport, "$Id$")

START_SECTION(OSW-backed OpenSWATH export readers and writers)
{
  String tmp_osw;
  NEW_TMP_FILE(tmp_osw);
  copySharedExportFixture_(tmp_osw);

  OSWFile osw(tmp_osw);

  OpenSwathExportFilterConfig filter_config;
  const auto result_rows = osw.readOpenSwathExportRows(filter_config);
  TEST_EQUAL(result_rows.empty(), false)

  OpenSwathExportFilterConfig empty_filter_config;
  empty_filter_config.max_rs_peakgroup_qvalue = 0.0;
  const auto empty_result_rows = osw.readOpenSwathExportRows(empty_filter_config);
  TEST_EQUAL(empty_result_rows.empty(), true)

  OpenSwathMatrixExportConfig matrix_config;
  matrix_config.level = OpenSwathMatrixLevel::Peptide;
  const auto matrix = OpenSwathMatrixExporter::buildMatrix(result_rows, matrix_config);
  TEST_EQUAL(matrix.identifier_rows.empty(), false)
  TEST_EQUAL(matrix.sample_column_names.empty(), false)

  String results_tsv;
  NEW_TMP_FILE(results_tsv);
  OpenSwathResultsExporter::write(results_tsv, result_rows, {});
  TEST_EQUAL(File::exists(results_tsv), true)
  {
    std::ifstream is(results_tsv.c_str());
    TEST_EQUAL(static_cast<bool>(is), true)
    std::string header;
    std::getline(is, header);
    TEST_EQUAL(header.find("GeneName") != std::string::npos, true)
    TEST_EQUAL(header.find("m_score_gene_global") != std::string::npos, true)
  }

  String empty_results_tsv;
  NEW_TMP_FILE(empty_results_tsv);
  OpenSwathResultsExporter::write(empty_results_tsv, empty_result_rows, {});
  TEST_EQUAL(File::exists(empty_results_tsv), true)
  {
    std::ifstream is(empty_results_tsv.c_str());
    TEST_EQUAL(static_cast<bool>(is), true)
    std::string header;
    std::getline(is, header);
    TEST_EQUAL(header.find("Sequence") != std::string::npos, true)
    std::string first_data_line;
    TEST_EQUAL(static_cast<bool>(std::getline(is, first_data_line)), false)
  }

  String matrix_tsv;
  NEW_TMP_FILE(matrix_tsv);
  OpenSwathMatrixExporter::writeMatrix(matrix_tsv, matrix, matrix_config);
  TEST_EQUAL(File::exists(matrix_tsv), true)
  {
    std::ifstream is(matrix_tsv.c_str());
    TEST_EQUAL(static_cast<bool>(is), true)
    std::string header;
    std::getline(is, header);
    TEST_EQUAL(header.find("Sequence") != std::string::npos, true)
  }

  OpenSwathParquetExportConfig parquet_config;
  parquet_config.filters.exclude_decoys = false;
  parquet_config.include_transition_data = true;

  const auto feature_table = osw.readOpenSwathFeatureScoreTable(parquet_config);
  TEST_EQUAL(feature_table.rows.empty(), false)
  TEST_EQUAL(feature_table.feature_ms2_column_names.empty(), false)
  TEST_EQUAL(feature_table.rows.front().score_ms2_qvalue.has_value(), true)

  const auto transition_table = osw.readOpenSwathTransitionScoreTable(parquet_config);
  TEST_EQUAL(transition_table.rows.empty(), false)
  TEST_EQUAL(transition_table.feature_transition_column_names.empty(), false)

  String feature_parquet;
  NEW_TMP_FILE(feature_parquet);
  feature_parquet += ".parquet";
  OpenSwathParquetExporter::writeFeatureScores(feature_parquet, feature_table);
  TEST_EQUAL(File::exists(feature_parquet), true)

  const auto feature_arrow = readParquetTable_(feature_parquet);
  TEST_EQUAL(feature_arrow != nullptr, true)
  if (feature_arrow)
  {
    TEST_EQUAL(feature_arrow->num_rows(), static_cast<Int64>(feature_table.rows.size()))
    TEST_NOT_EQUAL(feature_arrow->GetColumnByName("SCORE_MS2_QVALUE"), nullptr)
    TEST_NOT_EQUAL(feature_arrow->GetColumnByName("FEATURE_MS2_AREA_INTENSITY"), nullptr)
    TEST_NOT_EQUAL(feature_arrow->GetColumnByName("PROTEIN_ID"), nullptr)
  }

  String transition_parquet;
  NEW_TMP_FILE(transition_parquet);
  transition_parquet += ".parquet";
  OpenSwathParquetExporter::writeTransitionScores(transition_parquet, transition_table);
  TEST_EQUAL(File::exists(transition_parquet), true)

  const auto transition_arrow = readParquetTable_(transition_parquet);
  TEST_EQUAL(transition_arrow != nullptr, true)
  if (transition_arrow)
  {
    TEST_EQUAL(transition_arrow->num_rows(), static_cast<Int64>(transition_table.rows.size()))
    TEST_NOT_EQUAL(transition_arrow->GetColumnByName("IPF_PEPTIDE_ID"), nullptr)
    TEST_NOT_EQUAL(transition_arrow->GetColumnByName("SCORE_TRANSITION_QVALUE"), nullptr)
  }
}
END_SECTION

START_SECTION(OpenSwathMatrixExporter rejects invalid top_n and malformed matrix shapes)
{
  String tmp_osw;
  NEW_TMP_FILE(tmp_osw);
  copySharedExportFixture_(tmp_osw);

  OSWFile osw(tmp_osw);
  OpenSwathExportFilterConfig filter_config;
  const auto result_rows = osw.readOpenSwathExportRows(filter_config);

  OpenSwathMatrixExportConfig invalid_top_n_config;
  invalid_top_n_config.level = OpenSwathMatrixLevel::Peptide;
  invalid_top_n_config.top_n = 0;
  TEST_EXCEPTION(Exception::Precondition, OpenSwathMatrixExporter::buildMatrix(result_rows, invalid_top_n_config))

  OpenSwathQuantMatrix malformed_matrix;
  malformed_matrix.identifier_column_names = {"Sequence"};
  malformed_matrix.sample_column_names = {"run1"};
  malformed_matrix.identifier_rows = {{"PEPTIDE"}};
  malformed_matrix.values = {};
  TEST_EXCEPTION(Exception::Precondition, OpenSwathMatrixExporter::writeMatrix("ignored.tsv", malformed_matrix, {}))

  malformed_matrix.values = {{1.0}};
  malformed_matrix.identifier_rows = {{"PEPTIDE", "EXTRA"}};
  TEST_EXCEPTION(Exception::Precondition, OpenSwathMatrixExporter::writeMatrix("ignored.tsv", malformed_matrix, {}))

  malformed_matrix.identifier_rows = {{"PEPTIDE"}};
  malformed_matrix.values = {{1.0, 2.0}};
  TEST_EXCEPTION(Exception::Precondition, OpenSwathMatrixExporter::writeMatrix("ignored.tsv", malformed_matrix, {}))
}
END_SECTION

START_SECTION(OSWFile transition parquet reader validates TRANSITION_PRECURSOR_MAPPING presence)
{
  String tmp_osw;
  NEW_TMP_FILE(tmp_osw);
  copySharedExportFixture_(tmp_osw);
  dropTable_(tmp_osw, "TRANSITION_PRECURSOR_MAPPING");

  OSWFile osw(tmp_osw);
  OpenSwathParquetExportConfig parquet_config;
  parquet_config.include_transition_data = true;
  TEST_EXCEPTION(Exception::Precondition, osw.readOpenSwathTransitionScoreTable(parquet_config))
}
END_SECTION

END_TEST
