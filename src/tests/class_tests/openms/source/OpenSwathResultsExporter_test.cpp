// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathResultsExporter.h>
///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportData.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportConfig.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <fstream>

using namespace OpenMS;
using namespace std;

namespace
{
  OpenSwathExportRow mkRow(Int64 feature_id, const std::string& tg, const std::string& seq)
  {
    OpenSwathExportRow r;
    r.run_id = 99;
    r.filename = "run.mzML";
    r.run_name = "runA";
    r.feature_id = feature_id;
    r.transition_group_id = tg;
    r.sequence = seq;
    r.full_peptide_name = seq;
    r.protein_name = "PROT1";
    r.charge = 2;
    r.mz = 450.7;
    r.rt = 120.5;
    r.intensity = 1234.5;
    return r;
  }
}

START_TEST(OpenSwathResultsExporter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static void write(const std::string& filename, const std::vector<OpenSwathExportRow>& rows, const OpenSwathResultsExportConfig& config)))
{
  std::vector<OpenSwathExportRow> rows = {
    mkRow(7777, "TG1", "PEPTIDEK"),
    mkRow(8888, "TG2", "SAMPLER"),
  };

  // ----- TSV output -----
  {
    OpenSwathResultsExportConfig cfg;
    cfg.format = OpenSwathExportFileFormat::TSV;

    std::string fn;
    NEW_TMP_FILE(fn)
    OpenSwathResultsExporter::write(fn, rows, cfg);

    std::ifstream is(fn.c_str());
    std::string header, d0, d1;
    std::getline(is, header);
    std::getline(is, d0);
    std::getline(is, d1);

    StringList hcols, r0, r1;
    StringUtils::split(header, '\t', hcols);
    StringUtils::split(d0, '\t', r0);
    StringUtils::split(d1, '\t', r1);

    // header column layout
    TEST_STRING_EQUAL(hcols[0], "run_id")
    TEST_STRING_EQUAL(hcols[3], "feature_id")
    TEST_STRING_EQUAL(hcols[5], "transition_group_id")
    TEST_STRING_EQUAL(hcols[8], "Sequence")
    TEST_STRING_EQUAL(hcols[14], "RT")

    // first data row maps to the set fields
    TEST_STRING_EQUAL(r0[0], "99")        // run_id
    TEST_STRING_EQUAL(r0[3], "7777")      // feature_id
    TEST_STRING_EQUAL(r0[5], "TG1")       // transition_group_id
    TEST_STRING_EQUAL(r0[8], "PEPTIDEK")  // Sequence
    TEST_STRING_EQUAL(r0[14], "120.5")    // RT

    // second data row
    TEST_STRING_EQUAL(r1[3], "8888")
    TEST_STRING_EQUAL(r1[8], "SAMPLER")
  }

  // ----- Parquet output (same column names, read back via ParquetFile) -----
  {
    OpenSwathResultsExportConfig cfg;
    cfg.format = OpenSwathExportFileFormat::Parquet;

    std::string fn;
    NEW_TMP_FILE(fn)
    fn += ".parquet";
    OpenSwathResultsExporter::write(fn, rows, cfg);

    auto tbl = ParquetFile::readTable(fn);
    TEST_NOT_EQUAL(tbl, nullptr)
    TEST_EQUAL(tbl->num_rows(), 2)

    auto fid = ParquetFile::getColumn(tbl, "feature_id");
    TEST_EQUAL(ParquetFile::getInt64(fid, 0, -1, false), 7777)
    TEST_EQUAL(ParquetFile::getInt64(fid, 1, -1, false), 8888)
    TEST_EQUAL(ParquetFile::getString(ParquetFile::getColumn(tbl, "Sequence"), 0), "PEPTIDEK")
    TEST_EQUAL(ParquetFile::getString(ParquetFile::getColumn(tbl, "transition_group_id"), 1), "TG2")
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
