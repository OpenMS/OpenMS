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
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathMatrixExporter.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <fstream>

using namespace OpenMS;
using namespace std;

namespace
{
  // helper: one OpenSWATH export row with the fields the Precursor-level
  // matrix build cares about
  OpenSwathExportRow mkRow(Int64 run_id, const std::string& run_name,
                           const std::string& tg, const std::string& seq,
                           const std::string& prot, Int32 charge,
                           double intensity, double m_score, Int64 feature_id)
  {
    OpenSwathExportRow r;
    r.run_id = run_id;
    r.run_name = run_name;
    r.transition_group_id = tg;
    r.sequence = seq;
    r.full_peptide_name = seq;
    r.protein_name = prot;
    r.charge = charge;
    r.intensity = intensity;
    r.m_score = m_score;
    r.feature_id = feature_id;
    return r;
  }
}

START_TEST(OpenSwathMatrixExporter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static OpenSwathQuantMatrix buildMatrix(const std::vector<OpenSwathExportRow>& rows, const OpenSwathMatrixExportConfig& config)))
{
  // Two runs (runA=1, runB=2). TG1 is seen in both runs; in runA it has two
  // peak groups and the lower-m_score one (not the higher-intensity one) must
  // win. TG2 is only seen in runA -> its runB cell must be empty.
  std::vector<OpenSwathExportRow> rows = {
    mkRow(1, "runA", "TG1", "PEPTIDEK", "PROT1", 2, /*int*/1000, /*m*/0.01, 11),
    mkRow(1, "runA", "TG1", "PEPTIDEK", "PROT1", 2, /*int*/9999, /*m*/0.05, 12), // worse m_score
    mkRow(2, "runB", "TG1", "PEPTIDEK", "PROT1", 2, /*int*/2000, /*m*/0.02, 21),
    mkRow(1, "runA", "TG2", "SAMPLER",  "PROT2", 3, /*int*/300,  /*m*/0.03, 31),
  };

  OpenSwathMatrixExportConfig config;
  config.level = OpenSwathMatrixLevel::Precursor;
  config.normalization = OpenSwathMatrixNormalization::None;

  OpenSwathQuantMatrix m = OpenSwathMatrixExporter::buildMatrix(rows, config);

  // identifier columns are the precursor-level key
  TEST_EQUAL(m.identifier_column_names.size(), 5)
  TEST_STRING_EQUAL(m.identifier_column_names[0], "transition_group_id")
  TEST_STRING_EQUAL(m.identifier_column_names[4], "ProteinName")

  // samples are ordered by run_id
  TEST_EQUAL(m.sample_column_names.size(), 2)
  TEST_STRING_EQUAL(m.sample_column_names[0], "runA")
  TEST_STRING_EQUAL(m.sample_column_names[1], "runB")

  // one row per precursor, ordered lexicographically by the key (TG1 < TG2)
  TEST_EQUAL(m.identifier_rows.size(), 2)
  TEST_STRING_EQUAL(m.identifier_rows[0][0], "TG1")
  TEST_STRING_EQUAL(m.identifier_rows[0][1], "PEPTIDEK")
  TEST_STRING_EQUAL(m.identifier_rows[0][3], "2")          // charge
  TEST_STRING_EQUAL(m.identifier_rows[0][4], "PROT1")
  TEST_STRING_EQUAL(m.identifier_rows[1][0], "TG2")
  TEST_STRING_EQUAL(m.identifier_rows[1][4], "PROT2")

  // TG1: best peak group (m_score 0.01 -> intensity 1000) in runA, 2000 in runB
  TEST_EQUAL(m.values[0][0].has_value(), true)
  TEST_REAL_SIMILAR(m.values[0][0].value(), 1000.0)        // NOT 9999 (worse m_score)
  TEST_EQUAL(m.values[0][1].has_value(), true)
  TEST_REAL_SIMILAR(m.values[0][1].value(), 2000.0)

  // TG2: present only in runA -> runB cell is empty
  TEST_EQUAL(m.values[1][0].has_value(), true)
  TEST_REAL_SIMILAR(m.values[1][0].value(), 300.0)
  TEST_EQUAL(m.values[1][1].has_value(), false)
}
END_SECTION

START_SECTION([EXTRA] buildMatrix preconditions)
{
  OpenSwathMatrixExportConfig config;
  config.level = OpenSwathMatrixLevel::Precursor;

  // no rows -> Precondition
  std::vector<OpenSwathExportRow> empty;
  TEST_EXCEPTION(Exception::Precondition, OpenSwathMatrixExporter::buildMatrix(empty, config))

  // peptide/protein/gene level requires top_n > 0
  std::vector<OpenSwathExportRow> rows = {
    mkRow(1, "runA", "TG1", "PEPTIDEK", "PROT1", 2, 1000, 0.01, 11)
  };
  OpenSwathMatrixExportConfig bad;
  bad.level = OpenSwathMatrixLevel::Peptide;
  bad.top_n = 0;
  TEST_EXCEPTION(Exception::Precondition, OpenSwathMatrixExporter::buildMatrix(rows, bad))
}
END_SECTION

START_SECTION((static void writeMatrix(const std::string& filename, const OpenSwathQuantMatrix& matrix, const OpenSwathMatrixExportConfig& config)))
{
  std::vector<OpenSwathExportRow> rows = {
    mkRow(1, "runA", "TG1", "PEPTIDEK", "PROT1", 2, 1000, 0.01, 11),
    mkRow(2, "runB", "TG1", "PEPTIDEK", "PROT1", 2, 2000, 0.02, 21),
  };
  OpenSwathMatrixExportConfig config;
  config.level = OpenSwathMatrixLevel::Precursor;
  config.normalization = OpenSwathMatrixNormalization::None;
  config.format = OpenSwathExportFileFormat::TSV;

  OpenSwathQuantMatrix m = OpenSwathMatrixExporter::buildMatrix(rows, config);

  std::string fn;
  NEW_TMP_FILE(fn)
  OpenSwathMatrixExporter::writeMatrix(fn, m, config);

  std::ifstream is(fn.c_str());
  std::string header, data;
  std::getline(is, header);
  std::getline(is, data);

  StringList hcols, dcols;
  StringUtils::split(header, '\t', hcols);
  StringUtils::split(data, '\t', dcols);

  // header: 5 identifier columns + 2 sample columns
  TEST_EQUAL(hcols.size(), 7)
  TEST_STRING_EQUAL(hcols[0], "transition_group_id")
  TEST_STRING_EQUAL(hcols[5], "runA")
  TEST_STRING_EQUAL(hcols[6], "runB")

  // data row: precursor key + the two intensities
  TEST_EQUAL(dcols.size(), 7)
  TEST_STRING_EQUAL(dcols[0], "TG1")
  TEST_STRING_EQUAL(dcols[1], "PEPTIDEK")
  TEST_STRING_EQUAL(dcols[3], "2")
  TEST_STRING_EQUAL(dcols[4], "PROT1")
  TEST_STRING_EQUAL(dcols[5], "1000")
  TEST_STRING_EQUAL(dcols[6], "2000")

  // writing an empty matrix is rejected
  OpenSwathQuantMatrix empty_matrix;
  std::string fn2;
  NEW_TMP_FILE(fn2)
  TEST_EXCEPTION(Exception::Precondition, OpenSwathMatrixExporter::writeMatrix(fn2, empty_matrix, config))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
