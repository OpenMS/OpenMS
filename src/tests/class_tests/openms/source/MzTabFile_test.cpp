// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <algorithm>
///////////////////////////

using namespace OpenMS;
using namespace std;

class MzTabFile2 : public MzTabFile
{
  public:
    std::string generateMzTabPSMSectionRow2_(const MzTabPSMSectionRow& row, const vector<std::string>& optional_columns, const MzTabMetaData& meta) const
    {
      size_t n_columns = 0;
      return generateMzTabSectionRow_(row, optional_columns, meta, n_columns);
    }

    std::string generateMzTabProteinHeader2_(const MzTabProteinSectionRow& reference_row, Size n_best_search_engine_scores, const vector<std::string>& optional_columns, const MzTabMetaData& meta) const
    {
      size_t n_columns = 0;
      return generateMzTabProteinHeader_(reference_row, n_best_search_engine_scores, optional_columns, meta, n_columns);
    }

    std::string generateMzTabProteinSectionRow2_(const MzTabProteinSectionRow& row, const vector<std::string>& optional_columns, const MzTabMetaData& meta) const
    {
      size_t n_columns = 0;
      return generateMzTabSectionRow_(row, optional_columns, meta, n_columns);
    }

    std::string peptideHeader(Size runs, Size n_scores) const
    {
      size_t n_columns = 0;
      return generateMzTabPeptideHeader_(runs, 0, n_scores, 0, 0, {}, n_columns);
    }

    std::string smallMoleculeHeader(Size runs, Size n_scores) const
    {
      size_t n_columns = 0;
      return generateMzTabSmallMoleculeHeader_(runs, 0, n_scores, 0, 0, {}, n_columns);
    }

    std::string nucleicAcidHeader(Size runs, Size n_scores) const
    {
      size_t n_columns = 0;
      return generateMzTabNucleicAcidHeader_(runs, 0, n_scores, {}, n_columns);
    }

    std::string oligonucleotideHeader(Size runs, Size n_scores) const
    {
      size_t n_columns = 0;
      return generateMzTabOligonucleotideHeader_(runs, 0, n_scores, {}, n_columns);
    }

    template <typename RowT>
    std::string sectionRow(const RowT& row) const
    {
      size_t n_columns = 0;
      return generateMzTabSectionRow_(row, {}, MzTabMetaData(), n_columns);
    }
};

using ScoreRunMap = std::map<Size, std::map<Size, MzTabDouble>>;

// search_engine_score[s]_ms_run[r] = 10 * s + r, distinct for every (score, run) pair
ScoreRunMap distinctScores(Size n_scores, Size n_runs)
{
  ScoreRunMap m;
  for (Size s = 1; s <= n_scores; ++s)
  {
    for (Size r = 1; r <= n_runs; ++r)
    {
      m[s][r] = MzTabDouble(10.0 * s + r);
    }
  }
  return m;
}

START_TEST(MzTabFile, "$Id$")

/////////////////////////////////////////////////////////////

MzTabFile* ptr = nullptr;
MzTabFile* null_ptr = nullptr;

START_SECTION(MzTabFile())
{
  ptr = new MzTabFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(void load(const std::string& filename, MzTab& mzTab) )
  MzTab mzTab;
  MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_SILAC.mzTab"), mzTab);
END_SECTION

START_SECTION(void store(const std::string& filename, MzTab& mzTab) )
{
  std::vector<std::string> files_to_test;
  files_to_test.push_back("MzTabFile_SILAC.mzTab");
  files_to_test.push_back("MzTabFile_SILAC2.mzTab");
  files_to_test.push_back("MzTabFile_labelfree.mzTab");
  files_to_test.push_back("MzTabFile_iTRAQ.mzTab");
  files_to_test.push_back("MzTabFile_Cytidine.mzTab");

  for (std::vector<std::string>::const_iterator sit = files_to_test.begin(); sit != files_to_test.end(); ++sit)
  {
    // load mzTab
    MzTab mzTab;
    MzTabFile().load(OPENMS_GET_TEST_DATA_PATH(*sit), mzTab);

    // store mzTab
    std::string stored_mzTab;
    NEW_TMP_FILE(stored_mzTab)
    MzTabFile().store(stored_mzTab, mzTab);

    // compare original and stored mzTab (discarding row order and spaces)
    TextFile file1;
    TextFile file2;
    file1.load(stored_mzTab);
    file2.load(OPENMS_GET_TEST_DATA_PATH(*sit));
    std::sort(file1.begin(), file1.end());
    std::sort(file2.begin(), file2.end());

    for (TextFile::Iterator it = file1.begin(); it != file1.end(); ++it)
    {
      StringUtils::substitute(*it, " ", "");
    }

    for (TextFile::Iterator it = file2.begin(); it != file2.end(); ++it)
    {
      StringUtils::substitute(*it, " ", "");
    }

    std::string tmpfile1;
    std::string tmpfile2;
    NEW_TMP_FILE(tmpfile1)
    NEW_TMP_FILE(tmpfile2)
    file1.store(tmpfile1);
    file2.store(tmpfile2);
    TEST_FILE_SIMILAR(tmpfile1.c_str(), tmpfile2.c_str())
  }
}
END_SECTION

START_SECTION(~MzTabFile())
{
  delete ptr;
}
END_SECTION

START_SECTION(generateMzTabPSMSectionRow_(const MzTabPSMSectionRow& row, const vector<std::string>& optional_columns) const)
{
  MzTabFile2 mzTab;
  MzTabPSMSectionRow row;
  MzTabOptionalColumnEntry e;
  MzTabString s;
  
  row.sequence.fromCellString("NDYKAPPQPAPGK");
  row.PSM_ID.fromCellString("38");
  row.accession.fromCellString("IPI:B1");
  row.unique.fromCellString("1");
  row.database.fromCellString("null");
  row.database_version.fromCellString("null");
  row.search_engine.fromCellString("[, , Percolator, ]");
  row.search_engine_score[0].fromCellString("51.9678841193106");
  
  e.first = "Percolator_score";
  s.fromCellString("0.359083");
  e.second = s;
  row.opt_.push_back(e);
  
  e.first = "Percolator_qvalue";
  s.fromCellString("0.00649874");  
  e.second = s;
  row.opt_.push_back(e);
  
  e.first = "Percolator_PEP";
  s.fromCellString("0.0420992");  
  e.second = s;
  row.opt_.push_back(e);
  
  e.first = "search_engine_sequence";
  s.fromCellString("NDYKAPPQPAPGK");  
  e.second = s;
  row.opt_.push_back(e);
  
  // Tests ///////////////////////////////  
  vector<std::string> optional_columns;
  optional_columns.push_back("Percolator_score");
  optional_columns.push_back("Percolator_qvalue");
  optional_columns.push_back("EMPTY");
  optional_columns.push_back("Percolator_PEP");
  optional_columns.push_back("search_engine_sequence");
  optional_columns.push_back("AScore_1");

  MzTabMetaData m{};
    
  std::string strRow(mzTab.generateMzTabPSMSectionRow2_(row, optional_columns, m));
  std::vector<std::string> substrings;
  StringUtils::split(strRow, '\t', substrings);
  TEST_EQUAL(substrings[substrings.size() - 1],"null")
  TEST_EQUAL(substrings[substrings.size() - 2],"NDYKAPPQPAPGK")
  TEST_EQUAL(substrings[substrings.size() - 3],"0.0420992")
  TEST_EQUAL(substrings[substrings.size() - 4],"null")
}
END_SECTION

START_SECTION([EXTRA] store - every peptide identification needs its protein identification run)
{
  std::vector<ProteinIdentification> runs(1);
  runs[0].setIdentifier("search");
  PeptideIdentificationList peptides(1);
  peptides[0].setIdentifier("other");
  peptides[0].insertHit(PeptideHit(0.1, 1, 2, AASequence::fromString("PEPTIDE")));
  std::string file;
  NEW_TMP_FILE_EXT(file, ".mzTab")
  // fails before the file is opened, i.e. without leaving a truncated file
  TEST_EXCEPTION_WITH_MESSAGE(Exception::InvalidParameter, MzTabFile().store(file, runs, peptides, false, false),
    "Peptide identification has no matching protein run: 'other'. Every peptide identification needs the protein identification run (search run) with its identifier, which may have no protein hits.")
  TEST_FALSE(File::exists(file))

  ConsensusMap cmap;
  cmap.setProteinIdentifications(runs);
  cmap.getUnassignedPeptideIdentifications() = peptides;
  TEST_EXCEPTION(Exception::InvalidParameter, MzTabFile().store(file, cmap, false, false, true, false))
  TEST_FALSE(File::exists(file))
}
END_SECTION

START_SECTION([EXTRA] colunit metadata is written and read as specified by mzTab 1.0)
{
  MzTab mz_tab;
  MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_SILAC.mzTab"), mz_tab);

  MzTabMetaData md = mz_tab.getMetaData();
  md.colunit_protein.push_back("best_search_engine_score[1]=[MS, MS:1001171, Mascot:score, ]");
  md.colunit_psm.push_back("retention_time=[UO, UO:0000031, minute, ]");
  md.colunit_psm.push_back("exp_mass_to_charge=[UO, UO:0000221, dalton, ]");
  mz_tab.setMetaData(md);

  std::string stored;
  NEW_TMP_FILE(stored)
  MzTabFile().store(stored, mz_tab);

  // key and value are separate cells, and the PSM key is lower case
  TextFile stored_lines;
  stored_lines.load(stored, true);
  bool found_psm_line = false;
  for (const std::string& line : stored_lines)
  {
    if (line == "MTD\tcolunit-psm\tretention_time=[UO, UO:0000031, minute, ]") found_psm_line = true;
  }
  TEST_EQUAL(found_psm_line, true)

  // reading back does not fail on the unindexed keys and keeps every entry in order
  MzTab reloaded;
  MzTabFile().load(stored, reloaded);
  const MzTabMetaData& reloaded_md = reloaded.getMetaData();
  TEST_EQUAL(reloaded_md.colunit_protein.size(), 1)
  TEST_EQUAL(reloaded_md.colunit_psm.size(), 2)
  TEST_STRING_EQUAL(reloaded_md.colunit_psm[0], "retention_time=[UO, UO:0000031, minute, ]")
  TEST_STRING_EQUAL(reloaded_md.colunit_psm[1], "exp_mass_to_charge=[UO, UO:0000221, dalton, ]")
  TEST_EQUAL(reloaded_md.colunit_peptide.empty(), true)
}
END_SECTION

START_SECTION([EXTRA] optional PSM columns are read back)
{
  MzTab mz_tab;
  MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_SILAC2.mzTab"), mz_tab);
  MzTabPSMSectionRows psms = mz_tab.getPSMSectionRows();
  TEST_EQUAL(psms.empty(), false)
  for (Size i = 0; i < psms.size(); ++i)
  {
    MzTabString q;
    q.set(StringUtils::toStr(0.001 * (i + 1)));
    psms[i].opt_.emplace_back("opt_global_q-value", q);
  }
  mz_tab.setPSMSectionRows(psms);

  std::string stored;
  NEW_TMP_FILE(stored)
  MzTabFile().store(stored, mz_tab);

  MzTab reloaded;
  MzTabFile().load(stored, reloaded);
  const MzTabPSMSectionRows& reloaded_psms = reloaded.getPSMSectionRows();
  TEST_EQUAL(reloaded_psms.size(), psms.size())
  ABORT_IF(reloaded_psms.size() != psms.size())
  for (Size i = 0; i < reloaded_psms.size(); ++i)
  {
    TEST_EQUAL(reloaded_psms[i].opt_.size(), 1)
    ABORT_IF(reloaded_psms[i].opt_.size() != 1)
    TEST_STRING_EQUAL(reloaded_psms[i].opt_[0].first, "opt_global_q-value")
    TEST_STRING_EQUAL(reloaded_psms[i].opt_[0].second.toCellString(), StringUtils::toStr(0.001 * (i + 1)))
  }
}
END_SECTION

START_SECTION([EXTRA] protein section run-level search engine score columns line up with their header labels)
{
  // Two ms_runs and two search engine scores: search_engine_score_ms_run is
  // keyed [score_index][run_index]. Each value is distinct so a column/value
  // swap is easy to spot.
  MzTabProteinSectionRow row;
  row.search_engine_score_ms_run[1][1] = MzTabDouble(11.0); // score 1, run 1
  row.search_engine_score_ms_run[1][2] = MzTabDouble(12.0); // score 1, run 2
  row.search_engine_score_ms_run[2][1] = MzTabDouble(21.0); // score 2, run 1
  row.search_engine_score_ms_run[2][2] = MzTabDouble(22.0); // score 2, run 2

  MzTabFile2 f;
  MzTabMetaData meta;
  vector<std::string> optional_columns;

  std::string header_line = f.generateMzTabProteinHeader2_(row, 0, optional_columns, meta);
  std::string row_line = f.generateMzTabProteinSectionRow2_(row, optional_columns, meta);

  std::vector<std::string> header_cells;
  std::vector<std::string> row_cells;
  StringUtils::split(header_line, "\t", header_cells);
  StringUtils::split(row_line, "\t", row_cells);
  TEST_EQUAL(header_cells.size(), row_cells.size())
  ABORT_IF(header_cells.size() != row_cells.size())

  // every "search_engine_score[s]_ms_run[r]" column must carry exactly the
  // value that was stored under search_engine_score_ms_run[s][r]
  for (Size score = 1; score <= 2; ++score)
  {
    for (Size run = 1; run <= 2; ++run)
    {
      std::string label = "search_engine_score[" + StringUtils::toStr(score) + "]_ms_run[" + StringUtils::toStr(run) + "]";
      std::vector<std::string>::iterator pos = std::find(header_cells.begin(), header_cells.end(), label);
      ABORT_IF(pos == header_cells.end())
      Size idx = static_cast<Size>(std::distance(header_cells.begin(), pos));
      TEST_STRING_EQUAL(row_cells[idx], row.search_engine_score_ms_run[score][run].toCellString())
    }
  }
}
END_SECTION


// For each of the other four sections with run-level scores, every
// search_engine_score[s]_ms_run[r] column of the header must carry the value
// stored under search_engine_score_ms_run[s][r] (2 scores x 3 runs).
#define CHECK_SCORE_COLUMNS(HEADER_LINE, ROW_LINE, SCORES) \
  { \
    std::vector<std::string> header_cells, row_cells; \
    StringUtils::split(HEADER_LINE, "\t", header_cells); \
    StringUtils::split(ROW_LINE, "\t", row_cells); \
    for (const auto& score : SCORES) \
    { \
      for (const auto& run : score.second) \
      { \
        std::string label = "search_engine_score[" + StringUtils::toStr(score.first) + "]_ms_run[" + StringUtils::toStr(run.first) + "]"; \
        auto pos = std::find(header_cells.begin(), header_cells.end(), label); \
        TEST_EQUAL(pos != header_cells.end(), true) \
        ABORT_IF(pos == header_cells.end()) \
        Size idx = static_cast<Size>(std::distance(header_cells.begin(), pos)); \
        TEST_EQUAL(idx < row_cells.size(), true) \
        ABORT_IF(idx >= row_cells.size()) \
        TEST_STRING_EQUAL(row_cells[idx], run.second.toCellString()) \
      } \
    } \
  }

START_SECTION([EXTRA] peptide/small molecule/nucleic acid/oligonucleotide run-level score columns line up with their header labels)
{
  MzTabFile2 f;
  const ScoreRunMap scores = distinctScores(2, 3);

  MzTabPeptideSectionRow pep;
  pep.search_engine_score_ms_run = scores;
  CHECK_SCORE_COLUMNS(f.peptideHeader(3, 2), f.sectionRow(pep), scores)

  MzTabSmallMoleculeSectionRow sm;
  sm.search_engine_score_ms_run = scores;
  CHECK_SCORE_COLUMNS(f.smallMoleculeHeader(3, 2), f.sectionRow(sm), scores)

  MzTabNucleicAcidSectionRow nuc;
  nuc.search_engine_score_ms_run = scores;
  CHECK_SCORE_COLUMNS(f.nucleicAcidHeader(3, 2), f.sectionRow(nuc), scores)

  MzTabOligonucleotideSectionRow oli;
  oli.search_engine_score_ms_run = scores;
  CHECK_SCORE_COLUMNS(f.oligonucleotideHeader(3, 2), f.sectionRow(oli), scores)
}
END_SECTION

START_SECTION([EXTRA] store() and load() keep every run-level score of a file with two scores and several runs)
{
  // MzTabFile_SILAC.mzTab has 6 ms_runs and one search engine score in its
  // protein and peptide sections. Give every row a second score with a value
  // distinct for each (row, score, run), store, load and compare.
  MzTab mztab;
  MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_SILAC.mzTab"), mztab);

  MzTabProteinSectionRows prt = mztab.getProteinSectionRows();
  MzTabPeptideSectionRows pep = mztab.getPeptideSectionRows();
  TEST_EQUAL(prt.empty(), false)
  TEST_EQUAL(pep.empty(), false)
  ABORT_IF(prt.empty() || pep.empty())
  const Size n_runs = prt[0].search_engine_score_ms_run.at(1).size();
  TEST_EQUAL(n_runs, 6)

  for (Size i = 0; i != prt.size(); ++i)
  {
    for (Size r = 1; r <= n_runs; ++r)
    {
      prt[i].search_engine_score_ms_run[1][r] = MzTabDouble(1000.0 * i + 10 + r);
      prt[i].search_engine_score_ms_run[2][r] = MzTabDouble(1000.0 * i + 20 + r);
    }
  }
  for (Size i = 0; i != pep.size(); ++i)
  {
    for (Size r = 1; r <= n_runs; ++r)
    {
      pep[i].search_engine_score_ms_run[1][r] = MzTabDouble(1000.0 * i + 10 + r);
      pep[i].search_engine_score_ms_run[2][r] = MzTabDouble(1000.0 * i + 20 + r);
    }
  }
  mztab.setProteinSectionRows(prt);
  mztab.setPeptideSectionRows(pep);

  std::string stored;
  NEW_TMP_FILE(stored)
  MzTabFile().store(stored, mztab);

  MzTab reloaded;
  MzTabFile().load(stored, reloaded);
  const MzTabProteinSectionRows& prt2 = reloaded.getProteinSectionRows();
  const MzTabPeptideSectionRows& pep2 = reloaded.getPeptideSectionRows();
  TEST_EQUAL(prt2.size(), prt.size())
  TEST_EQUAL(pep2.size(), pep.size())
  ABORT_IF(prt2.size() != prt.size() || pep2.size() != pep.size())

  Size n_wrong = 0;
  for (Size i = 0; i != prt.size(); ++i)
  {
    for (const auto& score : prt[i].search_engine_score_ms_run)
    {
      for (const auto& run : score.second)
      {
        auto s_it = prt2[i].search_engine_score_ms_run.find(score.first);
        if (s_it == prt2[i].search_engine_score_ms_run.end()) { ++n_wrong; continue; }
        auto r_it = s_it->second.find(run.first);
        if (r_it == s_it->second.end() || r_it->second.toCellString() != run.second.toCellString()) { ++n_wrong; }
      }
    }
  }
  TEST_EQUAL(n_wrong, 0)

  n_wrong = 0;
  for (Size i = 0; i != pep.size(); ++i)
  {
    for (const auto& score : pep[i].search_engine_score_ms_run)
    {
      for (const auto& run : score.second)
      {
        auto s_it = pep2[i].search_engine_score_ms_run.find(score.first);
        if (s_it == pep2[i].search_engine_score_ms_run.end()) { ++n_wrong; continue; }
        auto r_it = s_it->second.find(run.first);
        if (r_it == s_it->second.end() || r_it->second.toCellString() != run.second.toCellString()) { ++n_wrong; }
      }
    }
  }
  TEST_EQUAL(n_wrong, 0)
}
END_SECTION

#undef CHECK_SCORE_COLUMNS

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
