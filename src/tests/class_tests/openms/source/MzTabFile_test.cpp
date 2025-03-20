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
///////////////////////////

using namespace OpenMS;
using namespace std;

class MzTabFile2 : public MzTabFile
{
  public:
    String generateMzTabPSMSectionRow2_(const MzTabPSMSectionRow& row, const vector<String>& optional_columns, const MzTabMetaData& meta) const
    {
      size_t n_columns = 0;
      return generateMzTabSectionRow_(row, optional_columns, meta, n_columns);
    }
};

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

START_SECTION(void load(const String& filename, MzTab& mzTab) )
  MzTab mzTab;
  MzTabFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabFile_SILAC.mzTab"), mzTab);
END_SECTION

START_SECTION(void store(const String& filename, MzTab& mzTab) )
{
  std::vector<String> files_to_test;
  files_to_test.push_back("MzTabFile_SILAC.mzTab");
  files_to_test.push_back("MzTabFile_SILAC2.mzTab");
  files_to_test.push_back("MzTabFile_labelfree.mzTab");
  files_to_test.push_back("MzTabFile_iTRAQ.mzTab");
  files_to_test.push_back("MzTabFile_Cytidine.mzTab");

  for (std::vector<String>::const_iterator sit = files_to_test.begin(); sit != files_to_test.end(); ++sit)
  {
    // load mzTab
    MzTab mzTab;
    MzTabFile().load(OPENMS_GET_TEST_DATA_PATH(*sit), mzTab);

    // store mzTab
    String stored_mzTab;
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
      it->substitute(" ","");
    }

    for (TextFile::Iterator it = file2.begin(); it != file2.end(); ++it)
    {
      it->substitute(" ","");
    }

    String tmpfile1;
    String tmpfile2;
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

START_SECTION(generateMzTabPSMSectionRow_(const MzTabPSMSectionRow& row, const vector<String>& optional_columns) const)
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
  vector<String> optional_columns;
  optional_columns.push_back("Percolator_score");
  optional_columns.push_back("Percolator_qvalue");
  optional_columns.push_back("EMPTY");
  optional_columns.push_back("Percolator_PEP");
  optional_columns.push_back("search_engine_sequence");
  optional_columns.push_back("AScore_1");

  MzTabMetaData m{};
    
  String strRow(mzTab.generateMzTabPSMSectionRow2_(row, optional_columns, m));
  std::vector<String> substrings;
  strRow.split('\t', substrings);
  TEST_EQUAL(substrings[substrings.size() - 1],"null")
  TEST_EQUAL(substrings[substrings.size() - 2],"NDYKAPPQPAPGK")
  TEST_EQUAL(substrings[substrings.size() - 3],"0.0420992")
  TEST_EQUAL(substrings[substrings.size() - 4],"null")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
