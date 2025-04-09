// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/FASTAContainer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FASTAContainer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

typedef FASTAContainer<TFI_Vector> FCVec;
typedef FASTAContainer<TFI_File> FCFile;

FCVec* ptr = nullptr;
FCVec* nullPointer = nullptr;
START_SECTION(FASTAContainer())
{
  ptr = new FCVec(std::vector<FASTAFile::FASTAEntry>());
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~FASTAContainer())
{
  delete ptr;
}
END_SECTION

std::vector<FASTAFile::FASTAEntry> fev = { {"id0", "desc0", "AAAA"},{ "id1", "desc1", "BBBB" },{ "id2", "desc2", "CCCC" },{ "id3", "desc3", "DDDD" } };

START_SECTION(FASTAContainer(const String& FASTA_file))
  FCFile f(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"));
  TEST_EQUAL(f.size(), 0)
END_SECTION

START_SECTION(FASTAContainer(std::vector<FASTAFile::FASTAEntry>& data))
  FCVec fv(fev);
  TEST_EQUAL(fv.size(), 4)
END_SECTION

START_SECTION(size_t getChunkOffset() const)
  // FCFile: tested below
  FCVec fv(fev);
  TEST_EQUAL(fv.getChunkOffset(), 0)
END_SECTION

START_SECTION(bool activateCache())
  // FCFile: tested below
  FCVec fv(fev);
  TEST_EQUAL(fv.activateCache(), 1)
  TEST_EQUAL(fv.activateCache(), 0)
END_SECTION

START_SECTION(void reset())
  // FCFile: tested below
  FCVec fv(fev);
  TEST_EQUAL(fv.activateCache(), 1)
  TEST_EQUAL(fv.activateCache(), 0)
  fv.reset();
  TEST_EQUAL(fv.activateCache(), 1)
END_SECTION

START_SECTION(bool cacheChunk(int suggested_size))
  // FCFile: tested below
  FCVec fv(fev);
  TEST_EQUAL(fv.cacheChunk(333), 1)
  TEST_EQUAL(fv.cacheChunk(333), 0)
END_SECTION

START_SECTION(size_t chunkSize() const)
  // FCFile: tested below
  FCVec fv(fev);
  TEST_EQUAL(fv.chunkSize(), 4)
END_SECTION

START_SECTION(const FASTAFile::FASTAEntry& chunkAt(size_t pos) const)
  // FCFile: tested below
  FCVec fv(fev);
  FASTAFile::FASTAEntry pe = fv.chunkAt(3);
  TEST_EQUAL(pe.identifier, "id3");
END_SECTION

START_SECTION(bool readAt(FASTAFile::FASTAEntry& protein, size_t pos))
  // FCFile: tested below
  FCVec fv(fev);
  FASTAFile::FASTAEntry pe;
  TEST_EQUAL(fv.readAt(pe, 3), true);
  TEST_EQUAL(pe.identifier, "id3");
END_SECTION

START_SECTION(bool empty() const)
  FCFile f(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"));
  TEST_EQUAL(f.empty(), false)
  FCFile f2(OPENMS_GET_TEST_DATA_PATH("degenerate_cases/empty.fasta"));
  TEST_EQUAL(f2.empty(), true)
  FCVec fv(fev);
  TEST_EQUAL(fv.empty(), false);
  std::vector<FASTAFile::FASTAEntry> feve;
  FCVec fv2(feve);
  TEST_EQUAL(fv2.empty(), true);
END_SECTION

START_SECTION(size_t size() const)
  FCFile f(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"));
  TEST_EQUAL(f.cacheChunk(2), true)
  TEST_EQUAL(f.size(), 2)
  TEST_EQUAL(f.activateCache(), true)
  TEST_EQUAL(f.size(), 2)
  FASTAFile::FASTAEntry pe, pe2;
  TEST_EQUAL(f.readAt(pe, 0), true);
  pe2 = f.chunkAt(0);
  TEST_TRUE(pe == pe2)
  TEST_EQUAL(pe.description, "This is the description of the first protein")
  pe2 = f.chunkAt(1);
  TEST_EQUAL(pe == pe2, false)
  TEST_EQUAL(pe2.description, "This is the description of the second protein")

  // read next chunk, and re-read from disk again, using byte offsets
  TEST_EQUAL(f.cacheChunk(1), true)
  TEST_EQUAL(f.activateCache(), true)
  TEST_EQUAL(f.readAt(pe, 0), true); // third global entry
  TEST_EQUAL(pe.identifier, "P68509|1433F_BOVIN")
  TEST_EQUAL(pe.description, "This is the description of the first protein")
  
  // read until end
  TEST_EQUAL(f.cacheChunk(3), true)  // only 2 can be read, but thats ok
  TEST_EQUAL(f.activateCache(), true)
  TEST_EQUAL(f.chunkSize(), 2)
  pe = f.chunkAt(1);
  TEST_EQUAL(pe.description, " ##0")
  TEST_EQUAL(f.readAt(pe2, 4), true);
  TEST_TRUE(pe == pe2)
  
  // reached the end after 5 entries
  TEST_EQUAL(f.cacheChunk(3), false)
  TEST_EQUAL(f.chunkSize(), 2)
  TEST_EQUAL(f.activateCache(), false)
  TEST_EQUAL(f.chunkSize(), 0)
  TEST_EQUAL(f.cacheChunk(3), false)
  TEST_EQUAL(f.activateCache(), false)

  // read from disk again after reaching EOF, using byte offsets
  TEST_EQUAL(f.readAt(pe, 0), true);
  TEST_EQUAL(pe.identifier, "P68509|1433F_BOVIN")
  TEST_EQUAL(pe.description, "This is the description of the first protein")
  TEST_EQUAL(f.readAt(pe, 4), true);
  TEST_EQUAL(pe.identifier, "test")
  TEST_EQUAL(pe.description, " ##0")

  FCVec fv(fev);
  TEST_EQUAL(fv.size(), 4);

  // read, then reset and start reading again
  f.reset();
  TEST_EQUAL(f.cacheChunk(2), true)
  TEST_EQUAL(f.size(), 2)
  TEST_EQUAL(f.activateCache(), true)
  TEST_EQUAL(f.size(), 2)
  FASTAFile::FASTAEntry pe3, pe4;
  TEST_EQUAL(f.readAt(pe3, 0), true);
  pe4 = f.chunkAt(0);
  TEST_TRUE(pe3 == pe4)
  TEST_EQUAL(pe3.description, "This is the description of the first protein")
  pe4 = f.chunkAt(1);
  TEST_EQUAL(pe3 == pe4, false)
  TEST_EQUAL(pe4.description, "This is the description of the second protein")

  f.reset();
  TEST_EQUAL(f.cacheChunk(2), true)
  TEST_EQUAL(f.size(), 2)
  TEST_EQUAL(f.activateCache(), true)
  TEST_EQUAL(f.size(), 2)
  FASTAFile::FASTAEntry pe5, pe6;
  TEST_EQUAL(f.readAt(pe5, 0), true);
  pe6 = f.chunkAt(0);
  TEST_TRUE(pe5 == pe6)
  TEST_EQUAL(pe5.description, "This is the description of the first protein")
  pe6 = f.chunkAt(1);
  TEST_EQUAL(pe5 == pe6, false)
  TEST_EQUAL(pe6.description, "This is the description of the second protein")

END_SECTION

START_SECTION(Result findDecoyString(FASTAContainer<T>& proteins))
// test without decoys in input
  FASTAContainer<TFI_File> f1{OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta")};
  DecoyHelper::Result r1 = {false, "?", true};
  TEST_EQUAL(DecoyHelper::findDecoyString(f1) == r1,true)
  // test with decoys in input
  FASTAContainer<TFI_File> f2{OPENMS_GET_TEST_DATA_PATH("FASTAContainer_test.fasta")};
  DecoyHelper::Result r2 = {true, "DECOY_", true};
  TEST_EQUAL(DecoyHelper::findDecoyString(f2) == r2, true);
END_SECTION

START_SECTION(Result countDecoys(FASTAContainer<T>& proteins))
  // test without decoys in input
  FCFile f1{OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta")};
  std::unordered_map<std::string, std::pair<Size, Size>> decoy_count;
  std::unordered_map<std::string, std::string> decoy_case_sensitive;
  DecoyHelper::DecoyStatistics ds1 = {decoy_count, decoy_case_sensitive,0,0,5};
  TEST_EQUAL(DecoyHelper::countDecoys(f1) == ds1, true)
  // test with decoys in input
  FCFile f2{OPENMS_GET_TEST_DATA_PATH("FASTAContainer_test.fasta")};
  decoy_case_sensitive["decoy_"] = "DECOY_";
  decoy_count["decoy_"] = std::make_pair(3,0);
  DecoyHelper::DecoyStatistics ds2 = { decoy_count, decoy_case_sensitive, 3, 0, 6};
  TEST_EQUAL(DecoyHelper::countDecoys(f2) == ds2, true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



