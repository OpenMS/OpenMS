// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/SYSTEM/File.h>

///////////////////////////

START_TEST(DTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

CsvFile* ptr = nullptr;
CsvFile* nullPointer = nullptr;
START_SECTION(CsvFile())
	ptr = new CsvFile;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~CsvFile())
	delete ptr;
END_SECTION

#if 0

// Something is terribly wrong here ... looks like an unintentional commit?

START_SECTION(CsvFile(const String& filename, char is = ',',bool ie = false, Int first_n = -1))
//tested in getRow
TEST_EXCEPTION(Exception::FileNotFound, CsvFile("CsvFile_1.csv"))
END_SECTION

START_SECTION(void load(const String& filename, char is = ',', bool ie = false, Int first_n = -1))
//tested in getRow
TEST_EXCEPTION(Exception::FileNotFound, f1.load("CsvFile_1.csv"))


END_SECTION

#endif

START_SECTION(bool getRow(Size row,StringList &list))
	TOLERANCE_ABSOLUTE(0.01)
	CsvFile f1,f3,f4;

	CsvFile f2(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"), '\t');
	StringList list;
	f2.getRow(0,list);
	TEST_EQUAL(list,ListUtils::create<String>("hello,world"))
	f2.getRow(1,list);
	TEST_EQUAL(list,ListUtils::create<String>("the,dude"))
	f2.getRow(2,list);
	TEST_EQUAL(list,ListUtils::create<String>("spectral,search"))

	f3.load(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"),'\t');
	f3.getRow(0,list);
	TEST_EQUAL(list,ListUtils::create<String>("hello,world"))
	f3.getRow(1,list);
	TEST_EQUAL(list,ListUtils::create<String>("the,dude"))
	f3.getRow(2,list);
	TEST_EQUAL(list,ListUtils::create<String>("spectral,search"))

	f4.load(OPENMS_GET_TEST_DATA_PATH("CsvFile_2.csv"),'\t',true);
	f4.getRow(0,list);
	TEST_EQUAL(list,ListUtils::create<String>("hello,world"))
	f4.getRow(1,list);
	TEST_EQUAL(list,ListUtils::create<String>("the,dude"))
	f4.getRow(2,list);
	TEST_EQUAL(list,ListUtils::create<String>("spectral,search"))

END_SECTION

START_SECTION(void store(const String& filename))
	CsvFile f1,f2;
	StringList list;

	f1.load(OPENMS_GET_TEST_DATA_PATH("CsvFile_2.csv"), '\t', true); // load from a file
	String tmpfile = File::getTemporaryFile();
  f1.store(tmpfile);          // store into a new one
	f2.load(tmpfile, '\t', true); // load the new one
	f2.getRow(0,list);
	TEST_EQUAL(list,ListUtils::create<String>("hello,world"))
	f2.getRow(1,list);
	TEST_EQUAL(list,ListUtils::create<String>("the,dude"))
	f2.getRow(2,list);
	TEST_EQUAL(list,ListUtils::create<String>("spectral,search"))
END_SECTION

START_SECTION(void addRow(const StringList& list))
	CsvFile f1, f2;
	StringList list;

	f1.addRow(ListUtils::create<String>("first,second,third"));
	f1.addRow(ListUtils::create<String>("4,5,6"));
  
  String tmpfile = File::getTemporaryFile();
	f1.store(tmpfile);
	f2.load(tmpfile, ',', false);
	f2.getRow(0,list);
	TEST_EQUAL(list, ListUtils::create<String>("first,second,third"))
	f2.getRow(1,list);
	TEST_EQUAL(list, ListUtils::create<String>("4,5,6"))
END_SECTION

START_SECTION(void clear())
	CsvFile f1;
	StringList list;

	f1.addRow(ListUtils::create<String>("hello,world"));
	f1.getRow(0, list);
	TEST_EQUAL(list, ListUtils::create<String>("hello,world"))
	f1.clear();
	TEST_EXCEPTION(Exception::InvalidIterator, f1.getRow(0, list))
END_SECTION

START_SECTION(std::vector<String>::size_type CsvFile::rowCount())
	CsvFile f1(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"), '\t');
	// file has 5 lines, 2 of them comments
	TEST_EQUAL(f1.rowCount(),3)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
