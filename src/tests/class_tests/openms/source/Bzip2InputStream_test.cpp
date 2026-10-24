// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/Bzip2InputStream.h>

using namespace OpenMS;


///////////////////////////

START_TEST(Bzip2InputStream, "$Id$")
Bzip2InputStream* ptr = nullptr;
Bzip2InputStream* nullPointer = nullptr;
START_SECTION(Bzip2InputStream(const   char* const     file_name))
	TEST_EXCEPTION(Exception::FileNotFound, Bzip2InputStream bzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	ptr = new Bzip2InputStream(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->getIsOpen(),true)
END_SECTION

START_SECTION((~Bzip2InputStream()))
	delete ptr;
END_SECTION

START_SECTION(Bzip2InputStream(const std::string& file_name))
	TEST_EXCEPTION(Exception::FileNotFound, Bzip2InputStream bzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	std::string filename = OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2");
	ptr = new Bzip2InputStream(filename);
	TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->getIsOpen(),true)
	delete ptr;
END_SECTION

START_SECTION(virtual Size readBytes(unsigned char *const to_fill, const Size max_to_read))
	
	Bzip2InputStream bzip(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	char buffer[31];
	buffer[30] = buffer[29] = '\0';
	unsigned char*  xml_buffer = reinterpret_cast<unsigned char* >(buffer);
		TEST_EQUAL(bzip.getIsOpen(),true)
	TEST_EQUAL(bzip.readBytes(xml_buffer,(Size)10),10)
	TEST_EQUAL(bzip.readBytes(&xml_buffer[10],(Size)10),10)
	TEST_EQUAL(bzip.readBytes(&xml_buffer[20],(Size)9),9)
	TEST_EQUAL(std::string(buffer),std::string("Was decompression successful?"))
	TEST_EQUAL(bzip.getIsOpen(),true)
	TEST_EQUAL(bzip.readBytes(&xml_buffer[30],(Size)10),1)
	TEST_EQUAL(bzip.getIsOpen(),false)

END_SECTION

START_SECTION(UInt64 curPos() const)
	Bzip2InputStream bzip(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
  TEST_EQUAL(bzip.curPos(), 0)
	char buffer[31];
	buffer[30] = buffer[29] = '\0';
	unsigned char*  xml_buffer = reinterpret_cast<unsigned char* >(buffer);
	bzip.readBytes(xml_buffer,(Size)10);
	TEST_EQUAL(bzip.curPos(),10)
	
END_SECTION

START_SECTION(bool getIsOpen() const)
	//test above
	NOT_TESTABLE
END_SECTION

START_SECTION(virtual const char16_t* getContentType() const)
	Bzip2InputStream bzip2(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
  char16_t* xmlch_nullPointer = nullptr;
  TEST_EQUAL(bzip2.getContentType(),xmlch_nullPointer)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION(Zero byte reads preserve stream position)
  Bzip2InputStream stream(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
  TEST_EQUAL(stream.readBytes(nullptr, 0), 0)
  TEST_EQUAL(stream.curPos(), 0)
  unsigned char first[3]{};
  TEST_EQUAL(stream.readBytes(first, 3), 3)
  TEST_EQUAL(stream.curPos(), 3)
END_SECTION
END_TEST
