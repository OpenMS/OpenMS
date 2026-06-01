// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>

using namespace OpenMS;
using namespace std;


///////////////////////////

START_TEST(Bzip2Ifstream_test, "$Id$")

Bzip2Ifstream* ptr = nullptr;
Bzip2Ifstream* nullPointer = nullptr;
START_SECTION((Bzip2Ifstream()))
	ptr = new Bzip2Ifstream;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~Bzip2Ifstream()))
	delete ptr;
END_SECTION

START_SECTION(Bzip2Ifstream(const char * filename))
	TEST_EXCEPTION(Exception::FileNotFound, Bzip2Ifstream bzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	
	Bzip2Ifstream bzip(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	
	TEST_EQUAL(bzip.streamEnd(), false)
	TEST_EQUAL(bzip.isOpen(),true)
	char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EQUAL(29, bzip.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))

END_SECTION

START_SECTION(void open(const char *filename))
	Bzip2Ifstream bzip;
	TEST_EXCEPTION(Exception::FileNotFound, bzip.open(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	
	bzip.open(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	
	TEST_EQUAL(bzip.streamEnd(), false)
	TEST_EQUAL(bzip.isOpen(),true)
	char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EQUAL(29, bzip.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
	
END_SECTION

START_SECTION(size_t read(char *s, size_t n))
	//tested in open(const char * filename)
	Bzip2Ifstream bzip(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1_corrupt.bz2"));
		char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EXCEPTION(Exception::ParseError, bzip.read(buffer,10))
	
	Bzip2Ifstream bzip2(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	bzip2.read(buffer, len);
	TEST_EQUAL(1, bzip2.read(buffer,10));
	TEST_EQUAL(bzip2.isOpen(), false)
	TEST_EQUAL(bzip2.streamEnd(),true)
	
	bzip2.open(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1_corrupt.bz2"));
	TEST_EXCEPTION(Exception::ParseError, bzip2.read(buffer,10))
	bzip2.close();
	TEST_EQUAL(bzip2.isOpen(), false)
	TEST_EQUAL(bzip2.streamEnd(),true)
	TEST_EXCEPTION(Exception::IllegalArgument, bzip2.read(buffer,10))
	bzip2.close();
	TEST_EQUAL(bzip2.isOpen(), false)
	TEST_EQUAL(bzip2.streamEnd(),true)
	TEST_EXCEPTION(Exception::IllegalArgument, bzip2.read(buffer,10))
	bzip2.open(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	TEST_EQUAL(29, bzip2.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
END_SECTION

START_SECTION(void close())
	//tested in read
	NOT_TESTABLE
END_SECTION
START_SECTION(bool streamEnd() const )
	//tested in open(const char * filename) and read
	NOT_TESTABLE
END_SECTION
START_SECTION(bool isOpen() const )
	//tested in open(const char * filename) and read
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
