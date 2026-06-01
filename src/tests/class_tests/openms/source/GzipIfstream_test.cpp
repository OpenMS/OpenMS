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
#include <OpenMS/FORMAT/GzipIfstream.h>
using namespace OpenMS;



///////////////////////////

START_TEST(GzipIfstream, "$Id$")

GzipIfstream* ptr = nullptr;
GzipIfstream* nullPointer = nullptr;
START_SECTION((GzipIfstream()))
	ptr = new GzipIfstream;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~GzipIfstream()))
	delete ptr;
END_SECTION

START_SECTION(GzipIfstream(const char * filename))
	TEST_EXCEPTION(Exception::FileNotFound, GzipIfstream gzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	
	GzipIfstream gzip(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	
	TEST_EQUAL(gzip.streamEnd(), false)
	TEST_EQUAL(gzip.isOpen(),true)
	char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EQUAL(29, gzip.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))

END_SECTION

START_SECTION(void open(const char *filename))
	GzipIfstream gzip;
	TEST_EXCEPTION(Exception::FileNotFound, gzip.open(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	
	gzip.open(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	
	TEST_EQUAL(gzip.streamEnd(), false)
	TEST_EQUAL(gzip.isOpen(),true)
	char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	TEST_EQUAL(29, gzip.read(buffer, len))
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
	
END_SECTION

START_SECTION(size_t read(char *s, size_t n))
	//tested in open(const char * filename)
	GzipIfstream gzip(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1_corrupt.gz"));
		char buffer[30];
	buffer[29] = '\0';
	size_t len = 29;
	//gzip.updateCRC32(buffer,10);
	//Why does that throw a "naked" Exception instead of a ConversionError?
	//~ TEST_EXCEPTION(Exception::BaseException,gzip.read(&buffer[9],10));
//	gzip.updateCRC32(&buffer[9],19);
//	TEST_EQUAL(gzip.isCorrupted(),true)
	
	GzipIfstream gzip2(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	TEST_EQUAL(gzip2.isOpen(),true)
	gzip2.read(buffer, len);
	TEST_EQUAL(1, gzip2.read(buffer,10));
	TEST_EQUAL(gzip2.isOpen(), false)
	TEST_EQUAL(gzip2.streamEnd(),true)
	
	gzip2.open(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1_corrupt.gz"));

	 //gzip2.updateCRC32(buffer,(size_t)30);
	//TEST_EQUAL(gzip2.isCorrupted(),true )
	gzip2.close();
	TEST_EQUAL(gzip2.isOpen(), false)
	TEST_EQUAL(gzip2.streamEnd(),true)
	TEST_EXCEPTION(Exception::IllegalArgument, gzip2.read(buffer,10))
	gzip2.close();			
	TEST_EQUAL(gzip2.isOpen(), false)
	TEST_EQUAL(gzip2.streamEnd(),true)
	TEST_EXCEPTION(Exception::IllegalArgument, gzip2.read(buffer,10))
	gzip2.open(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	
		TEST_EQUAL(5, gzip2.read(buffer, 5))
					//			gzip2.updateCRC32(buffer,5);
		TEST_EQUAL(5, gzip2.read(&buffer[5], 5))
						//gzip2.updateCRC32(&buffer[5],5);
		TEST_EQUAL(5, gzip2.read(&buffer[10], 5))
				//		gzip2.updateCRC32(&buffer[10],5);
		TEST_EQUAL(5, gzip2.read(&buffer[15], 5))
					//				gzip2.updateCRC32(&buffer[15],5);
		TEST_EQUAL(5, gzip2.read(&buffer[20], 5))
						//	gzip2.updateCRC32(&buffer[20],5);
					TEST_EQUAL(4, gzip2.read(&buffer[25], 4))
						//			gzip2.updateCRC32(&buffer[25],4);
					char end_of_file[1];
					TEST_EQUAL(1,gzip2.read(end_of_file,2))
		//							gzip2.updateCRC32(end_of_file,1);
					TEST_EQUAL(gzip2.streamEnd(),true)				
					buffer[29]= '\0';				
//	TEST_EQUAL(gzip2.isCorrupted(),false)
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
END_SECTION

START_SECTION(void close())
	//tested in read
	NOT_TESTABLE
END_SECTION
START_SECTION(bool streamEnd() const )
	//!!!tested in open(const char * filename) and read
	NOT_TESTABLE
END_SECTION
START_SECTION(bool isOpen() const)
	//tested in open(const char * filename) and read
	NOT_TESTABLE
END_SECTION
/*
(updateCRC32(char* s, size_t n))
	//tested in open(const char * filename) and read
	_TESTABLE
_SECTION
_SECTION(isCorrupted())
	//tested in open(const char * filename) and read
	TESTABLE
_SECTION*/

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
