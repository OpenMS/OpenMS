// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: David Wojnar$
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>

using namespace OpenMS;
using namespace std;


///////////////////////////

START_TEST(Bzip2Ifstream_test, "$Id$")

Bzip2Ifstream* ptr = 0;
START_SECTION((Bzip2Ifstream()))
	ptr = new Bzip2Ifstream;
	TEST_NOT_EQUAL(ptr, 0)
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
