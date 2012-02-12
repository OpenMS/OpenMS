// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/FORMAT/Bzip2InputStream.h>

using namespace OpenMS;


///////////////////////////

START_TEST(Bzip2InputStream, "$Id$")

Bzip2InputStream* ptr = 0;
Bzip2InputStream* nullPointer = 0;
START_SECTION(Bzip2InputStream(const   char* const     file_name))
	TEST_EXCEPTION(Exception::FileNotFound, Bzip2InputStream bzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	ptr = new Bzip2InputStream(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->getIsOpen(),true)
END_SECTION

START_SECTION((~Bzip2InputStream()))
	delete ptr;
END_SECTION

START_SECTION(Bzip2InputStream(const String& file_name))
	TEST_EXCEPTION(Exception::FileNotFound, Bzip2InputStream bzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	String filename = OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2");
	ptr = new Bzip2InputStream(filename);
	TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->getIsOpen(),true)
	delete ptr;
END_SECTION

START_SECTION(virtual XMLSize_t readBytes(XMLByte *const to_fill, const XMLSize_t max_to_read))
	
	Bzip2InputStream bzip(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
	char buffer[31];
	buffer[30] = buffer[29] = '\0';
	XMLByte*  xml_buffer = reinterpret_cast<XMLByte* >(buffer);
		TEST_EQUAL(bzip.getIsOpen(),true)
	TEST_EQUAL(bzip.readBytes(xml_buffer,(XMLSize_t)10),10)
	TEST_EQUAL(bzip.readBytes(&xml_buffer[10],(XMLSize_t)10),10)
	TEST_EQUAL(bzip.readBytes(&xml_buffer[20],(XMLSize_t)9),9)
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
	TEST_EQUAL(bzip.getIsOpen(),true)
	TEST_EQUAL(bzip.readBytes(&xml_buffer[30],(XMLSize_t)10),1)
	TEST_EQUAL(bzip.getIsOpen(),false)

END_SECTION

START_SECTION(XMLFilePos curPos() const)
	Bzip2InputStream bzip(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
  TEST_EQUAL(bzip.curPos(), 0)
	char buffer[31];
	buffer[30] = buffer[29] = '\0';
	XMLByte*  xml_buffer = reinterpret_cast<XMLByte* >(buffer);
	bzip.readBytes(xml_buffer,(XMLSize_t)10);
	TEST_EQUAL(bzip.curPos(),10)
	
END_SECTION

START_SECTION(bool getIsOpen() const)
	//test above
	NOT_TESTABLE
END_SECTION

START_SECTION(virtual const XMLCh* getContentType() const)
	Bzip2InputStream bzip2(OPENMS_GET_TEST_DATA_PATH("Bzip2IfStream_1.bz2"));
  XMLCh* xmlch_nullPointer = 0;
  TEST_EQUAL(bzip2.getContentType(),xmlch_nullPointer)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
