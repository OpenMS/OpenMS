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
#include <OpenMS/FORMAT/GzipInputStream.h>
#include <OpenMS/DATASTRUCTURES/String.h>
using namespace OpenMS;


///////////////////////////

START_TEST(GzipInputStream, "$Id$")

GzipInputStream* ptr = 0;
GzipInputStream* nullPointer = 0;
START_SECTION(GzipInputStream(const char *const file_name))
	TEST_EXCEPTION(Exception::FileNotFound, GzipInputStream gzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	ptr = new GzipInputStream(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->getIsOpen(),true)
END_SECTION

START_SECTION((~GzipInputStream()))
	delete ptr;
END_SECTION

START_SECTION(GzipInputStream(const String& file_name))
	TEST_EXCEPTION(Exception::FileNotFound, GzipInputStream gzip2(OPENMS_GET_TEST_DATA_PATH("ThisFileDoesNotExist")))
	String filename(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	ptr = new GzipInputStream(filename);
	TEST_NOT_EQUAL(ptr, nullPointer)
	TEST_EQUAL(ptr->getIsOpen(),true)
	delete ptr;
END_SECTION

START_SECTION(virtual XMLSize_t readBytes(XMLByte *const to_fill, const XMLSize_t max_to_read))
	
	GzipInputStream gzip(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
	char buffer[31];
	buffer[30] = buffer[29] = '\0';
	XMLByte*  xml_buffer = reinterpret_cast<XMLByte* >(buffer);
		TEST_EQUAL(gzip.getIsOpen(),true)
	TEST_EQUAL(gzip.readBytes(xml_buffer,(XMLSize_t)10),10)
	TEST_EQUAL(gzip.readBytes(&xml_buffer[10],(XMLSize_t)10),10)
	TEST_EQUAL(gzip.readBytes(&xml_buffer[20],(XMLSize_t)9),9)
	TEST_EQUAL(String(buffer), String("Was decompression successful?"))
	TEST_EQUAL(gzip.getIsOpen(),true)
	TEST_EQUAL(gzip.readBytes(&xml_buffer[30],(XMLSize_t)10),1)
	TEST_EQUAL(gzip.getIsOpen(),false)

END_SECTION

START_SECTION(XMLFilePos curPos() const)
	GzipInputStream gzip(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
  TEST_EQUAL(gzip.curPos(), 0)
	char buffer[31];
	buffer[30] = buffer[29] = '\0';
	XMLByte*  xml_buffer = reinterpret_cast<XMLByte* >(buffer);
	gzip.readBytes(xml_buffer,(XMLSize_t)10);
	TEST_EQUAL(gzip.curPos(),10)
	
END_SECTION

START_SECTION(bool getIsOpen() const )
	//test above
	NOT_TESTABLE
END_SECTION

START_SECTION(virtual const XMLCh* getContentType() const )
	GzipInputStream gzip2(OPENMS_GET_TEST_DATA_PATH("GzipIfStream_1.gz"));
  XMLCh* xmlch_nullPointer = 0;
  TEST_EQUAL(gzip2.getContentType(),xmlch_nullPointer)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
