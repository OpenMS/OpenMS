// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

///////////////////////////

START_TEST(DTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

CsvFile* ptr = 0;
START_SECTION(CsvFile())
	ptr = new CsvFile;
	TEST_NOT_EQUAL(ptr, 0)
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

START_SECTION(void fload(const String& filename, char is = ',', bool ie = false, Int first_n = -1))
//tested in getRow
TEST_EXCEPTION(Exception::FileNotFound, f1.fload("CsvFile_1.csv"))	


END_SECTION

#endif

START_SECTION(bool getRow(Size row,StringList &list))
	TOLERANCE_ABSOLUTE(0.01)
	CsvFile f1,f3,f4;
	
	
	CsvFile f2(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"), '\t');
	StringList list;
	f2.getRow(0,list);
	TEST_EQUAL(list,StringList::create("hello,world"))
	f2.getRow(1,list);
	TEST_EQUAL(list,StringList::create("the,dude"))
	f2.getRow(2,list);
	TEST_EQUAL(list,StringList::create("spectral,search"))
	
	f3.fload(OPENMS_GET_TEST_DATA_PATH("CsvFile_1.csv"),'\t');
	f3.getRow(0,list);
	TEST_EQUAL(list,StringList::create("hello,world"))
	f3.getRow(1,list);
	TEST_EQUAL(list,StringList::create("the,dude"))
	f3.getRow(2,list);
	TEST_EQUAL(list,StringList::create("spectral,search"))
	
	f4.fload(OPENMS_GET_TEST_DATA_PATH("CsvFile_2.csv"),'\t',true);
	f4.getRow(0,list);
	TEST_EQUAL(list,StringList::create("hello,world"))
	f4.getRow(1,list);
	TEST_EQUAL(list,StringList::create("the,dude"))
	f4.getRow(2,list);
	TEST_EQUAL(list,StringList::create("spectral,search"))	
	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
