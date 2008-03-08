// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/StringList.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(StringList, "$Id: String_test.C 2734 2008-02-12 12:01:25Z marc_sturm $")

/////////////////////////////////////////////////////////////

StringList* ptr = 0;
CHECK(StringList())
	ptr = new StringList;
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~StringList())
	delete ptr;
RESULT

CHECK((template<typename StringType> StringList& operator<<(const StringType& string)))
	StringList list;
	list << "a" << "b" << "c" << "a";
	TEST_EQUAL(list.size(),4);
	TEST_STRING_EQUAL(list[0],"a");
	TEST_STRING_EQUAL(list[1],"b");
	TEST_STRING_EQUAL(list[2],"c");
	TEST_STRING_EQUAL(list[3],"a");
RESULT

CHECK(static StringList getYesNoList())
	StringList list = StringList::getYesNoList();
	TEST_EQUAL(list.size(),2);
	TEST_STRING_EQUAL(list[0],"yes");
	TEST_STRING_EQUAL(list[1],"no");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
