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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/IntList.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IntList, "$$")

/////////////////////////////////////////////////////////////

IntList* ptr = 0;
IntList* nullPointer = 0;
START_SECTION(IntList())
	ptr = new IntList;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~IntList())
	delete ptr;
END_SECTION

START_SECTION(static IntList create(const String& list))
	IntList list = IntList::create("1,5");
	TEST_EQUAL(list.size(),2);
	TEST_EQUAL(list[0],1);
	TEST_EQUAL(list[1],5);

	IntList list2 = IntList::create("2");
	TEST_EQUAL(list2.size(),1);
	TEST_EQUAL(list2[0],2);

	IntList list3 = IntList::create("");
	TEST_EQUAL(list3.size(),0);
END_SECTION

START_SECTION(static IntList create(const StringList& list))
	IntList list = IntList::create(StringList::create("1,5"));
	TEST_EQUAL(list.size(),2);
	TEST_EQUAL(list[0],1);
	TEST_EQUAL(list[1],5);
	
	IntList list2 = IntList::create(StringList::create("2"));
	TEST_EQUAL(list2.size(),1);
	TEST_EQUAL(list2[0],2);

	IntList list3 = IntList::create(StringList::create(""));
	TEST_EQUAL(list3.size(),0);
	
	TEST_EXCEPTION(Exception::ConversionError,IntList::create(StringList::create("ein,exception")));
END_SECTION


START_SECTION(IntList(const IntList& rhs))
	IntList list = IntList::create("1,3");
	IntList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);
END_SECTION

START_SECTION(IntList(const std::vector<Int>& rhs))
	std::vector<Int> list;
	list.push_back(1);
	list.push_back(3);
	IntList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);
END_SECTION

START_SECTION(IntList(const std::vector<UInt>& rhs))
	std::vector<UInt> list;
	list.push_back(1);
	list.push_back(2);
	IntList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],2);

END_SECTION

START_SECTION(IntList& operator=(const IntList& rhs))
	IntList list = IntList::create("1,3");
	IntList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);

END_SECTION

START_SECTION(IntList& operator=(const std::vector<Int>& rhs))
	std::vector<Int> list;
	list.push_back(1);
	list.push_back(3);
	IntList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);

END_SECTION

START_SECTION(IntList& operator=(const std::vector<UInt>& rhs))
	std::vector<UInt> list;
	list.push_back(1);
	list.push_back(3);
	IntList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_EQUAL(list2[0],1);
	TEST_EQUAL(list2[1],3);
END_SECTION

START_SECTION((template<typename IntType> IntList& operator<<(IntType value)))
	IntList list;
	list << 1 << 2 << 3 << 1;
	TEST_EQUAL(list.size(),4);
	TEST_EQUAL(list[0],1);
	TEST_EQUAL(list[1],2);
	TEST_EQUAL(list[2],3);
	TEST_EQUAL(list[3],1);
END_SECTION

START_SECTION(bool contains(Int s) const)
	IntList list = IntList::create("1,3");
	TEST_EQUAL(list.contains(1),true)
	TEST_EQUAL(list.contains(3),true)
	TEST_EQUAL(list.contains(4),false)	
	TEST_EQUAL(list.contains(2),false)
	TEST_EQUAL(list.contains(0),false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

