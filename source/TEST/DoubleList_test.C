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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DoubleList.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(DoubleList, "$Id$")

/////////////////////////////////////////////////////////////

DoubleList* ptr = 0;
DoubleList* nullPointer = 0;
START_SECTION(DoubleList())
	ptr = new DoubleList;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~DoubleList())
	delete ptr;
END_SECTION

START_SECTION(static DoubleList create(const String& list))
	DoubleList list = DoubleList::create("1.222,5.33789");
	TEST_EQUAL(list.size(),2);
	TEST_REAL_SIMILAR(list[0],1.222);
	TEST_REAL_SIMILAR(list[1],5.33789);

	DoubleList list2 = DoubleList::create("2.33334");
	TEST_EQUAL(list2.size(),1);
	TEST_REAL_SIMILAR(list2[0],2.33334);

	DoubleList list3 = DoubleList::create("");
	TEST_EQUAL(list3.size(),0);
END_SECTION

START_SECTION(static DoubleList create(const StringList& list))
	DoubleList list = DoubleList::create(StringList::create("1.222,5.33789"));
	TEST_EQUAL(list.size(),2);
	TEST_REAL_SIMILAR(list[0],1.222);
	TEST_REAL_SIMILAR(list[1],5.33789);

	DoubleList list2 = DoubleList::create(StringList::create("2.33334"));
	TEST_EQUAL(list2.size(),1);
	TEST_REAL_SIMILAR(list2[0],2.33334);

	DoubleList list3 = DoubleList::create(StringList::create(""));
	TEST_EQUAL(list3.size(),0);
	TEST_EXCEPTION(Exception::ConversionError,DoubleList::create(StringList::create("ein,exception")));
END_SECTION

START_SECTION(DoubleList(const DoubleList& rhs))
	DoubleList list = DoubleList::create("1.2,3.4");
	DoubleList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.2);
	TEST_REAL_SIMILAR(list2[1],3.4);
END_SECTION

START_SECTION(DoubleList(const std::vector<DoubleReal>& rhs))
	std::vector<DoubleReal> list;
	list.push_back(1.2345);
	list.push_back(3.45678);
	DoubleList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.2345);
	TEST_REAL_SIMILAR(list2[1],3.45678);
END_SECTION

START_SECTION(DoubleList(const std::vector<Real>& rhs))
	std::vector<Real> list;
	list.push_back(1.234f);
	list.push_back(2.345f);
	DoubleList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.234);
	TEST_REAL_SIMILAR(list2[1],2.345);

END_SECTION

START_SECTION(DoubleList& operator=(const DoubleList& rhs))
	DoubleList list = DoubleList::create("1.22,3.33");
	DoubleList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.22);
	TEST_REAL_SIMILAR(list2[1],3.33);

END_SECTION

START_SECTION(DoubleList& operator=(const std::vector<DoubleReal>& rhs))
	std::vector<DoubleReal> list;
	list.push_back(1.22);
	list.push_back(3.67);
	DoubleList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.22);
	TEST_REAL_SIMILAR(list2[1],3.67);
END_SECTION

START_SECTION(DoubleList& operator=(const std::vector<Real>& rhs))
	std::vector<Real> list;
	list.push_back(1.22f);
	list.push_back(3.67f);
	DoubleList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_REAL_SIMILAR(list2[0],1.22);
	TEST_REAL_SIMILAR(list2[1],3.67);
END_SECTION

START_SECTION((template<typename DoubleType> DoubleList& operator<<(DoubleType value)))
	DoubleList list;
	list << 1.2 << 2.3456 << 3.5678999 << 1.2;
	TEST_EQUAL(list.size(),4);
	TEST_EQUAL(list[0],1.2);
	TEST_EQUAL(list[1],2.3456);
	TEST_EQUAL(list[2],3.5678999);
	TEST_EQUAL(list[3],1.2);
END_SECTION

START_SECTION(bool contains(DoubleReal s, DoubleReal tolerance=0.00001) const)
	DoubleList list = DoubleList::create("1.2,3.4");
	TEST_EQUAL(list.contains(1.2),true)
	TEST_EQUAL(list.contains(1.21),false)
	TEST_EQUAL(list.contains(1.19),false)
	TEST_EQUAL(list.contains(1.21,0.02),true)
	TEST_EQUAL(list.contains(1.19,0.02),true)
	TEST_EQUAL(list.contains(3.4),true)
	TEST_EQUAL(list.contains(4.2),false)
	TEST_EQUAL(list.contains(2),false)
	TEST_EQUAL(list.contains(0),false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

