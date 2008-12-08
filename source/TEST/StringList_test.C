// -*- mode: C++; tab-width: 2; -*-
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

START_TEST(StringList, "$Id$")

/////////////////////////////////////////////////////////////

StringList* ptr = 0;
START_SECTION(StringList())
	ptr = new StringList;
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~StringList())
	delete ptr;
END_SECTION

START_SECTION(static StringList create(const String& list))
	StringList list = StringList::create("yes,no");
	TEST_EQUAL(list.size(),2);
	TEST_STRING_EQUAL(list[0],"yes");
	TEST_STRING_EQUAL(list[1],"no");

	StringList list2 = StringList::create("no");
	TEST_EQUAL(list2.size(),1);
	TEST_STRING_EQUAL(list2[0],"no");

	StringList list3 = StringList::create("");
	TEST_EQUAL(list3.size(),0);
END_SECTION

START_SECTION(StringList(const StringList& rhs))
	StringList list = StringList::create("yes,no");
	StringList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_STRING_EQUAL(list2[0],"yes");
	TEST_STRING_EQUAL(list2[1],"no");
END_SECTION

START_SECTION(StringList(const std::vector<String>& rhs))
	std::vector<String> list;
	list.push_back("yes");
	list.push_back("no");
	StringList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_STRING_EQUAL(list2[0],"yes");
	TEST_STRING_EQUAL(list2[1],"no");
END_SECTION

START_SECTION(StringList(const std::vector<std::string>& rhs))
	std::vector<string> list;
	list.push_back("yes");
	list.push_back("no");
	StringList list2(list);
	TEST_EQUAL(list2.size(),2);
	TEST_STRING_EQUAL(list2[0],"yes");
	TEST_STRING_EQUAL(list2[1],"no");

END_SECTION

START_SECTION(StringList& operator=(const StringList& rhs))
	StringList list = StringList::create("yes,no");
	StringList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_STRING_EQUAL(list2[0],"yes");
	TEST_STRING_EQUAL(list2[1],"no");

END_SECTION

START_SECTION(StringList& operator=(const std::vector<String>& rhs))
	std::vector<String> list;
	list.push_back("yes");
	list.push_back("no");
	StringList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_STRING_EQUAL(list2[0],"yes");
	TEST_STRING_EQUAL(list2[1],"no");

END_SECTION

START_SECTION(StringList& operator=(const std::vector<std::string>& rhs))
	std::vector<string> list;
	list.push_back("yes");
	list.push_back("no");
	StringList list2;
	list2 = list;
	TEST_EQUAL(list2.size(),2);
	TEST_STRING_EQUAL(list2[0],"yes");
	TEST_STRING_EQUAL(list2[1],"no");

END_SECTION

START_SECTION((template<typename StringType> StringList& operator<<(const StringType& string)))
	StringList list;
	list << "a" << "b" << "c" << "a";
	TEST_EQUAL(list.size(),4);
	TEST_STRING_EQUAL(list[0],"a");
	TEST_STRING_EQUAL(list[1],"b");
	TEST_STRING_EQUAL(list[2],"c");
	TEST_STRING_EQUAL(list[3],"a");
END_SECTION

START_SECTION(bool contains(const String& s) const)
	StringList list = StringList::create("yes,no");
	TEST_EQUAL(list.contains("yes"),true)
	TEST_EQUAL(list.contains("no"),true)
	TEST_EQUAL(list.contains("jup"),false)	
	TEST_EQUAL(list.contains(""),false)
	TEST_EQUAL(list.contains("noe"),false)
END_SECTION

START_SECTION(void toUpper())
	StringList list = StringList::create("yes,no");
	list.toUpper();
	TEST_EQUAL(list[0],"YES")
	TEST_EQUAL(list[1],"NO")
END_SECTION

START_SECTION(void toLower())
	StringList list = StringList::create("yES,nO");
	list.toLower();
	TEST_EQUAL(list[0],"yes")
	TEST_EQUAL(list[1],"no")
END_SECTION

START_SECTION(String concatenate(const String& glue=""))
	StringList list = StringList::create("1,2,3,4,5");
	TEST_STRING_EQUAL(list.concatenate("g"),"1g2g3g4g5");
	TEST_STRING_EQUAL(list.concatenate(""),"12345");
	
	list.clear();
	TEST_STRING_EQUAL(list.concatenate("g"),"");
	TEST_STRING_EQUAL(list.concatenate(""),"");
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
