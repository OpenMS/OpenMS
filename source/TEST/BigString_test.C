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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <iostream>
///////////////////////////
#include <OpenMS/DATASTRUCTURES/BigString.h>
//////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BigString, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
typedef std::pair <String, String> FASTAEntry;

BigString* ptr = 0;

CHECK(BigString())
	ptr = new BigString();
	TEST_EQUAL (ptr->getBigString(),"$");
	TEST_EQUAL (ptr->size(),1);
	TEST_EQUAL (ptr->length(),1);
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~BigString())
	delete ptr;
RESULT

CHECK(void add(FASTAEntry const &new_entry))
	ptr = new BigString();
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	ptr->add(fe);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$");
	const FASTAEntry fe2 ("ENTRY 2","BBBBB");
	ptr->add(fe2);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$BBBBB$");
	TEST_EQUAL (ptr->size(),3);
	TEST_EQUAL (ptr->length(),13);
RESULT

CHECK(void setSeparator(const char sep))
	ptr = new BigString();
	ptr->setSeparator('&');
	TEST_EQUAL(ptr->getSeparator(),'&');
RESULT

CHECK(char getSeparator())
	ptr = new BigString();
	TEST_EQUAL(ptr->getSeparator(),'$');
	ptr->setSeparator('&');
	TEST_EQUAL(ptr->getSeparator(),'&');
RESULT

CHECK(BigString(const BigString &bs))
	ptr = new BigString();
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	ptr->add(fe);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$");
	const FASTAEntry fe2 ("ENTRY 2","BBBBB");
	ptr->add(fe2);
	ptr->setSeparator('&');
	BigString * new_big_string = new BigString (*ptr);
	TEST_EQUAL (ptr->getSeparator(),new_big_string->getSeparator());
	TEST_EQUAL (ptr->getBigString(),new_big_string->getBigString());
	TEST_EQUAL (ptr->size(),new_big_string->size());
	TEST_EQUAL (ptr->length(),new_big_string->length());	
	pair<String, String> result, ptr_result;
	new_big_string->getPeptide(result, 2, 2);
	ptr->getPeptide(ptr_result, 2, 2);
	TEST_EQUAL (ptr_result.first, result.first);
	TEST_EQUAL (ptr_result.second, result.second);
RESULT

CHECK(const String& getBigString() const )
	ptr = new BigString();
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	ptr->add(fe);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$");
	const FASTAEntry fe2 ("ENTRY 2","BBBBB");
	ptr->add(fe2);
	TEST_EQUAL (ptr->getBigString(),"$AAAAA$BBBBB$");
RESULT

CHECK(unsigned int size())
	ptr = new BigString();
	TEST_EQUAL (ptr->size(),1);
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	for (UInt i= 1; i < 10; i++)
	{
		ptr->add(fe);
		TEST_EQUAL (ptr->size(),i+1);
	}
RESULT

CHECK(unsigned int length())
	ptr = new BigString();
	TEST_EQUAL (ptr->length(),1);
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	for (UInt i= 1; i < 10; i++)
	{
		ptr->add(fe);
		TEST_EQUAL (ptr->length(), i * 6 + 1);
	}
RESULT

CHECK(void getPeptide(FASTAEntry& entry, unsigned int start, unsigned int length))
	ptr = new BigString();
	const FASTAEntry fe ("ENTRY 1","AAAAA");
	ptr->add(fe);
	const FASTAEntry fe2 ("ENTRY 2","BBBBB");
	ptr->add(fe2);
	const FASTAEntry fe3 ("ENTRY 3","CCCCC");
	ptr->add(fe3);
	const FASTAEntry fe4 ("ENTRY 4","DDDDD");
	ptr->add(fe4);
	const FASTAEntry fe5 ("ENTRY 5","EEEEE");
	ptr->add(fe5);
	FASTAEntry res;
	ptr->getPeptide(res, 1, 3);
	TEST_EQUAL(res.first,"ENTRY 1");
	TEST_EQUAL(res.second,"AAA");
	ptr->getPeptide(res, 1, 5);
	TEST_EQUAL(res.first,"ENTRY 1");
	TEST_EQUAL(res.second,"AAAAA");
	ptr->getPeptide(res, 3, 2);
	TEST_EQUAL(res.first,"ENTRY 1");
	TEST_EQUAL(res.second,"AA");
	ptr->getPeptide(res, 7, 2);
	TEST_EQUAL(res.first,"ENTRY 2");
	TEST_EQUAL(res.second,"BB");
	ptr->getPeptide(res, 19, 2);
	TEST_EQUAL(res.first,"ENTRY 4");
	TEST_EQUAL(res.second,"DD");
	TEST_EXCEPTION(Exception::InvalidValue, ptr->getPeptide(res, 1,10));
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
