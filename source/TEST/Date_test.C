// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/Date.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Date, "$Id$")

/////////////////////////////////////////////////////////////

Date* s_ptr = 0;
CHECK((Date()))
	s_ptr = new Date();
	TEST_NOT_EQUAL(s_ptr, 0)
RESULT

CHECK((~Date()))
	delete s_ptr;
RESULT

CHECK((bool isLeapYear(UnsignedInt year) const))
  Date d;
  TEST_EQUAL(d.isLeapYear(1999),false);
  TEST_EQUAL(d.isLeapYear(2000),true);
  TEST_EQUAL(d.isLeapYear(2001),false);
  TEST_EQUAL(d.isLeapYear(1800),false);
  TEST_EQUAL(d.isLeapYear(1900),false);
  TEST_EQUAL(d.isLeapYear(2100),false);
  TEST_EQUAL(d.isLeapYear(2400),true);
RESULT

CHECK((void get(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year) const))
  Date date;
  UnsignedInt d,m,y;
  date.get(m,d,y);
  TEST_EQUAL(m,0);
  TEST_EQUAL(d,0);
  TEST_EQUAL(y,0);
RESULT

CHECK((void set(UnsignedInt month, UnsignedInt day, UnsignedInt year) throw(Exception::ParseError)))
  Date date;
  UnsignedInt d,m,y;
  date.set(12,1,1977);
  date.get(m,d,y);
  TEST_EQUAL(m,12);
  TEST_EQUAL(d,1);
  TEST_EQUAL(y,1977);
  
  //exceptions
  TEST_EXCEPTION(Exception::ParseError,date.set(0,12,1977));
  TEST_EXCEPTION(Exception::ParseError,date.set(12,0,1977));
  TEST_EXCEPTION(Exception::ParseError,date.set(1,32,1977));
  TEST_EXCEPTION(Exception::ParseError,date.set(13,1,1977));
	TEST_EXCEPTION(Exception::ParseError,date.set(02,29,2100));
RESULT

CHECK((bool operator == (const Date& rhs) const))
  Date date, date2;
  TEST_EQUAL(date==date2,true);
  date.set(12,1,1977);
  TEST_EQUAL(date==date2,false);
  date2.set(12,1,1977);
  TEST_EQUAL(date==date2,true);
RESULT

CHECK((bool operator != (const Date& rhs) const))
  Date date, date2;
  TEST_EQUAL(date!=date2,false);
  date.set(12,1,1977);
  TEST_EQUAL(date!=date2,true);
  date2.set(12,1,1977);
  TEST_EQUAL(date!=date2,false);
RESULT

CHECK((Date& operator= (const Date& source)))
  Date date, date2;
  date.set(12,1,1977);
  TEST_EQUAL(date==date2,false);
	date2 = date;
	TEST_EQUAL(date==date2,true);
RESULT

CHECK((Date(const Date& date)))
  Date date;
  date.set(12,1,1977);
	Date date2(date);
	TEST_EQUAL(date==date2,true);
RESULT

CHECK((void set(const String& date) throw(Exception::ParseError)))
  Date date;
  //german
  date.set("01.12.1977");
  UnsignedInt d,m,y;
  date.get(m,d,y);
  TEST_EQUAL(m,12);
  TEST_EQUAL(d,1);
  TEST_EQUAL(y,1977);  

  //english
  date.set("12/01/1977");
  date.get(m,d,y);
  TEST_EQUAL(m,12);
  TEST_EQUAL(d,1);
  TEST_EQUAL(y,1977);

  //iso/ansi
  date.set("1967-12-23");
  date.get(m,d,y);
  TEST_EQUAL(d,23);
  TEST_EQUAL(m,12);
  TEST_EQUAL(y,1967);
    
   //german short
  date.set("6.1.888");
  date.get(m,d,y);
  TEST_EQUAL(m,1);
  TEST_EQUAL(d,6);
  TEST_EQUAL(y,888);

	//exceptions
  TEST_EXCEPTION(Exception::ParseError,date.set("bla"));
  TEST_EXCEPTION(Exception::ParseError,date.set("01.01.01.2005"));
  TEST_EXCEPTION(Exception::ParseError,date.set("f1.01.1977"));
  TEST_EXCEPTION(Exception::ParseError,date.set("01.1x.1977"));
  TEST_EXCEPTION(Exception::ParseError,date.set("01.12.i135"));
  TEST_EXCEPTION(Exception::ParseError,date.set("1135-64-3"));
RESULT

CHECK((void get(String& date) const))
  Date d;
  String s;
  d.get(s);
  TEST_EQUAL(s,"0000-00-00");
  d.set("11.12.1977");
  d.get(s);
  TEST_EQUAL(s,"1977-12-11");
	d.set("02.01.888");
  d.get(s);
  TEST_EQUAL(s,"0888-01-02");
RESULT

CHECK((void clear()))
  Date d;
  String s;
  d.set("11.12.1977");
  d.get(s);
  TEST_EQUAL(s,"1977-12-11");
	d.clear();
  d.get(s);
  TEST_EQUAL(s,"0000-00-00");
RESULT

CHECK(void today())
	// not testable
RESULT

CHECK(static std::string now())
	// not testable
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
