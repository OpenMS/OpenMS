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

CHECK(([EXTRA]~Date()))
	delete s_ptr;
RESULT

CHECK(Date(const QDate &date))
	QDate qd(1999,12,24);
	Date d(qd);
	TEST_EQUAL(d.year(),1999)
	TEST_EQUAL(d.month(),12)
	TEST_EQUAL(d.day(),24)
RESULT

CHECK((void get(UInt& month, UInt& day, UInt& year) const))
  Date date;
  UInt d,m,y;
  date.set("2007-12-03");
  date.get(m,d,y);
  TEST_EQUAL(m,12);
  TEST_EQUAL(d,3);
  TEST_EQUAL(y,2007);
RESULT

CHECK((void set(UInt month, UInt day, UInt year) ))
  Date date;
  UInt d,m,y;
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

CHECK((void set(const String& date) ))
  Date date;
  //german
  date.set("01.12.1977");
  UInt d,m,y;
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
  date.set("06.01.1688");
  date.get(m,d,y);
  TEST_EQUAL(m,1);
  TEST_EQUAL(d,6);
  TEST_EQUAL(y,1688);

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
	d.set("02.01.1999");
  d.get(s);
  TEST_EQUAL(s,"1999-01-02");
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

CHECK((void today()))
NOT_TESTABLE
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
