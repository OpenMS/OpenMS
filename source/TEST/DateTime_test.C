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
// $Maintainer: Nico Pfeifer $
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(DateTime, "$Id$")

/////////////////////////////////////////////////////////////

DateTime* ptr = 0;
CHECK((DateTime& operator= (const DateTime& source)))
  // ???
RESULT

CHECK((DateTime()))
	ptr = new DateTime();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((DateTime(const DateTime& date)))
	DateTime date1;
	DateTime date2;
	date1.set("2006-12-12 11:59:59");
	date2 = date1;
	TEST_EQUAL(date1 == date2, true)
RESULT

CHECK((bool operator != (const DateTime& rhs) const))
	DateTime date1;
	DateTime date2;
	date1.set("2006-12-12 11:59:59");
	TEST_EQUAL(date1 != date2, true)
	
RESULT

CHECK((bool operator == (const DateTime& rhs) const))
	DateTime date1;
	DateTime date2;
	date1.set("2006-12-12 11:59:59");
	date2 = date1;
	TEST_EQUAL(date1 == date2, true)
RESULT

CHECK((void clear()))
	DateTime date1;
	DateTime date2;
	date1.set("2006-12-12 11:59:59");
	date1.clear();
	TEST_EQUAL(date1 == date2, true)
RESULT

CHECK((void get(String& date) const))
	DateTime date_time;
	String date_time_string("1999-11-24 14:24:31");
	String output;
	
	date_time.set(date_time_string);
	
	date_time.get(output);
	TEST_EQUAL(date_time_string, output)

RESULT

CHECK((void get(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year, UnsignedInt& hour, UnsignedInt& minute, UnsignedInt& second) const))
	DateTime date;
	UnsignedInt month;
	UnsignedInt day;
	UnsignedInt year;
	UnsignedInt hour; 
	UnsignedInt minute; 
	UnsignedInt second;
	
	date.set("2006-12-14 11:59:58");
	date.get(month, day, year, hour, minute, second);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		
RESULT

CHECK((void getDate(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year) const))
	DateTime date;
	UnsignedInt month;
	UnsignedInt day;
	UnsignedInt year;
	
	date.set("2006-12-14");

	date.getDate(month, day, year);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	

RESULT

CHECK((void getDate(String& date) const))
	DateTime date;
	String temp;
	
	date.set("2006-12-14");
	date.getDate(temp);
	TEST_EQUAL(temp, String("2006-12-14"))	

RESULT

CHECK((void getTime(UnsignedInt& hour, UnsignedInt& minute, UnsignedInt& second) const))
	DateTime date;
	UnsignedInt hour; 
	UnsignedInt minute; 
	UnsignedInt second;
	
	date.set("2006-12-14 11:59:58");

	date.getTime(hour, minute, second);
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		

RESULT

CHECK((void getTime(String& time) const))
	DateTime date;
	String tmp;
		
	date.set("2006-12-14 11:59:58");

	date.getTime(tmp);
	TEST_EQUAL(tmp, "11:59:58")		

RESULT

CHECK((void set(UnsignedInt month, UnsignedInt day, UnsignedInt year, UnsignedInt hour, UnsignedInt minute, UnsignedInt second) throw(Exception::ParseError)))
	DateTime date;
	UnsignedInt month = 12;
	UnsignedInt day = 14;
	UnsignedInt year = 2006;
	UnsignedInt hour = 11; 
	UnsignedInt minute = 59; 
	UnsignedInt second = 58;
	
	date.set(month, day, year, hour, minute, second);
	date.get(month, day, year, hour, minute, second);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		
RESULT

CHECK((void set(const String& date) throw(Exception::ParseError)))
	DateTime date_time;
	String date_time_string("1999-11-24 14:24:31");
	String output;
	
	date_time.set(date_time_string);
	
	date_time.get(output);
	TEST_EQUAL(date_time_string, output)
RESULT

CHECK((void setDate(UnsignedInt month, UnsignedInt day, UnsignedInt year) throw(Exception::ParseError)))
	DateTime date;
	UnsignedInt month = 12;
	UnsignedInt day = 14;
	UnsignedInt year = 2006;

	date.setDate(month, day, year);

	date.getDate(month, day, year);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	

RESULT

CHECK((void setDate(const String& date) throw(Exception::ParseError)))
	DateTime date;
	UnsignedInt month;
	UnsignedInt day;
	UnsignedInt year;
	
	date.set("2006-12-14 11:59:58");

	date.getDate(month, day, year);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	

RESULT

CHECK((void setTime(UnsignedInt hour, UnsignedInt minute, UnsignedInt second) throw(Exception::ParseError)))
	DateTime date;
	UnsignedInt hour; 
	UnsignedInt minute; 
	UnsignedInt second;
	
	date.setTime(11, 59, 58);

	date.getTime(hour, minute, second);
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		

RESULT

CHECK((void setTime(const String& date) throw(Exception::ParseError)))
	DateTime date;
	UnsignedInt hour; 
	UnsignedInt minute; 
	UnsignedInt second;
	
	date.setTime("11:59:58");

	date.getTime(hour, minute, second);
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		

RESULT

CHECK((void now()))
	DateTime date1;
	DateTime date2;
	
	date1.now();
	TEST_EQUAL(date1 != date2, true)

RESULT

CHECK((~DateTime()))
	ptr = new DateTime();
	delete ptr;
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
