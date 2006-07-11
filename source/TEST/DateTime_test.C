// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

START_TEST(DateTime, "$Id: DateTime_test.C,v 1.1 2006/06/09 23:47:35 nicopfeifer Exp $")

/////////////////////////////////////////////////////////////

DateTime* s_ptr = 0;
CHECK(DateTime& operator= (const DateTime& source))
  // ???
RESULT

CHECK(DateTime())
  // ???
RESULT

CHECK(DateTime(const DateTime& date))
  // ???
RESULT

CHECK(bool operator != (const DateTime& rhs) const)
  // ???
RESULT

CHECK(bool operator == (const DateTime& rhs) const)
  // ???
RESULT

CHECK(void clear())
  // ???
RESULT

CHECK(void get(String& date) const)
	DateTime date_time;
	String date_time_string("1999-11-24 14:24:31");
	String output;
	
	date_time.set(date_time_string);
	
	date_time.get(output);
	TEST_EQUAL(date_time_string, output)

RESULT

CHECK((void get(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year, UnsignedInt& hour, UnsignedInt& minute, UnsignedInt& second) const))
  // ???
RESULT

CHECK((void getDate(UnsignedInt& month, UnsignedInt& day, UnsignedInt& year) const))
  // ???
RESULT

CHECK((void getTime(UnsignedInt& hour, UnsignedInt& minute, UnsignedInt& second) const))
  // ???
RESULT

CHECK((void set(UnsignedInt month, UnsignedInt day, UnsignedInt year, UnsignedInt hour, UnsignedInt minute, UnsignedInt second) throw(Exception::ParseError)))
  // ???
RESULT

CHECK(void set(const String& date) throw(Exception::ParseError))
	DateTime date_time;
	String date_time_string("1999-11-24 14:24:31");
	String output;
	
	date_time.set(date_time_string);
	
	date_time.get(output);
	TEST_EQUAL(date_time_string, output)
RESULT

CHECK((void setDate(UnsignedInt month, UnsignedInt day, UnsignedInt year) throw(Exception::ParseError)))
  // ???
RESULT

CHECK(void setDate(const String& date) throw(Exception::ParseError))
  // ???
RESULT

CHECK((void setTime(UnsignedInt hour, UnsignedInt minute, UnsignedInt second) throw(Exception::ParseError)))
  // ???
RESULT

CHECK(void setTime(const String& date) throw(Exception::ParseError))
  // ???
RESULT

CHECK(void today())
  // ???
RESULT

CHECK(~DateTime())
	s_ptr = new DateTime();
	delete s_ptr;
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
