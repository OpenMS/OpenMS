// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(DateTime, "$Id$")

/////////////////////////////////////////////////////////////

DateTime* ptr = nullptr;
DateTime* nullPointer = nullptr;
START_SECTION((DateTime& operator= (const DateTime& source)))
  DateTime date, date2;
  date.set("2006-12-12 11:59:59");
  TEST_EQUAL(date==date2,false);
	date2 = date;
	TEST_EQUAL(date==date2,true);
END_SECTION

START_SECTION((DateTime()))
	ptr = new DateTime();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((DateTime(const DateTime& date)))
	DateTime date1;
	DateTime date2;
	DateTime date3;
	
	date1.set("2006-12-12 11:59:59");
	date2 = DateTime(date1);
	TEST_EQUAL(date1 == date2, true)
END_SECTION

START_SECTION((DateTime(const QDateTime& date)))
	QDateTime date3;
	DateTime date2;
	QTime time;
	QDate date(2006, 12, 12);	
	time.setHMS(11, 59, 59);

	date3.setTime(time);
	date3.setDate(date);
	
	date2.set("2006-12-12 11:59:59");
	DateTime date1(date3);
	
	TEST_EQUAL(date1 == date2, true)
END_SECTION

START_SECTION((void clear()))
	DateTime date1;
	DateTime date2;
	date1.set("2006-12-12 11:59:59");
	date1.clear();
	TEST_EQUAL(date1 == date2, true)
END_SECTION

START_SECTION((String get() const))
	DateTime date_time;
	date_time.set("1999-11-24 14:24:31");
	TEST_EQUAL(date_time.get(),"1999-11-24 14:24:31")
END_SECTION

START_SECTION((void get(UInt& month, UInt& day, UInt& year, UInt& hour, UInt& minute, UInt& second) const))
	DateTime date;
	UInt month;
	UInt day;
	UInt year;
	UInt hour; 
	UInt minute; 
	UInt second;
	
	date.set("2006-12-14 11:59:58");
	date.get(month, day, year, hour, minute, second);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		
END_SECTION

START_SECTION((void getDate(UInt& month, UInt& day, UInt& year) const))
	DateTime date;
	UInt month;
	UInt day;
	UInt year;
	
	date.set("2006-12-14 21:12:02");

	date.getDate(month, day, year);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	

END_SECTION

START_SECTION((String getDate() const))
	DateTime date;
	date.set("2006-12-14 21:12:02");
	TEST_STRING_EQUAL(date.getDate(), String("2006-12-14"))	

END_SECTION

START_SECTION((void getTime(UInt& hour, UInt& minute, UInt& second) const))
	DateTime date;
	UInt hour; 
	UInt minute; 
	UInt second;
	
	date.set("2006-12-14 11:59:58");

	date.getTime(hour, minute, second);
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		

END_SECTION

START_SECTION((String getTime() const))
	DateTime date;
	date.set("2006-12-14 11:59:58");
	TEST_STRING_EQUAL(date.getTime(), "11:59:58")		
END_SECTION

START_SECTION((void set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second)))
	DateTime date;
	UInt month = 12;
	UInt day = 14;
	UInt year = 2006;
	UInt hour = 11; 
	UInt minute = 59; 
	UInt second = 58;
	
	date.set(month, day, year, hour, minute, second);
	date.get(month, day, year, hour, minute, second);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		
END_SECTION

START_SECTION((void set(const String &date)))
	DateTime date_time;
	date_time.set("1999-11-24 14:24:31");
	TEST_EQUAL(date_time.get(), "1999-11-24 14:24:31")

	date_time.set("01.02.2000 14:24:32");
	TEST_EQUAL(date_time.get(), "2000-02-01 14:24:32")

	date_time.set("01/02/2000 14:24:32");
	TEST_EQUAL(date_time.get(), "2000-01-02 14:24:32")

	date_time.set("2005-11-13T10:58:57");
	TEST_EQUAL(date_time.get(), "2005-11-13 10:58:57")

	date_time.set("2008-11-13 10:59:57");
	TEST_EQUAL(date_time.get(), "2008-11-13 10:59:57")

	date_time.set("2006-12-14Z");
	TEST_EQUAL(date_time.get(), "2006-12-14 00:00:00")

	date_time.set("2006-12-14+11:00");
	TEST_EQUAL(date_time.get(), "2006-12-14 11:00:00")

	// test if get is able to ignore the +02:00 timezone part / with and without milliseconds
	// this test is due to #209
	date_time.set("2011-08-05T15:32:07.468+02:00");
	TEST_EQUAL(date_time.get(), "2011-08-05 15:32:07")

	date_time.set("2011-08-05T15:32:07+02:00");
	TEST_EQUAL(date_time.get(), "2011-08-05 15:32:07")

	TEST_EXCEPTION(Exception::ParseError, date_time.set("2006ff-12-14+11:00"))
	TEST_EXCEPTION(Exception::ParseError, date_time.set("2006-12-14-11:00"))
	TEST_EXCEPTION(Exception::ParseError, date_time.set("2006-12-14Z11:00"))
	TEST_EXCEPTION(Exception::ParseError, date_time.set("-2006-12-14Z11:00"))
	
END_SECTION

START_SECTION((void setDate(UInt month, UInt day, UInt year)))
	DateTime date;
	UInt month = 12;
	UInt day = 14;
	UInt year = 2006;

	date.setDate(month, day, year);

	date.getDate(month, day, year);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	

END_SECTION

START_SECTION((void setDate(const String &date)))
	DateTime date;
	UInt month;
	UInt day;
	UInt year;
	
	date.set("2006-12-14 11:59:58");

	date.getDate(month, day, year);
	TEST_EQUAL(month, 12)	
	TEST_EQUAL(day, 14)	
	TEST_EQUAL(year, 2006)	
END_SECTION

START_SECTION((void setTime(UInt hour, UInt minute, UInt second)))
	DateTime date;
	UInt hour; 
	UInt minute; 
	UInt second;
	
	date.setTime(11, 59, 58);

	date.getTime(hour, minute, second);
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		

END_SECTION

START_SECTION((void setTime(const String &date)))
	DateTime date;
	UInt hour; 
	UInt minute; 
	UInt second;
	
	date.setTime("11:59:58");

	date.getTime(hour, minute, second);
	TEST_EQUAL(hour, 11)	
	TEST_EQUAL(minute, 59)	
	TEST_EQUAL(second, 58)		

END_SECTION

START_SECTION(([EXTRA] Three digit year should get leading zero according to Qt 4.4.3 documentation ))
	// This is a regression test.  Leave it here even if the issue gets hacked away in DateTime.
	DateTime one_moment_in_time;
  one_moment_in_time.set(5,4,666,3,2,1);

	// this behaviour is not critical and does not work on Qt 4.3 machines
	// so the leading zero is not checked! (who really needs dates before the year 1000 in this library?) 
  TEST_EQUAL(one_moment_in_time.get().hasSubstring("666-05-04 03:02:01"), true);
	
END_SECTION

START_SECTION((static DateTime now()))
  TEST_EQUAL(DateTime::now().isValid(), true)
END_SECTION

START_SECTION((~DateTime()))
	ptr = new DateTime();
	delete ptr;
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
