// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

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

START_SECTION((DateTime()))
{
  ptr = new DateTime();
  TEST_NOT_EQUAL(ptr, nullPointer)
  delete ptr;
}
END_SECTION

START_SECTION((~DateTime()))
{
  ptr = new DateTime();
  delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
// Copy constructor, move constructor, assignment operator, move assignment operator, equality

START_SECTION((DateTime(const DateTime& date)))
{
  DateTime date1;
  DateTime date2;
  DateTime date3;

  date1.set("2006-12-12 11:59:59");
  date2 = DateTime(date1);
  TEST_TRUE(date1 == date2)
}
END_SECTION

START_SECTION((DateTime(const DateTime&& date)))
{
  // Ensure that DateTime has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(DateTime(std::declval<DateTime&&>())), true)

  DateTime date1;
  DateTime date2;
  DateTime date3;

  date1.set("2006-12-12 11:59:59");
  date2 = DateTime(std::move(date1));
  TEST_EQUAL(date2.get(), "2006-12-12 11:59:59")
}
END_SECTION

START_SECTION((DateTime& operator= (const DateTime& source)))
{
  DateTime date, date2;
  date.set("2006-12-12 11:59:59");
  TEST_EQUAL(date==date2,false);
  date2 = date;
  TEST_EQUAL(date==date2,true);
}
END_SECTION

START_SECTION((DateTime& operator= (DateTime&& source)))
{
  DateTime date, date2;
  date.set("2006-12-12 11:59:59");
  TEST_EQUAL(date==date2,false);
  date2 = std::move(date);
  TEST_EQUAL(date2.get(), "2006-12-12 11:59:59")
}
END_SECTION

START_SECTION((void clear()))
{
  DateTime date1;
  DateTime date2;
  date1.set("2006-12-12 11:59:59");
  date1.clear();
  TEST_TRUE(date1 == date2)
  TEST_EQUAL(date1.isNull(), true)
}
END_SECTION

START_SECTION((String get() const))
{
  DateTime date_time;
  date_time.set("1999-11-24 14:24:31");
  TEST_EQUAL(date_time.get(),"1999-11-24 14:24:31")
}
END_SECTION

START_SECTION((void get(UInt& month, UInt& day, UInt& year, UInt& hour, UInt& minute, UInt& second) const))
{
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
}
END_SECTION

START_SECTION((void getDate(UInt& month, UInt& day, UInt& year) const))
{
  DateTime date;
  UInt month;
  UInt day;
  UInt year;

  date.set("2006-12-14 21:12:02");

  date.getDate(month, day, year);
  TEST_EQUAL(month, 12)  
  TEST_EQUAL(day, 14)  
  TEST_EQUAL(year, 2006)  
}
END_SECTION

START_SECTION((String getDate() const))
{
  DateTime date;
  date.set("2006-12-14 21:12:02");
  TEST_STRING_EQUAL(date.getDate(), String("2006-12-14"))  
}
END_SECTION

START_SECTION((void getTime(UInt& hour, UInt& minute, UInt& second) const))
{
  DateTime date;
  UInt hour; 
  UInt minute; 
  UInt second;

  date.set("2006-12-14 11:59:58");

  date.getTime(hour, minute, second);
  TEST_EQUAL(hour, 11)  
  TEST_EQUAL(minute, 59)  
  TEST_EQUAL(second, 58)    
}
END_SECTION

START_SECTION((String getTime() const))
{
  DateTime date;
  date.set("2006-12-14 11:59:58");
  TEST_STRING_EQUAL(date.getTime(), "11:59:58")    
}
END_SECTION

START_SECTION((void set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second)))
{
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
}
END_SECTION

START_SECTION((void set(const String &date)))
{
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

}
END_SECTION

START_SECTION((void setDate(UInt month, UInt day, UInt year)))
{
  DateTime date;
  UInt month = 12;
  UInt day = 14;
  UInt year = 2006;

  date.setDate(month, day, year);

  date.getDate(month, day, year);
  TEST_EQUAL(month, 12)  
  TEST_EQUAL(day, 14)  
  TEST_EQUAL(year, 2006)  
}
END_SECTION

START_SECTION((void setDate(const String &date)))
{
  DateTime date;
  UInt month;
  UInt day;
  UInt year;

  date.set("2006-12-14 11:59:58");

  date.getDate(month, day, year);
  TEST_EQUAL(month, 12)  
  TEST_EQUAL(day, 14)  
  TEST_EQUAL(year, 2006)  
}
END_SECTION

START_SECTION((void setTime(UInt hour, UInt minute, UInt second)))
{
  DateTime date;
  UInt hour; 
  UInt minute; 
  UInt second;

  date.setTime(11, 59, 58);

  date.getTime(hour, minute, second);
  TEST_EQUAL(hour, 11)  
  TEST_EQUAL(minute, 59)  
  TEST_EQUAL(second, 58)    
}
END_SECTION

START_SECTION((void setTime(const String &date)))
{
  DateTime date;
  UInt hour; 
  UInt minute; 
  UInt second;

  date.setTime("11:59:58");

  date.getTime(hour, minute, second);
  TEST_EQUAL(hour, 11)  
  TEST_EQUAL(minute, 59)  
  TEST_EQUAL(second, 58)    
}
END_SECTION

START_SECTION(([EXTRA] Three digit year should get leading zero according to Qt 4.4.3 documentation ))
{
  // This is a regression test.  Leave it here even if the issue gets hacked away in DateTime.
  DateTime one_moment_in_time;
  one_moment_in_time.set(5,4,666,3,2,1);

  // this behaviour is not critical and does not work on Qt 4.3 machines
  // so the leading zero is not checked! (who really needs dates before the year 1000 in this library?) 
  TEST_EQUAL(one_moment_in_time.get().hasSubstring("666-05-04 03:02:01"), true);
}
END_SECTION

START_SECTION((static DateTime now()))
{
  TEST_EQUAL(DateTime::now().isValid(), true)
}
END_SECTION

START_SECTION((DateTime& addSecs(int s) - backward traversal past year 0))
{
  // Test that subtracting seconds that would result in year < 0 throws an exception
  DateTime date;
  date.set(1, 1, 1, 0, 0, 0);  // January 1, year 1, 00:00:00

  // Subtracting 1 second should move us to year 0, which should throw
  TEST_EXCEPTION(Exception::ParseError, date.addSecs(-1))

  // Test with a larger negative offset
  DateTime date2;
  date2.set(1, 1, 5, 0, 0, 0);  // January 1, year 5, 00:00:00

  // Subtracting enough seconds to go back more than 5 years should throw
  int days_to_subtract = 365 * 6;  // More than 5 years
  int seconds_to_subtract = -1 * days_to_subtract * 24 * 3600;
  TEST_EXCEPTION(Exception::ParseError, date2.addSecs(seconds_to_subtract))
}
END_SECTION

START_SECTION((DateTime& addSecs(int s) - backward traversal near year 0))
{
  // Test that we can safely traverse backward to year 1
  DateTime date;
  date.set(1, 5, 2, 12, 0, 0);  // January 5, year 2, 12:00:00

  // Subtract 5 days (should be safe)
  date.addSecs(-5 * 24 * 3600);
  UInt month, day, year, hour, minute, second;
  date.get(month, day, year, hour, minute, second);
  TEST_EQUAL(year, 1)  // Should be in year 1
  TEST_EQUAL(month, 12)  // December
  TEST_EQUAL(day, 31)  // December 31

  // Test that we can move to the earliest valid moment (year 1, month 1, day 1)
  DateTime date2;
  date2.set(1, 1, 1, 0, 0, 1);  // January 1, year 1, 00:00:01
  date2.addSecs(-1);  // Go back 1 second to January 1, year 1, 00:00:00
  date2.get(month, day, year, hour, minute, second);
  TEST_EQUAL(year, 1)
  TEST_EQUAL(month, 1)
  TEST_EQUAL(day, 1)
  TEST_EQUAL(hour, 0)
  TEST_EQUAL(minute, 0)
  TEST_EQUAL(second, 0)
}
END_SECTION

START_SECTION((DateTime& addSecs(int s) - forward and backward traversal))
{
  // Test normal forward and backward operations still work correctly
  DateTime date;
  date.set(6, 15, 2020, 14, 30, 45);  // June 15, 2020, 14:30:45

  // Add 1 day
  date.addSecs(24 * 3600);
  TEST_EQUAL(date.get(), "2020-06-16 14:30:45")

  // Subtract 2 days
  date.addSecs(-2 * 24 * 3600);
  TEST_EQUAL(date.get(), "2020-06-14 14:30:45")

  // Test month boundary crossing
  DateTime date2;
  date2.set(3, 1, 2020, 0, 0, 0);  // March 1, 2020, 00:00:00 (leap year)
  date2.addSecs(-1);  // Should go to Feb 29, 2020, 23:59:59
  UInt month, day, year, hour, minute, second;
  date2.get(month, day, year, hour, minute, second);
  TEST_EQUAL(year, 2020)
  TEST_EQUAL(month, 2)
  TEST_EQUAL(day, 29)  // Leap year
  TEST_EQUAL(hour, 23)
  TEST_EQUAL(minute, 59)
  TEST_EQUAL(second, 59)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
