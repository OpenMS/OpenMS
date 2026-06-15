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
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <iostream>
#include <vector>
#include <clocale>
#include <cstdlib>
#include <ctime>

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

START_SECTION((std::string get() const))
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

START_SECTION((std::string getDate() const))
{
  DateTime date;
  date.set("2006-12-14 21:12:02");
  TEST_STRING_EQUAL(date.getDate(),std::string("2006-12-14"))  
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

START_SECTION((std::string getTime() const))
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

START_SECTION((void set(const std::string &date)))
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

START_SECTION((void setDate(const std::string &date)))
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

START_SECTION((void setTime(const std::string &date)))
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
  TEST_EQUAL(StringUtils::hasSubstring(one_moment_in_time.get(), "666-05-04 03:02:01"), true);
}
END_SECTION

START_SECTION((static DateTime now()))
{
  TEST_EQUAL(DateTime::now().isValid(), true)
}
END_SECTION

START_SECTION([EXTRA] locale- and timezone-independent parsing and formatting)
{
  // DateTime is naive: set()/get()/toString()/fromString() parse and format via
  // sscanf("%d", ...) (locale-independent integer parsing), strncmp month matching,
  // and timegm/gmtime (UTC) for arithmetic -- never localtime or locale-sensitive
  // number/string functions. The output must therefore be byte-identical regardless
  // of the C locale (e.g. tr_TR's dotted-i, de_DE's decimal comma) and the TZ
  // environment variable. The other sections never vary either; this pins that
  // contract. (now() is excluded -- it legitimately returns local wall-clock time.)
  auto check_invariant = []()
  {
    DateTime dt;

    dt.set("1999-11-24 14:24:31");
    TEST_EQUAL(dt.get(), "1999-11-24 14:24:31")

    dt.set("01.02.2000 14:24:32");           // dd.mm.yyyy
    TEST_EQUAL(dt.get(), "2000-02-01 14:24:32")

    dt.set("2005-11-13T10:58:57");           // ISO 'T' separator
    TEST_EQUAL(dt.get(), "2005-11-13 10:58:57")

    dt.set("2011-08-05T15:32:07.468+02:00"); // millisecond + timezone offset, both dropped
    TEST_EQUAL(dt.get(), "2011-08-05 15:32:07")

    dt.set("Wed Dec 14 11:59:58 2006");      // asctime-style; exercises month-name matching
    TEST_EQUAL(dt.get(), "2006-12-14 11:59:58")

    DateTime dt2 = DateTime::fromString("2020-03-09T08:07:06", "yyyy-MM-ddThh:mm:ss");
    TEST_EQUAL(dt2.toString("yyyy-MM-dd hh:mm:ss"), "2020-03-09 08:07:06")
  };

  // Baseline: whatever locale/TZ the test process started in.
  check_invariant();

#ifndef OPENMS_WINDOWSPLATFORM
  // --- vary TZ (POSIX setenv/tzset; tzdata is essentially always present) ---
  const char* old_tz = std::getenv("TZ");
  const std::string saved_tz = (old_tz != nullptr) ? std::string(old_tz) : std::string();
  const bool had_tz = (old_tz != nullptr);

  for (const char* tz : {"UTC", "America/New_York", "Asia/Kolkata", "Pacific/Chatham", "Etc/GMT-14"})
  {
    setenv("TZ", tz, 1);
    tzset();
    check_invariant();
  }

  if (had_tz) { setenv("TZ", saved_tz.c_str(), 1); }
  else        { unsetenv("TZ"); }
  tzset();

  // --- vary the C locale (guarded: unavailable locales are simply skipped, so the
  // section still passes on minimal images that lack de_DE / tr_TR) ---
  const char* cur = setlocale(LC_ALL, nullptr);
  const std::string saved_loc = (cur != nullptr) ? std::string(cur) : std::string("C");

  for (const char* loc : {"C", "POSIX", "de_DE.UTF-8", "tr_TR.UTF-8", "de_DE.utf8", "tr_TR.utf8"})
  {
    if (setlocale(LC_ALL, loc) != nullptr)
    {
      check_invariant();
    }
  }
  setlocale(LC_ALL, saved_loc.c_str());
#endif
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
