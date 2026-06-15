// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/Date.h>
#include <iostream>
#include <vector>
#include <clocale>
#include <string>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Date, "$Id$")

/////////////////////////////////////////////////////////////

Date* s_ptr = nullptr;
Date* s_nullPointer= nullptr;

START_SECTION((Date()))
	s_ptr = new Date();
  TEST_NOT_EQUAL(s_ptr, s_nullPointer)
END_SECTION

START_SECTION(([EXTRA]~Date()))
	delete s_ptr;
END_SECTION

START_SECTION(([EXTRA] Date constructed from set()))
	Date d;
	d.set(12, 24, 1999);
	TEST_EQUAL(d.year(),1999)
	TEST_EQUAL(d.month(),12)
	TEST_EQUAL(d.day(),24)
END_SECTION

START_SECTION((void get(UInt& month, UInt& day, UInt& year) const))
  Date date;
  UInt d,m,y;
  date.set("2007-12-03");
  date.get(m,d,y);
  TEST_EQUAL(m,12);
  TEST_EQUAL(d,3);
  TEST_EQUAL(y,2007);
END_SECTION

START_SECTION((void set(UInt month, UInt day, UInt year) ))
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
END_SECTION

START_SECTION((Date& operator= (const Date& source)))
  Date date, date2;
  date.set(12,1,1977);
  TEST_EQUAL(date==date2,false);
	date2 = date;
	TEST_EQUAL(date==date2,true);
END_SECTION

START_SECTION((Date(const Date& date)))
  Date date;
  date.set(12,1,1977);
	Date date2(date);
	TEST_EQUAL(date==date2,true);
END_SECTION

START_SECTION((void set(const std::string& date) ))
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
END_SECTION

START_SECTION((std::string get() const))
  Date d;
  TEST_EQUAL(d.get(),"0000-00-00");
  d.set("11.12.1977");
  TEST_EQUAL(d.get(),"1977-12-11");
	d.set("02.01.1999");
  TEST_EQUAL(d.get(),"1999-01-02");
END_SECTION

START_SECTION((void clear()))
  Date d;
  d.set("11.12.1977");
  TEST_EQUAL(d.get(),"1977-12-11");
	d.clear();
  TEST_EQUAL(d.get(),"0000-00-00");
END_SECTION

START_SECTION((static Date today()))
  TEST_EQUAL(Date::today().isValid(), true)
END_SECTION

START_SECTION([EXTRA] locale-independent date parsing and formatting)
{
  // Date is naive and date-only: set() parses via sscanf("%d.%d.%d", ...) and get()
  // formats via snprintf("%04d-%02d-%02d", ...) -- both locale-independent integer
  // operations, with no locale-sensitive number/string functions. So parsing and
  // formatting must be byte-identical across C locales (e.g. tr_TR's dotted-i,
  // de_DE's decimal comma). The other sections never vary the locale. (There is no
  // TZ axis -- Date carries no time-of-day; today(), which reads local time, is
  // excluded.)
  auto check_invariant = []()
  {
    Date d;
    UInt mo, dy, yr;

    d.set("01.12.1977");          // dd.mm.yyyy (German)
    d.get(mo, dy, yr);
    TEST_EQUAL(mo, 12) TEST_EQUAL(dy, 1) TEST_EQUAL(yr, 1977)
    TEST_EQUAL(d.get(), "1977-12-01")

    d.set("12/01/1977");          // mm/dd/yyyy (English)
    TEST_EQUAL(d.get(), "1977-12-01")

    d.set("1967-12-23");          // ISO yyyy-mm-dd
    TEST_EQUAL(d.get(), "1967-12-23")
  };

  // Baseline: whatever locale the test process started in.
  check_invariant();

#ifndef OPENMS_WINDOWSPLATFORM
  // Vary the C locale (guarded: unavailable locales are simply skipped, so the
  // section still passes on minimal images that lack de_DE / tr_TR).
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
