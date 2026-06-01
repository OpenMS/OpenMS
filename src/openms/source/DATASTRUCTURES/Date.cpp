// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <chrono>
#include <cstdio>
#include <ctime>
#include <tuple>

using namespace std;

namespace OpenMS
{
  // helper: validate date fields
  static bool isValidDate_(int year, int month, int day)
  {
    if (month < 1 || month > 12 || day < 1 || year < 1)
    {
      return false;
    }
    static const int days_in_month[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    int max_day = days_in_month[month];
    if (month == 2)
    {
      bool leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
      if (leap) max_day = 29;
    }
    return day <= max_day;
  }

  Date::Date() = default;

  bool Date::operator==(const Date& rhs) const
  {
    return std::tie(fields_.year, fields_.month, fields_.day, fields_.valid)
        == std::tie(rhs.fields_.year, rhs.fields_.month, rhs.fields_.day, rhs.fields_.valid);
  }

  bool Date::operator!=(const Date& rhs) const
  {
    return !(*this == rhs);
  }

  bool Date::operator<(const Date& rhs) const
  {
    return std::tie(fields_.year, fields_.month, fields_.day)
         < std::tie(rhs.fields_.year, rhs.fields_.month, rhs.fields_.day);
  }

  void Date::set(const String& date)
  {
    clear();

    int year = 0, month = 0, day = 0;
    bool parsed = false;

    int n = 0; // track how many chars consumed

    if (date.has('.'))
    {
      if (sscanf(date.c_str(), "%d.%d.%d%n", &day, &month, &year, &n) == 3
          && n == (int)date.size())
      {
        parsed = true;
      }
    }
    else if (date.has('/'))
    {
      if (sscanf(date.c_str(), "%d/%d/%d%n", &month, &day, &year, &n) == 3
          && n == (int)date.size())
      {
        parsed = true;
      }
    }
    else if (date.has('-'))
    {
      if (sscanf(date.c_str(), "%d-%d-%d%n", &year, &month, &day, &n) == 3
          && n == (int)date.size())
      {
        parsed = true;
      }
    }

    if (!parsed || !isValidDate_(year, month, day))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Is no valid german, english or iso date");
    }

    fields_.year = year;
    fields_.month = month;
    fields_.day = day;
    fields_.valid = true;
  }

  void Date::set(UInt month, UInt day, UInt year)
  {
    if (!isValidDate_((int)year, (int)month, (int)day))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  String(year) + "-" + String(month) + "-" + String(day), "Invalid date");
    }

    fields_.year = (int)year;
    fields_.month = (int)month;
    fields_.day = (int)day;
    fields_.valid = true;
  }

  Date Date::today()
  {
    auto tp = chrono::system_clock::now();
    time_t tt = chrono::system_clock::to_time_t(tp);
    struct tm local{};
#ifdef _WIN32
    localtime_s(&local, &tt);
#else
    localtime_r(&tt, &local);
#endif

    Date d;
    d.fields_.year = local.tm_year + 1900;
    d.fields_.month = local.tm_mon + 1;
    d.fields_.day = local.tm_mday;
    d.fields_.valid = true;
    return d;
  }

  String Date::get() const
  {
    if (fields_.valid)
    {
      char buf[16];
      snprintf(buf, sizeof(buf), "%04d-%02d-%02d", fields_.year, fields_.month, fields_.day);
      return String(buf);
    }
    return "0000-00-00";
  }

  void Date::get(UInt& month, UInt& day, UInt& year) const
  {
    day = fields_.day;
    month = fields_.month;
    year = fields_.year;
  }

  void Date::clear()
  {
    fields_ = Fields{};
  }

  bool Date::isValid() const
  {
    return fields_.valid;
  }

  bool Date::isNull() const
  {
    return !fields_.valid;
  }

  int Date::year() const
  {
    return fields_.year;
  }

  int Date::month() const
  {
    return fields_.month;
  }

  int Date::day() const
  {
    return fields_.day;
  }

} // namespace OpenMS
