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
#include <ctime>
#include <sstream>
#include <iomanip>

using namespace std;

namespace OpenMS
{

  Date::Date() :
    year_(0),
    month_(0),
    day_(0)
  {
  }

  bool Date::operator==(const Date& rhs) const
  {
    return year_ == rhs.year_ && month_ == rhs.month_ && day_ == rhs.day_;
  }

  namespace
  {
    // Helper function to check if a year is a leap year
    bool isLeapYear(UInt year)
    {
      return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    }

    // Helper function to get days in a month
    UInt daysInMonth(UInt month, UInt year)
    {
      static const UInt days[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
      if (month == 2 && isLeapYear(year))
      {
        return 29;
      }
      return days[month];
    }

    // Helper function to validate a date
    bool isValidDate(UInt year, UInt month, UInt day)
    {
      if (year == 0 && month == 0 && day == 0)
      {
        return false; // undefined date
      }
      if (month < 1 || month > 12)
      {
        return false;
      }
      if (day < 1 || day > daysInMonth(month, year))
      {
        return false;
      }
      return true;
    }

    // Helper function to parse date from string
    bool parseDateString(const String& date, const String& format, UInt& year, UInt& month, UInt& day)
    {
      std::istringstream iss(date);
      char sep1, sep2;

      if (format == "dd.MM.yyyy")
      {
        iss >> day >> sep1 >> month >> sep2 >> year;
        // Verify parse succeeded, separators are correct, and entire string was consumed
        return !iss.fail() && sep1 == '.' && sep2 == '.' && iss.peek() == EOF;
      }
      else if (format == "MM/dd/yyyy")
      {
        iss >> month >> sep1 >> day >> sep2 >> year;
        return !iss.fail() && sep1 == '/' && sep2 == '/' && iss.peek() == EOF;
      }
      else if (format == "yyyy-MM-dd")
      {
        iss >> year >> sep1 >> month >> sep2 >> day;
        return !iss.fail() && sep1 == '-' && sep2 == '-' && iss.peek() == EOF;
      }
      return false;
    }
  }

  void Date::set(const String& date)
  {
    clear();

    UInt y, m, d;
    bool parsed = false;

    //check for format (german/english/iso)
    if (date.has('.'))
    {
      parsed = parseDateString(date, "dd.MM.yyyy", y, m, d);
    }
    else if (date.has('/'))
    {
      parsed = parseDateString(date, "MM/dd/yyyy", y, m, d);
    }
    else if (date.has('-'))
    {
      parsed = parseDateString(date, "yyyy-MM-dd", y, m, d);
    }

    if (!parsed || !isValidDate(y, m, d))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Is no valid german, english or iso date");
    }

    year_ = y;
    month_ = m;
    day_ = d;
  }

  void Date::set(UInt month, UInt day, UInt year)
  {
    if (!isValidDate(year, month, day))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(year) + "-" + String(month) + "-" + String(day), "Invalid date");
    }

    year_ = year;
    month_ = month;
    day_ = day;
  }

  Date Date::today()
  {
    Date d;

    // Use C++11 chrono to get current date
    auto now = std::chrono::system_clock::now();
    auto now_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm local_tm;

    #ifdef _WIN32
    localtime_s(&local_tm, &now_time_t);
    #else
    localtime_r(&now_time_t, &local_tm);
    #endif

    d.year_ = local_tm.tm_year + 1900;
    d.month_ = local_tm.tm_mon + 1;
    d.day_ = local_tm.tm_mday;

    return d;
  }

  String Date::get() const
  {
    if (isValid())
    {
      std::ostringstream oss;
      oss << std::setfill('0') << std::setw(4) << year_ << "-"
          << std::setw(2) << month_ << "-"
          << std::setw(2) << day_;
      return oss.str();
    }
    return "0000-00-00";
  }

  void Date::get(UInt& month, UInt& day, UInt& year) const
  {
    day = day_;
    month = month_;
    year = year_;
  }

  void Date::clear()
  {
    year_ = 0;
    month_ = 0;
    day_ = 0;
  }

  bool Date::isValid() const
  {
    return isValidDate(year_, month_, day_);
  }

} // namespace OpenMS
