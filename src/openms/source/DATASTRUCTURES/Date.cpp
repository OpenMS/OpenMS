// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <ctime>
#include <iomanip>
#include <sstream>

namespace OpenMS
{
  Date::Date()
  {
    clear();
  }

  bool Date::isLeap_(int y)
  {
    return ((y % 4 == 0) && (y % 100 != 0)) || (y % 400 == 0);
  }

  bool Date::valid_(int y, int m, int d)
  {
    if (y <= 0 || m < 1 || m > 12 || d < 1) return false;
    static const int mdays[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    int dim = mdays[m - 1] + ((m == 2 && isLeap_(y)) ? 1 : 0);
    return d <= dim;
  }

  void Date::set(const String& date)
  {
    clear();
    // dd.MM.yyyy
    if (date.has('.'))
    {
      int d = 0, m = 0, y = 0;
      char extra; // to catch any extra characters
      int parsed = std::sscanf(date.c_str(), "%d.%d.%d%c", &d, &m, &y, &extra);
      if (parsed == 3 && valid_(y, m, d))
      {
        day_ = d; month_ = m; year_ = y; valid_flag_ = true;
        return;
      }
    }
    // MM/dd/yyyy
    else if (date.has('/'))
    {
      int d = 0, m = 0, y = 0;
      char extra; // to catch any extra characters
      int parsed = std::sscanf(date.c_str(), "%d/%d/%d%c", &m, &d, &y, &extra);
      if (parsed == 3 && valid_(y, m, d))
      {
        day_ = d; month_ = m; year_ = y; valid_flag_ = true;
        return;
      }
    }
    // yyyy-MM-dd
    else if (date.has('-'))
    {
      int d = 0, m = 0, y = 0;
      char extra; // to catch any extra characters
      int parsed = std::sscanf(date.c_str(), "%d-%d-%d%c", &y, &m, &d, &extra);
      if (parsed == 3 && valid_(y, m, d))
      {
        day_ = d; month_ = m; year_ = y; valid_flag_ = true;
        return;
      }
    }

    throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Invalid date format (supported: dd.MM.yyyy, MM/dd/yyyy, yyyy-MM-dd)");
  }

  void Date::set(UInt month, UInt day, UInt year)
  {
    if (!valid_(static_cast<int>(year), static_cast<int>(month), static_cast<int>(day)))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  String(year) + "-" + String(month) + "-" + String(day), "Invalid date");
    }
    year_ = static_cast<int>(year);
    month_ = static_cast<int>(month);
    day_ = static_cast<int>(day);
    valid_flag_ = true;
  }

  Date Date::today()
  {
    Date d;
    std::time_t t = std::time(nullptr);
    std::tm lt{};
#if defined(_WIN32)
    localtime_s(&lt, &t);
#else
    lt = *std::localtime(&t);
#endif
    d.year_ = lt.tm_year + 1900;
    d.month_ = lt.tm_mon + 1;
    d.day_ = lt.tm_mday;
    d.valid_flag_ = true;
    return d;
  }

  String Date::get() const
  {
    if (!valid_flag_) return "0000-00-00";
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(4) << year_ << '-'
        << std::setw(2) << month_ << '-'
        << std::setw(2) << day_;
    return String(oss.str());
  }

  void Date::get(UInt& month, UInt& day, UInt& year) const
  {
    month = static_cast<UInt>(valid_flag_ ? month_ : 0);
    day   = static_cast<UInt>(valid_flag_ ? day_   : 0);
    year  = static_cast<UInt>(valid_flag_ ? year_  : 0);
  }

  void Date::clear()
  {
    year_ = month_ = day_ = 0;
    valid_flag_ = false;
  }
  int Date::year() const
  {
    return valid_flag_ ? year_ : 0;
  }

  int Date::month() const
  {
    return valid_flag_ ? month_ : 0;
  }

  int Date::day() const
  {
    return valid_flag_ ? day_ : 0;
  }

  bool Date::isValid() const
  {
    return valid_flag_;
  }

  bool Date::operator==(const Date& rhs) const
  {
    if (!valid_flag_ && !rhs.valid_flag_) return true;
    return year_ == rhs.year_ && month_ == rhs.month_ && day_ == rhs.day_ && valid_flag_ == rhs.valid_flag_;
  }

  bool Date::operator!=(const Date& rhs) const
  {
    return !(*this == rhs);
  }

} // namespace OpenMS
