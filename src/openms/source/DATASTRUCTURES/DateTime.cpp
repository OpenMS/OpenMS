// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>

using namespace std;

namespace OpenMS
{
  namespace
  {
    // Helper functions (same as in Date.cpp)
    bool isLeapYear(UInt year)
    {
      return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
    }

    UInt daysInMonth(UInt month, UInt year)
    {
      if (month < 1 || month > 12)
      {
        throw std::out_of_range("Month must be in range 1-12, got: " + std::to_string(month));
      }
      static const UInt days[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
      if (month == 2 && isLeapYear(year))
      {
        return 29;
      }
      return days[month];
    }

    bool isValidDate(UInt year, UInt month, UInt day)
    {
      if (year == 0 && month == 0 && day == 0)
      {
        return false;
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

    bool isValidTime(UInt hour, UInt minute, UInt second)
    {
      return hour < 24 && minute < 60 && second < 60;
    }

    bool isValidDateTime(UInt year, UInt month, UInt day, UInt hour, UInt minute, UInt second)
    {
      return isValidDate(year, month, day) && isValidTime(hour, minute, second);
    }
  }

  DateTime::DateTime() :
    year_(0),
    month_(0),
    day_(0),
    hour_(0),
    minute_(0),
    second_(0)
  {
  }

  DateTime::DateTime(const DateTime& date) = default;

  DateTime::DateTime(DateTime&& rhs) noexcept = default;

  DateTime& DateTime::operator=(const DateTime& source) = default;

  DateTime& DateTime::operator=(DateTime&& source) & noexcept = default;

  DateTime::~DateTime() = default;

  bool DateTime::operator==(const DateTime& rhs) const
  {
    return year_ == rhs.year_ && month_ == rhs.month_ && day_ == rhs.day_ &&
           hour_ == rhs.hour_ && minute_ == rhs.minute_ && second_ == rhs.second_;
  }

  bool DateTime::operator!=(const DateTime& rhs) const
  {
    return !(*this == rhs);
  }

  bool DateTime::operator<(const DateTime& rhs) const
  {
    if (year_ != rhs.year_) return year_ < rhs.year_;
    if (month_ != rhs.month_) return month_ < rhs.month_;
    if (day_ != rhs.day_) return day_ < rhs.day_;
    if (hour_ != rhs.hour_) return hour_ < rhs.hour_;
    if (minute_ != rhs.minute_) return minute_ < rhs.minute_;
    return second_ < rhs.second_;
  }

  bool DateTime::isValid() const
  {
    return isValidDateTime(year_, month_, day_, hour_, minute_, second_);
  }

  String DateTime::toString(const std::string& format) const
  {
    String fmt = format;
    String result = fmt;

    // Replace format tokens
    String year_str = String(year_);
    while (year_str.size() < 4) year_str = "0" + year_str;
    result.substitute("yyyy", year_str);

    String month_str = String(month_);
    if (month_str.size() < 2) month_str = "0" + month_str;
    result.substitute("MM", month_str);

    String day_str = String(day_);
    if (day_str.size() < 2) day_str = "0" + day_str;
    result.substitute("dd", day_str);

    String hour_str = String(hour_);
    if (hour_str.size() < 2) hour_str = "0" + hour_str;
    result.substitute("hh", hour_str);

    String minute_str = String(minute_);
    if (minute_str.size() < 2) minute_str = "0" + minute_str;
    result.substitute("mm", minute_str);

    String second_str = String(second_);
    if (second_str.size() < 2) second_str = "0" + second_str;
    result.substitute("ss", second_str);

    return result;
  }

  void DateTime::set(const String& date)
  {
    clear();

    UInt y = 0, mo = 0, d = 0, h = 0, mi = 0, s = 0;
    bool parsed = false;

    // Try different formats
    std::istringstream iss(date);
    char sep1, sep2, sep3, sep4;

    // Format: dd.MM.yyyy hh:mm:ss
    if (date.has('.') && !date.has('T'))
    {
      iss >> d >> sep1 >> mo >> sep2 >> y >> h >> sep3 >> mi >> sep4 >> s;
      parsed = !iss.fail() && sep1 == '.' && sep2 == '.' && sep3 == ':' && sep4 == ':';
    }
    // Format: MM/dd/yyyy hh:mm:ss
    else if (date.has('/'))
    {
      iss >> mo >> sep1 >> d >> sep2 >> y >> h >> sep3 >> mi >> sep4 >> s;
      parsed = !iss.fail() && sep1 == '/' && sep2 == '/' && sep3 == ':' && sep4 == ':';
    }
    // Format: yyyy-MM-ddThh:mm:ss or yyyy-MM-dd hh:mm:ss
    else if (date.has('-'))
    {
      if (date.has('T'))
      {
        // Handle timezone: yyyy-MM-ddThh:mm:ss+hh:mm
        String date_trimmed = date;
        if (date.has('+'))
        {
          date_trimmed = date.prefix('+');
        }
        // Handle milliseconds: yyyy-MM-ddThh:mm:ss.zzz
        if (date_trimmed.has('.'))
        {
          size_t dot_pos = date_trimmed.find('.');
          date_trimmed = date_trimmed.substr(0, dot_pos);
        }

        iss.str(date_trimmed);
        char sep5;
        iss >> y >> sep1 >> mo >> sep2 >> d >> sep3 >> h >> sep4 >> mi >> sep5 >> s;
        parsed = !iss.fail() && sep1 == '-' && sep2 == '-' && sep3 == 'T' && sep4 == ':' && sep5 == ':';
      }
      else if (date.has('Z'))
      {
        // Format: yyyy-MM-ddZ (date only with Z timezone)
        char extra;
        iss >> y >> sep1 >> mo >> sep2 >> d >> sep3;
        // Check that sep3 is 'Z' and there's nothing after it
        parsed = !iss.fail() && sep1 == '-' && sep2 == '-' && sep3 == 'Z';
        if (parsed && iss >> extra)
        {
          // There's something after 'Z', which is invalid
          parsed = false;
        }
        h = mi = s = 0;
      }
      else if (date.has('+') && !date.has('T'))
      {
        // Format: yyyy-MM-dd+hh:mm ('+' acts as time separator)
        iss >> y >> sep1 >> mo >> sep2 >> d >> sep3 >> h >> sep4 >> mi;
        parsed = !iss.fail() && sep1 == '-' && sep2 == '-' && sep3 == '+' && sep4 == ':';
        s = 0;
      }
      else
      {
        // Format: yyyy-MM-dd hh:mm:ss
        iss >> y >> sep1 >> mo >> sep2 >> d >> h >> sep3 >> mi >> sep4 >> s;
        parsed = !iss.fail() && sep1 == '-' && sep2 == '-' && sep3 == ':' && sep4 == ':';
      }
    }

    if (!parsed || !isValidDateTime(y, mo, d, h, mi, s))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Invalid date time string");
    }

    year_ = y;
    month_ = mo;
    day_ = d;
    hour_ = h;
    minute_ = mi;
    second_ = s;
  }

  void DateTime::set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second)
  {
    if (!isValidDateTime(year, month, day, hour, minute, second))
    {
      String date_time = String(year) + "-" + String(month) + "-" + String(day)
                         + " " + String(hour) + ":" + String(minute) + ":" + String(second);
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date_time, "Invalid date time");
    }

    year_ = year;
    month_ = month;
    day_ = day;
    hour_ = hour;
    minute_ = minute;
    second_ = second;
  }

  DateTime DateTime::now()
  {
    DateTime d;

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
    d.hour_ = local_tm.tm_hour;
    d.minute_ = local_tm.tm_min;
    d.second_ = local_tm.tm_sec;

    return d;
  }

  String DateTime::get() const
  {
    if (isValid())
    {
      std::ostringstream oss;
      oss << std::setfill('0')
          << std::setw(4) << year_ << "-"
          << std::setw(2) << month_ << "-"
          << std::setw(2) << day_ << " "
          << std::setw(2) << hour_ << ":"
          << std::setw(2) << minute_ << ":"
          << std::setw(2) << second_;
      return oss.str();
    }
    return "0000-00-00 00:00:00";
  }

  void DateTime::get(UInt& month, UInt& day, UInt& year,
                     UInt& hour, UInt& minute, UInt& second) const
  {
    year = year_;
    month = month_;
    day = day_;
    hour = hour_;
    minute = minute_;
    second = second_;
  }

  void DateTime::clear()
  {
    year_ = 0;
    month_ = 0;
    day_ = 0;
    hour_ = 0;
    minute_ = 0;
    second_ = 0;
  }

  void DateTime::setDate(const String& date)
  {
    UInt y = 0, mo = 0, d = 0;
    bool parsed = false;
    std::istringstream iss(date);
    char sep1, sep2;

    if (date.has('-'))
    {
      iss >> y >> sep1 >> mo >> sep2 >> d;
      parsed = !iss.fail() && sep1 == '-' && sep2 == '-';
    }
    else if (date.has('.'))
    {
      iss >> d >> sep1 >> mo >> sep2 >> y;
      parsed = !iss.fail() && sep1 == '.' && sep2 == '.';
    }
    else if (date.has('/'))
    {
      iss >> mo >> sep1 >> d >> sep2 >> y;
      parsed = !iss.fail() && sep1 == '/' && sep2 == '/';
    }
    else
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Could not set date");
    }

    if (!parsed || !isValidDate(y, mo, d))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Could not set date");
    }

    year_ = y;
    month_ = mo;
    day_ = d;
  }

  void DateTime::setTime(const String& time)
  {
    UInt h = 0, mi = 0, s = 0;
    std::istringstream iss(time);
    char sep1, sep2;

    iss >> h >> sep1 >> mi >> sep2 >> s;
    if (iss.fail() || sep1 != ':' || sep2 != ':' || !isValidTime(h, mi, s))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, time, "Could not set time");
    }

    hour_ = h;
    minute_ = mi;
    second_ = s;
  }

  void DateTime::setDate(UInt month, UInt day, UInt year)
  {
    if (!isValidDate(year, month, day))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(year) + "-" + String(month) + "-" + String(day), "Could not set date");
    }

    year_ = year;
    month_ = month;
    day_ = day;
  }

  void DateTime::setTime(UInt hour, UInt minute, UInt second)
  {
    if (!isValidTime(hour, minute, second))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(hour) + ":" + String(minute) + ":" + String(second), "Could not set time");
    }

    hour_ = hour;
    minute_ = minute;
    second_ = second;
  }

  void DateTime::getDate(UInt& month, UInt& day, UInt& year) const
  {
    month = month_;
    day = day_;
    year = year_;
  }

  String DateTime::getDate() const
  {
    if (isValidDate(year_, month_, day_))
    {
      std::ostringstream oss;
      oss << std::setfill('0')
          << std::setw(4) << year_ << "-"
          << std::setw(2) << month_ << "-"
          << std::setw(2) << day_;
      return oss.str();
    }
    return "0000-00-00";
  }

  void DateTime::getTime(UInt& hour, UInt& minute, UInt& second) const
  {
    hour = hour_;
    minute = minute_;
    second = second_;
  }

  String DateTime::getTime() const
  {
    if (isValidTime(hour_, minute_, second_))
    {
      std::ostringstream oss;
      oss << std::setfill('0')
          << std::setw(2) << hour_ << ":"
          << std::setw(2) << minute_ << ":"
          << std::setw(2) << second_;
      return oss.str();
    }
    return "00:00:00";
  }

  DateTime& DateTime::addSecs(int s)
  {
    // Convert to total seconds since start of day
    int total_seconds = hour_ * 3600 + minute_ * 60 + second_ + s;

    // Handle negative values and days overflow
    int days_offset = 0;
    while (total_seconds < 0)
    {
      total_seconds += 24 * 3600;
      days_offset--;
    }
    while (total_seconds >= 24 * 3600)
    {
      total_seconds -= 24 * 3600;
      days_offset++;
    }

    // Update time
    hour_ = total_seconds / 3600;
    minute_ = (total_seconds % 3600) / 60;
    second_ = total_seconds % 60;

    // Update date with proper month/year overflow handling
    if (days_offset != 0)
    {
      // Use signed arithmetic to properly handle negative day values
      int day_signed = static_cast<int>(day_) + days_offset;

      // Handle day overflow (forward in time)
      while (day_signed > static_cast<int>(daysInMonth(month_, year_)))
      {
        day_signed -= daysInMonth(month_, year_);
        month_++;
        if (month_ > 12)
        {
          month_ = 1;
          year_++;
        }
      }

      // Handle day underflow (backward in time)
      while (day_signed < 1)
      {
        month_--;
        if (month_ < 1)
        {
          // Check for year underflow before decrementing
          if (year_ <= 1)
          {
            throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Date arithmetic resulted in year underflow (year < 1)",
                                      "Cannot represent dates before year 1");
          }
          month_ = 12;
          year_--;
        }
        day_signed += daysInMonth(month_, year_);
      }

      // Assign back to unsigned after all arithmetic is done
      day_ = static_cast<UInt>(day_signed);
    }

    return *this;
  }

  bool DateTime::isNull() const
  {
    return year_ == 0 && month_ == 0 && day_ == 0 && hour_ == 0 && minute_ == 0 && second_ == 0;
  }

  // static
  DateTime DateTime::fromString(const std::string& date, const std::string& format)
  {
    DateTime d;

    // Simple implementation using the format string as a template
    // This is simplified - real implementation would need full format parsing
    if (format == "yyyy-MM-ddThh:mm:ss")
    {
      d.set(date);
    }
    else
    {
      // Fallback
      d.set(date);
    }

    return d;
  }

} // namespace OpenMS
