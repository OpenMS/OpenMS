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

#include <cctype>
#include <chrono>
#include <cstdio>
#include <cstring>
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

  // helper: validate time fields
  static bool isValidTime_(int hour, int minute, int second)
  {
    return hour >= 0 && hour <= 23 && minute >= 0 && minute <= 59 && second >= 0 && second <= 59;
  }

  // helper: convert fields to a time_point, add seconds, convert back
  // Uses std::chrono year_month_day and hh_mm_ss (C++20)
  static void addSecsToFields_(int& year, int& month, int& day,
                               int& hour, int& minute, int& second, int secs)
  {
    // Build a tm struct and use mktime/gmtime approach via chrono
    struct tm t{};
    t.tm_year = year - 1900;
    t.tm_mon = month - 1;
    t.tm_mday = day;
    t.tm_hour = hour;
    t.tm_min = minute;
    t.tm_sec = second;
    t.tm_isdst = -1; // let mktime figure it out -- actually we want naive arithmetic

    // Use timegm (POSIX) or _mkgmtime (MSVC) to avoid timezone issues
    // since our DateTime is naive (no timezone state)
#ifdef _WIN32
    time_t tt = _mkgmtime(&t);
#else
    time_t tt = timegm(&t);
#endif

    tt += secs;

    struct tm result{};
#ifdef _WIN32
    _gmtime64_s(&result, &tt);
#else
    gmtime_r(&tt, &result);
#endif

    year = result.tm_year + 1900;
    month = result.tm_mon + 1;
    day = result.tm_mday;
    hour = result.tm_hour;
    minute = result.tm_min;
    second = result.tm_sec;
  }

  // helper: parse 3-letter month abbreviation to month number (1-12), returns 0 on failure
  static int parseMonthAbbrev_(const char* s)
  {
    static const char* months[] = {
      "Jan", "Feb", "Mar", "Apr", "May", "Jun",
      "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };
    for (int i = 0; i < 12; ++i)
    {
      if (strncmp(s, months[i], 3) == 0)
      {
        return i + 1;
      }
    }
    return 0;
  }

  DateTime::DateTime() = default;

  bool DateTime::operator==(const DateTime& rhs) const
  {
    return std::tie(fields_.year, fields_.month, fields_.day,
                    fields_.hour, fields_.minute, fields_.second, fields_.millisecond, fields_.valid)
        == std::tie(rhs.fields_.year, rhs.fields_.month, rhs.fields_.day,
                    rhs.fields_.hour, rhs.fields_.minute, rhs.fields_.second, rhs.fields_.millisecond, rhs.fields_.valid);
  }

  bool DateTime::operator!=(const DateTime& rhs) const
  {
    return !(*this == rhs);
  }

  bool DateTime::operator<(const DateTime& rhs) const
  {
    return std::tie(fields_.year, fields_.month, fields_.day,
                    fields_.hour, fields_.minute, fields_.second, fields_.millisecond)
         < std::tie(rhs.fields_.year, rhs.fields_.month, rhs.fields_.day,
                    rhs.fields_.hour, rhs.fields_.minute, rhs.fields_.second, rhs.fields_.millisecond);
  }

  bool DateTime::isValid() const
  {
    return fields_.valid;
  }

  String DateTime::toString(const std::string& format) const
  {
    if (!fields_.valid) return String();

    char buf[64];

    if (format == "yyyy-MM-ddThh:mm:ss")
    {
      snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02d",
               fields_.year, fields_.month, fields_.day,
               fields_.hour, fields_.minute, fields_.second);
    }
    else if (format == "yyyy-MM-ddThh:mm:ss.zzz")
    {
      snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02d.%03d",
               fields_.year, fields_.month, fields_.day,
               fields_.hour, fields_.minute, fields_.second, fields_.millisecond);
    }
    else if (format == "yyyy-MM-dd hh:mm:ss")
    {
      snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d",
               fields_.year, fields_.month, fields_.day,
               fields_.hour, fields_.minute, fields_.second);
    }
    else if (format == "yyyy-MM-dd+hh:mm")
    {
      snprintf(buf, sizeof(buf), "%04d-%02d-%02d+%02d:%02d",
               fields_.year, fields_.month, fields_.day,
               fields_.hour, fields_.minute);
    }
    else if (format == "yyyy-MM-ddThh:mm:ssZ")
    {
      snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02dZ",
               fields_.year, fields_.month, fields_.day,
               fields_.hour, fields_.minute, fields_.second);
    }
    else if (format == "yyyy-MM-dd")
    {
      snprintf(buf, sizeof(buf), "%04d-%02d-%02d",
               fields_.year, fields_.month, fields_.day);
    }
    else if (format == "hh:mm:ss")
    {
      snprintf(buf, sizeof(buf), "%02d:%02d:%02d",
               fields_.hour, fields_.minute, fields_.second);
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Unknown DateTime format string", format);
    }

    return String(buf);
  }

  void DateTime::set(const String& date)
  {
    clear();

    int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0, millisecond = 0;
    bool parsed = false;

    if (date.has('.') && !date.has('T'))
    {
      // German format: dd.MM.yyyy hh:mm:ss
      if (sscanf(date.c_str(), "%d.%d.%d %d:%d:%d", &day, &month, &year, &hour, &minute, &second) == 6)
      {
        parsed = true;
      }
    }
    else if (date.has('/'))
    {
      // US format: MM/dd/yyyy hh:mm:ss
      if (sscanf(date.c_str(), "%d/%d/%d %d:%d:%d", &month, &day, &year, &hour, &minute, &second) == 6)
      {
        parsed = true;
      }
    }
    else if (date.has('-'))
    {
      if (date.has('T'))
      {
        if (date.has('+'))
        {
          // ISO with timezone suffix -- strip at '+'
          String stripped = date.prefix('+');
          if (stripped.has('.'))
          {
            // with milliseconds
            if (sscanf(stripped.c_str(), "%d-%d-%dT%d:%d:%d.%d", &year, &month, &day, &hour, &minute, &second, &millisecond) == 7)
            {
              // Normalize fractional seconds: e.g. ".4" -> 400, ".46" -> 460, ".468" -> 468
              auto dot_pos = stripped.find('.');
              if (dot_pos != String::npos)
              {
                int frac_digits = 0;
                for (size_t i = dot_pos + 1; i < stripped.size() && std::isdigit(static_cast<unsigned char>(stripped[i])); ++i) ++frac_digits;
                for (int i = frac_digits; i < 3; ++i) millisecond *= 10;
                millisecond %= 1000;
              }
              parsed = true;
            }
          }
          else
          {
            if (sscanf(stripped.c_str(), "%d-%d-%dT%d:%d:%d", &year, &month, &day, &hour, &minute, &second) == 6)
            {
              parsed = true;
            }
          }
        }
        else if (date.has('.'))
        {
          // ISO 8601 with milliseconds, no timezone: yyyy-MM-ddThh:mm:ss.zzz
          if (sscanf(date.c_str(), "%d-%d-%dT%d:%d:%d.%d", &year, &month, &day, &hour, &minute, &second, &millisecond) == 7)
          {
            // Normalize fractional seconds: e.g. ".4" -> 400, ".46" -> 460, ".468" -> 468
            auto dot_pos = date.find('.');
            if (dot_pos != String::npos)
            {
              int frac_digits = 0;
              for (size_t i = dot_pos + 1; i < date.size() && std::isdigit(static_cast<unsigned char>(date[i])); ++i) ++frac_digits;
              for (int i = frac_digits; i < 3; ++i) millisecond *= 10;
              millisecond %= 1000;
            }
            parsed = true;
          }
        }
        else
        {
          // ISO 8601: yyyy-MM-ddThh:mm:ss
          if (sscanf(date.c_str(), "%d-%d-%dT%d:%d:%d", &year, &month, &day, &hour, &minute, &second) == 6)
          {
            parsed = true;
          }
        }
      }
      else if (date.has('Z'))
      {
        // Legacy literal: yyyy-MM-ddZ  (Z is discarded, NOT UTC)
        if (sscanf(date.c_str(), "%d-%d-%dZ", &year, &month, &day) == 3)
        {
          // Verify there's nothing unexpected after Z
          // Check: the format is exactly "yyyy-MM-ddZ" -- if there's more after Z, it's invalid
          // Find Z position
          auto zpos = date.find('Z');
          if (zpos != std::string::npos && zpos + 1 == date.size())
          {
            parsed = true;
          }
        }
      }
      else if (date.has('+'))
      {
        // Legacy literal: yyyy-MM-dd+hh:mm  (+ is separator, NOT timezone)
        if (sscanf(date.c_str(), "%d-%d-%d+%d:%d", &year, &month, &day, &hour, &minute) == 5)
        {
          parsed = true;
        }
      }
      else
      {
        // ISO with space: yyyy-MM-dd hh:mm:ss
        if (sscanf(date.c_str(), "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &minute, &second) == 6)
        {
          parsed = true;
        }
      }
    }

    // Fallback: Legacy protXML format "ddd MMM d YYYY" e.g. "Tue Jul 13 2004"
    // or "ddd MMM dd hh:mm:ss yyyy" e.g. "Thu Dec 14 11:45:26 2006"
    if (!parsed)
    {
      // Try: "Xxx Xxx dd hh:mm:ss yyyy"
      char day_name[4] = {};
      char month_name[4] = {};
      int d2 = 0, h2 = 0, m2 = 0, s2 = 0, y2 = 0;
      if (sscanf(date.c_str(), "%3s %3s %d %d:%d:%d %d", day_name, month_name, &d2, &h2, &m2, &s2, &y2) == 7)
      {
        int mon = parseMonthAbbrev_(month_name);
        if (mon > 0)
        {
          year = y2; month = mon; day = d2;
          hour = h2; minute = m2; second = s2;
          parsed = true;
        }
      }
      else if (sscanf(date.c_str(), "%3s %3s %d %d", day_name, month_name, &d2, &y2) == 4)
      {
        int mon = parseMonthAbbrev_(month_name);
        if (mon > 0)
        {
          year = y2; month = mon; day = d2;
          parsed = true;
        }
      }
    }

    if (!parsed || !isValidDate_(year, month, day) || !isValidTime_(hour, minute, second))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Invalid date time string");
    }

    fields_.year = year;
    fields_.month = month;
    fields_.day = day;
    fields_.hour = hour;
    fields_.minute = minute;
    fields_.second = second;
    fields_.millisecond = millisecond;
    fields_.valid = true;
  }

  void DateTime::set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second)
  {
    if (!isValidDate_((int)year, (int)month, (int)day) || !isValidTime_((int)hour, (int)minute, (int)second))
    {
      String date_time = String(year) + "-" + String(month) + "-" + String(day)
                         + " " + String(hour) + ":" + String(minute) + ":" + String(second);
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date_time, "Invalid date time");
    }

    fields_.year = (int)year;
    fields_.month = (int)month;
    fields_.day = (int)day;
    fields_.hour = (int)hour;
    fields_.minute = (int)minute;
    fields_.second = (int)second;
    fields_.millisecond = 0;
    fields_.valid = true;
  }

  DateTime DateTime::now()
  {
    auto tp = chrono::system_clock::now();
    time_t tt = chrono::system_clock::to_time_t(tp);
    struct tm local{};
#ifdef _WIN32
    localtime_s(&local, &tt);
#else
    localtime_r(&tt, &local);
#endif

    DateTime d;
    d.fields_.year = local.tm_year + 1900;
    d.fields_.month = local.tm_mon + 1;
    d.fields_.day = local.tm_mday;
    d.fields_.hour = local.tm_hour;
    d.fields_.minute = local.tm_min;
    d.fields_.second = local.tm_sec;
    d.fields_.millisecond = 0;
    d.fields_.valid = true;
    return d;
  }

  DateTime DateTime::nowUTC()
  {
    auto tp = chrono::system_clock::now();
    time_t tt = chrono::system_clock::to_time_t(tp);
    struct tm utc{};
#ifdef _WIN32
    _gmtime64_s(&utc, &tt);
#else
    gmtime_r(&tt, &utc);
#endif

    DateTime d;
    d.fields_.year = utc.tm_year + 1900;
    d.fields_.month = utc.tm_mon + 1;
    d.fields_.day = utc.tm_mday;
    d.fields_.hour = utc.tm_hour;
    d.fields_.minute = utc.tm_min;
    d.fields_.second = utc.tm_sec;
    d.fields_.millisecond = 0;
    d.fields_.valid = true;
    return d;
  }

  String DateTime::get() const
  {
    if (fields_.valid)
    {
      return toString("yyyy-MM-dd hh:mm:ss");
    }
    return "0000-00-00 00:00:00";
  }

  void DateTime::get(UInt& month, UInt& day, UInt& year,
                     UInt& hour, UInt& minute, UInt& second) const
  {
    year = fields_.year;
    month = fields_.month;
    day = fields_.day;
    hour = fields_.hour;
    minute = fields_.minute;
    second = fields_.second;
  }

  void DateTime::clear()
  {
    fields_ = Fields{};
  }

  void DateTime::setDate(const String& date)
  {
    int year = 0, month = 0, day = 0;
    bool parsed = false;

    if (date.has('-'))
    {
      if (sscanf(date.c_str(), "%d-%d-%d", &year, &month, &day) == 3)
      {
        parsed = true;
      }
    }
    else if (date.has('.'))
    {
      if (sscanf(date.c_str(), "%d.%d.%d", &day, &month, &year) == 3)
      {
        parsed = true;
      }
    }
    else if (date.has('/'))
    {
      if (sscanf(date.c_str(), "%d/%d/%d", &month, &day, &year) == 3)
      {
        parsed = true;
      }
    }

    if (!parsed || !isValidDate_(year, month, day))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Could not set date");
    }

    fields_.year = year;
    fields_.month = month;
    fields_.day = day;
    fields_.valid = true;
  }

  void DateTime::setTime(const String& time)
  {
    int hour = 0, minute = 0, second = 0;

    if (sscanf(time.c_str(), "%d:%d:%d", &hour, &minute, &second) != 3 || !isValidTime_(hour, minute, second))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, time, "Could not set time");
    }

    fields_.hour = hour;
    fields_.minute = minute;
    fields_.second = second;
    // If we're setting time on a previously null DateTime, mark as valid
    // (matches Qt behavior where setTime on a null QDateTime makes it valid with date 0)
    fields_.valid = true;
  }

  void DateTime::setDate(UInt month, UInt day, UInt year)
  {
    if (!isValidDate_((int)year, (int)month, (int)day))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  String(year) + "-" + String(month) + "-" + String(day), "Could not set date");
    }

    fields_.year = (int)year;
    fields_.month = (int)month;
    fields_.day = (int)day;
    fields_.valid = true;
  }

  void DateTime::setTime(UInt hour, UInt minute, UInt second)
  {
    if (!isValidTime_((int)hour, (int)minute, (int)second))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  String(hour) + ":" + String(minute) + ":" + String(second), "Could not set time");
    }

    fields_.hour = (int)hour;
    fields_.minute = (int)minute;
    fields_.second = (int)second;
    fields_.valid = true;
  }

  void DateTime::getDate(UInt& month, UInt& day, UInt& year) const
  {
    month = fields_.month;
    day = fields_.day;
    year = fields_.year;
  }

  String DateTime::getDate() const
  {
    if (fields_.valid)
    {
      return toString("yyyy-MM-dd");
    }
    return "0000-00-00";
  }

  void DateTime::getTime(UInt& hour, UInt& minute, UInt& second) const
  {
    hour = fields_.hour;
    minute = fields_.minute;
    second = fields_.second;
  }

  String DateTime::getTime() const
  {
    if (fields_.valid)
    {
      return toString("hh:mm:ss");
    }
    return "00:00:00";
  }

  DateTime& DateTime::addSecs(int s)
  {
    addSecsToFields_(fields_.year, fields_.month, fields_.day,
                     fields_.hour, fields_.minute, fields_.second, s);
    return *this;
  }

  bool DateTime::isNull() const
  {
    return !fields_.valid;
  }

  // static
  DateTime DateTime::fromString(const std::string& date, const std::string& format)
  {
    DateTime d;
    int year = 0, month = 0, day = 0, hour = 0, minute = 0, second = 0, millisecond = 0;
    int n = 0;

    if (format == "yyyy-MM-ddThh:mm:ss")
    {
      n = sscanf(date.c_str(), "%d-%d-%dT%d:%d:%d", &year, &month, &day, &hour, &minute, &second);
      if (n != 6) return d; // return invalid
    }
    else if (format == "yyyy-MM-ddThh:mm:ss.zzz")
    {
      n = sscanf(date.c_str(), "%d-%d-%dT%d:%d:%d.%d", &year, &month, &day, &hour, &minute, &second, &millisecond);
      if (n != 7) return d;
      // Normalize fractional seconds: e.g. ".4" -> 400, ".46" -> 460, ".468" -> 468
      auto dot_pos = date.find('.');
      if (dot_pos != std::string::npos)
      {
        int frac_digits = 0;
        for (size_t i = dot_pos + 1; i < date.size() && std::isdigit(static_cast<unsigned char>(date[i])); ++i) ++frac_digits;
        for (int i = frac_digits; i < 3; ++i) millisecond *= 10;
        millisecond %= 1000;
      }
    }
    else if (format == "yyyy-MM-dd hh:mm:ss")
    {
      n = sscanf(date.c_str(), "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &minute, &second);
      if (n != 6) return d;
    }
    else if (format == "yyyy-MM-dd+hh:mm")
    {
      n = sscanf(date.c_str(), "%d-%d-%d+%d:%d", &year, &month, &day, &hour, &minute);
      if (n != 5) return d;
    }
    else if (format == "yyyy-MM-ddThh:mm:ssZ")
    {
      n = sscanf(date.c_str(), "%d-%d-%dT%d:%d:%dZ", &year, &month, &day, &hour, &minute, &second);
      if (n != 6) return d;
    }
    else if (format == "yyyy-MM-dd")
    {
      n = sscanf(date.c_str(), "%d-%d-%d", &year, &month, &day);
      if (n != 3) return d;
    }
    else if (format == "hh:mm:ss")
    {
      n = sscanf(date.c_str(), "%d:%d:%d", &hour, &minute, &second);
      if (n != 3) return d;
    }
    else
    {
      return d; // return invalid for unknown format
    }

    if (!isValidDate_(year, month, day) && format != "hh:mm:ss")
    {
      return d;
    }
    if (!isValidTime_(hour, minute, second))
    {
      return d;
    }

    d.fields_.year = year;
    d.fields_.month = month;
    d.fields_.day = day;
    d.fields_.hour = hour;
    d.fields_.minute = minute;
    d.fields_.second = second;
    d.fields_.millisecond = millisecond;
    d.fields_.valid = true;
    return d;
  }

} // namespace OpenMS
