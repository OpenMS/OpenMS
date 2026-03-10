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

#include <ctime>
#include <iomanip>
#include <sstream>
#include <tuple>
#include <cctype>

namespace OpenMS
{
  // small helpers for parsing fixed width integers
  static bool parse_int2_(const char* p, int& v)
  {
    if (!std::isdigit(static_cast<unsigned char>(p[0])) || !std::isdigit(static_cast<unsigned char>(p[1]))) return false;
    v = (p[0] - '0') * 10 + (p[1] - '0');
    return true;
  }

  static bool parse_int4_(const char* p, int& v)
  {
    for (int i = 0; i < 4; ++i)
    {
      if (!std::isdigit(static_cast<unsigned char>(p[i]))) return false;
    }
    v = (p[0]-'0')*1000 + (p[1]-'0')*100 + (p[2]-'0')*10 + (p[3]-'0');
    return true;
  }

  DateTime::DateTime()
  {
    clear();
  }

  DateTime::DateTime(const DateTime& date) = default;
  DateTime::DateTime(DateTime&&) noexcept = default;
  DateTime& DateTime::operator=(const DateTime& source) = default;
  DateTime& DateTime::operator=(DateTime&&) & noexcept = default;
  DateTime::~DateTime() = default;

  bool DateTime::operator==(const DateTime& rhs) const
  {
    return std::tie(year_, month_, day_, hour_, minute_, second_, valid_)
         == std::tie(rhs.year_, rhs.month_, rhs.day_, rhs.hour_, rhs.minute_, rhs.second_, rhs.valid_);
  }

  bool DateTime::operator!=(const DateTime& rhs) const
  {
    return !(*this == rhs);
  }

  bool DateTime::operator<(const DateTime& rhs) const
  {
    return std::tie(year_, month_, day_, hour_, minute_, second_)
         < std::tie(rhs.year_, rhs.month_, rhs.day_, rhs.hour_, rhs.minute_, rhs.second_);
  }

  bool DateTime::validDate_(int y, int m, int d)
  {
    if (y <= 0 || m < 1 || m > 12 || d < 1) return false;
    static const int mdays[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    const bool isLeap = ((y % 4 == 0) && (y % 100 != 0)) || (y % 400 == 0);
    const int dim = mdays[m-1] + ((m == 2 && isLeap) ? 1 : 0);
    return d <= dim;
  }

  bool DateTime::validTime_(int hh, int mm, int ss)
  {
    return (hh >= 0 && hh < 24) && (mm >= 0 && mm < 60) && (ss >= 0 && ss < 60);
  }

  void DateTime::parseDate_(const String& s, int& y, int& m, int& d)
  {
    // yyyy-MM-dd (strict length and separators)
    if (s.has('-'))
    {
      if (s.size() != 10) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid date string");
      if (!parse_int4_(s.c_str(), y)) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid year");
      if (s[4] != '-') throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid date separator");
      if (!parse_int2_(s.c_str()+5, m)) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid month");
      if (s[7] != '-') throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid date separator");
      if (!parse_int2_(s.c_str()+8, d)) throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid day");
    }
    // MM/dd/yyyy (strict length and separators)
    else if (s.has('/'))
    {
      if (s.size() != 10 || s[2] != '/' || s[5] != '/')
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid date");
      int mm=0, dd=0, yy=0;
      if (std::sscanf(s.c_str(), "%d/%d/%d", &mm, &dd, &yy) != 3)
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid date");
      m = mm; d = dd; y = yy;
    }
    // dd.MM.yyyy (strict length and separators)
    else if (s.has('.'))
    {
      if (s.size() != 10 || s[2] != '.' || s[5] != '.')
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid date");
      int mm=0, dd=0, yy=0;
      if (std::sscanf(s.c_str(), "%d.%d.%d", &dd, &mm, &yy) != 3)
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid date");
      m = mm; d = dd; y = yy;
    }
    else
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Unsupported date format");
    }

    if (!validDate_(y, m, d))
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid date");
  }

  void DateTime::parseTime_(const String& s, int& hh, int& mm, int& ss)
  {
    if (std::sscanf(s.c_str(), "%d:%d:%d", &hh, &mm, &ss) != 3)
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid time");
    if (!validTime_(hh, mm, ss))
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, s, "Invalid time");
  }

  void DateTime::setDate(const String& date)
  {
    // Treat sentinel "0000-00-00" as 'unset' (keep invalid state)
    if (date == "0000-00-00")
    {
      clear();
      return;
    }
    int y=0, m=0, d=0;
    parseDate_(date, y, m, d);
    year_ = y; month_ = m; day_ = d;
    valid_ = true;
  }

  void DateTime::setTime(const String& time)
  {
    int hh=0, mm=0, ss=0;
    parseTime_(time, hh, mm, ss);
    hour_ = hh; minute_ = mm; second_ = ss;
    valid_ = true;
  }

  void DateTime::setDate(UInt month, UInt day, UInt year)
  {
    int y = static_cast<int>(year);
    int m = static_cast<int>(month);
    int d = static_cast<int>(day);
    if (!validDate_(y, m, d))
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(year) + "-" + String(month) + "-" + String(day), "Could not set date");
    year_ = y; month_ = m; day_ = d; valid_ = true;
  }

  void DateTime::setTime(UInt hour, UInt minute, UInt second)
  {
    int hh = static_cast<int>(hour);
    int mm = static_cast<int>(minute);
    int ss = static_cast<int>(second);
    if (!validTime_(hh, mm, ss))
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(hour) + ":" + String(minute) + ":" + String(second), "Could not set time");
    hour_ = hh; minute_ = mm; second_ = ss; valid_ = true;
  }

  void DateTime::set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second)
  {
    setDate(month, day, year);
    setTime(hour, minute, second);
  }

  void DateTime::get(UInt& month, UInt& day, UInt& year, UInt& hour, UInt& minute, UInt& second) const
  {
    month = static_cast<UInt>(valid_ ? month_ : 0);
    day   = static_cast<UInt>(valid_ ? day_   : 0);
    year  = static_cast<UInt>(valid_ ? year_  : 0);
    hour  = static_cast<UInt>(valid_ ? hour_  : 0);
    minute= static_cast<UInt>(valid_ ? minute_: 0);
    second= static_cast<UInt>(valid_ ? second_: 0);
  }

  void DateTime::getDate(UInt& month, UInt& day, UInt& year) const
  {
    month = static_cast<UInt>(valid_ ? month_ : 0);
    day   = static_cast<UInt>(valid_ ? day_   : 0);
    year  = static_cast<UInt>(valid_ ? year_  : 0);
  }

  String DateTime::getDate() const
  {
    if (!valid_) return "0000-00-00";
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(4) << year_ << '-'
        << std::setw(2) << month_ << '-'
        << std::setw(2) << day_;
    return String(oss.str());
  }

  void DateTime::getTime(UInt& hour, UInt& minute, UInt& second) const
  {
    hour   = static_cast<UInt>(valid_ ? hour_   : 0);
    minute = static_cast<UInt>(valid_ ? minute_ : 0);
    second = static_cast<UInt>(valid_ ? second_ : 0);
  }

  String DateTime::getTime() const
  {
    if (!valid_) return "00:00:00";
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << hour_ << ':'
        << std::setw(2) << minute_ << ':'
        << std::setw(2) << second_;
    return String(oss.str());
  }

  DateTime& DateTime::addSecs(int s)
  {
    // Build tm
    std::tm tm{};
    tm.tm_year = (year_ > 0 ? year_ - 1900 : 0);
    tm.tm_mon  = (month_ > 0 ? month_ - 1 : 0);
    tm.tm_mday = (day_ > 0 ? day_ : 1);
    tm.tm_hour = hour_;
    tm.tm_min  = minute_;
    tm.tm_sec  = second_;
    std::time_t t = std::mktime(&tm);
    if (t == static_cast<std::time_t>(-1))
    {
      // if invalid current state, start from epoch
      t = 0;
    }
    t += s;
#if defined(_WIN32)
    std::tm ntm{};
    localtime_s(&ntm, &t);
    year_  = ntm.tm_year + 1900;
    month_ = ntm.tm_mon + 1;
    day_   = ntm.tm_mday;
    hour_  = ntm.tm_hour;
    minute_= ntm.tm_min;
    second_= ntm.tm_sec;
#else
    std::tm* ntm = std::localtime(&t);
    year_  = ntm->tm_year + 1900;
    month_ = ntm->tm_mon + 1;
    day_   = ntm->tm_mday;
    hour_  = ntm->tm_hour;
    minute_= ntm->tm_min;
    second_= ntm->tm_sec;
#endif
    valid_ = true;
    return *this;
  }

  DateTime DateTime::now()
  {
    DateTime d;
    std::time_t t = std::time(nullptr);
#if defined(_WIN32)
    std::tm lt{};
    localtime_s(&lt, &t);
    d.year_   = lt.tm_year + 1900;
    d.month_  = lt.tm_mon + 1;
    d.day_    = lt.tm_mday;
    d.hour_   = lt.tm_hour;
    d.minute_ = lt.tm_min;
    d.second_ = lt.tm_sec;
#else
    std::tm* lt = std::localtime(&t);
    d.year_   = lt->tm_year + 1900;
    d.month_  = lt->tm_mon + 1;
    d.day_    = lt->tm_mday;
    d.hour_   = lt->tm_hour;
    d.minute_ = lt->tm_min;
    d.second_ = lt->tm_sec;
#endif
    d.valid_ = true;
    return d;
  }

  bool DateTime::isValid() const { return valid_; }
  bool DateTime::isNull() const  { return !valid_; }

  void DateTime::clear()
  {
    year_ = month_ = day_ = hour_ = minute_ = second_ = 0;
    valid_ = false;
  }

  String DateTime::toString(const std::string& format) const
  {
    // Minimal support: "yyyy-MM-ddThh:mm:ss" and default "yyyy-MM-dd hh:mm:ss"
    // Note: "yyyy-MM-dd+hh:mm" is NOT a valid Qt format, so we fall back to default like Qt did
    if (!valid_) return "0000-00-00 00:00:00";
    std::ostringstream oss;
    if (format == "yyyy-MM-ddThh:mm:ss")
    {
      oss << std::setfill('0') << std::setw(4) << year_ << '-'
          << std::setw(2) << month_ << '-'
          << std::setw(2) << day_ << 'T'
          << std::setw(2) << hour_ << ':'
          << std::setw(2) << minute_ << ':'
          << std::setw(2) << second_;
      return String(oss.str());
    }
    else if (format == "yyyy-MM-dd+hh:mm") // modified ISO 8601 format used in legacy OpenMS
    {
      // For any unsupported format, fall back to default like Qt did
      oss << std::setfill('0') << std::setw(4) << year_ << '-'
          << std::setw(2) << month_ << '-'
          << std::setw(2) << day_ << '+'
          << std::setw(2) << hour_ << ':'
          << std::setw(2) << minute_ ;
      return String(oss.str());
    }
    else // default "yyyy-MM-dd hh:mm:ss"
    {
      oss << std::setfill('0') << std::setw(4) << year_ << '-'
          << std::setw(2) << month_ << '-'
          << std::setw(2) << day_ << ' '
          << std::setw(2) << hour_ << ':'
          << std::setw(2) << minute_ << ':'
          << std::setw(2) << second_;
      return String(oss.str());
    }
  }

  DateTime DateTime::fromString(const std::string& date, const std::string& format)
  {
    DateTime d;
    if (format == "yyyy-MM-ddThh:mm:ss")
    {
      if (date.size() < 19) return d;
      d.set(String(date.substr(0, 10) + " " + date.substr(11, 8)));
      return d;
    }
    // default try to set() directly
    d.set(String(date));
    return d;
  }

  DateTime DateTime::nowUTC()
  {
    DateTime d;
    std::time_t t = std::time(nullptr);
#if defined(_WIN32)
    std::tm gt{};
    gmtime_s(&gt, &t);
    d.year_   = gt.tm_year + 1900;
    d.month_  = gt.tm_mon + 1;
    d.day_    = gt.tm_mday;
    d.hour_   = gt.tm_hour;
    d.minute_ = gt.tm_min;
    d.second_ = gt.tm_sec;
#else
    std::tm* gt = std::gmtime(&t);
    d.year_   = gt->tm_year + 1900;
    d.month_  = gt->tm_mon + 1;
    d.day_    = gt->tm_mday;
    d.hour_   = gt->tm_hour;
    d.minute_ = gt->tm_min;
    d.second_ = gt->tm_sec;
#endif
    d.valid_ = true;
    return d;
  }

  String DateTime::get() const
  {
    return toString("yyyy-MM-dd hh:mm:ss");
  }

  void DateTime::set(const String& date)
  {
    clear();
    // Normalize and handle sentinel "0000-00-00" forms as 'unset'
    String d = date;
    d = d.trim();
    if (d == "0000-00-00" || d == "0000-00-00 00:00:00")
    {
      // keep invalid state after clear()
      return;
    }
    // Supported patterns:
    // - "MM/dd/yyyy hh:mm:ss"
    // - "dd.MM.yyyy hh:mm:ss"
    // - "yyyy-MM-dd hh:mm:ss"
    // - "yyyy-MM-ddThh:mm:ss[.sss][+hh:mm]" (timezone/milliseconds ignored)
    // - "yyyy-MM-ddZ" (date-only)
    // - "yyyy-MM-dd+hh:mm" (treat '+hh:mm' as time-of-day)
    // - "EEE MMM d hh:mm:ss yyyy" (Unix ctime-like; weekday ignored; also accepts 3+ digit years)

    // Try Unix ctime-like format early to avoid failing on the generic ' ' split
    {
      const char* p = d.c_str();
      // Skip leading spaces
      while (*p == ' ') ++p;
      // Parse weekday token (>=3 letters)
      int letters = 0;
      while (std::isalpha(static_cast<unsigned char>(*p))) { ++p; ++letters; }
      bool has_weekday = (letters >= 3);

      // Proceed if weekday token and next non-space seems like a month token
      const char* q = p;
      while (*q == ' ') ++q;
      auto tolower_c = [](char ch) -> char { return static_cast<char>(std::tolower(static_cast<unsigned char>(ch))); };
      auto month_from = [&](char a, char b, char c) -> int
      {
        a = tolower_c(a); b = tolower_c(b); c = tolower_c(c);
        if (a=='j' && b=='a' && c=='n') return 1;
        if (a=='f' && b=='e' && c=='b') return 2;
        if (a=='m' && b=='a' && c=='r') return 3;
        if (a=='a' && b=='p' && c=='r') return 4;
        if (a=='m' && b=='a' && c=='y') return 5;
        if (a=='j' && b=='u' && c=='n') return 6;
        if (a=='j' && b=='u' && c=='l') return 7;
        if (a=='a' && b=='u' && c=='g') return 8;
        if (a=='s' && b=='e' && c=='p') return 9;
        if (a=='o' && b=='c' && c=='t') return 10;
        if (a=='n' && b=='o' && c=='v') return 11;
        if (a=='d' && b=='e' && c=='c') return 12;
        return 0;
      };

      if (has_weekday &&
          std::isalpha(static_cast<unsigned char>(q[0])) &&
          std::isalpha(static_cast<unsigned char>(q[1])) &&
          std::isalpha(static_cast<unsigned char>(q[2])) )
      {
        int mmonth = month_from(q[0], q[1], q[2]);
        if (mmonth != 0)
        {
          // Move p to after month
          p = q + 3;
          // Skip spaces
          while (*p == ' ') ++p;
          // Day of month (1-2 digits)
          int dd = 0;
          int nd = 0;
          while (std::isdigit(static_cast<unsigned char>(*p)))
          {
            dd = dd * 10 + (*p - '0');
            ++p; ++nd;
          }
          if (nd > 0)
          {
            // Skip spaces
            while (*p == ' ') ++p;
            // Time hh:mm:ss
            int hh = 0, mm = 0, ss = 0;
            int n = 0;
            while (std::isdigit(static_cast<unsigned char>(*p))) { hh = hh*10 + (*p - '0'); ++p; ++n; }
            if (n > 0 && *p == ':')
            {
              ++p; n = 0;
              while (std::isdigit(static_cast<unsigned char>(*p))) { mm = mm*10 + (*p - '0'); ++p; ++n; }
              if (n > 0 && *p == ':')
              {
                ++p; n = 0;
                while (std::isdigit(static_cast<unsigned char>(*p))) { ss = ss*10 + (*p - '0'); ++p; ++n; }
                if (n > 0)
                {
                  // Skip spaces
                  while (*p == ' ') ++p;
                  // Year (accept >=1 digits)
                  int yy = 0; n = 0;
                  while (std::isdigit(static_cast<unsigned char>(*p)))
                  {
                    yy = yy*10 + (*p - '0');
                    ++p; ++n;
                  }
                  if (n > 0)
                  {
                    // Validate and set
                    if (validDate_(yy, mmonth, dd) && validTime_(hh, mm, ss))
                    {
                      year_ = yy; month_ = mmonth; day_ = dd;
                      hour_ = hh; minute_ = mm; second_ = ss;
                      valid_ = true;
                      return;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Supported patterns:
    // - "MM/dd/yyyy hh:mm:ss"
    // - "dd.MM.yyyy hh:mm:ss"
    // - "yyyy-MM-dd hh:mm:ss"
    // - "yyyy-MM-ddThh:mm:ss[.sss][+hh:mm]" (timezone/milliseconds ignored)
    // - "yyyy-MM-ddZ" (date-only)
    // - "yyyy-MM-dd+hh:mm" (treat '+hh:mm' as time-of-day)

    if (d.has('T'))
    {
      // ISO "yyyy-MM-ddThh:mm:ss[...]" - ignore trailing fractional seconds and timezone
      String ds = d.prefix('T');
      String ts = String(d.c_str() + static_cast<int>(ds.size()) + 1);
      setDate(ds);
      setTime(ts);
      return;
    }

    // 'Z' indicates date-only, must be the final character
    if (d.has('Z'))
    {
      if (d.size() > 0 && d[static_cast<int>(d.size()) - 1] == 'Z')
      {
        setDate(d.prefix('Z'));
        setTime(String("00:00:00"));
        return;
      }
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, d, "Invalid timezone/date suffix");
    }

    // '+' between date and time-of-day: "yyyy-MM-dd+hh:mm"
    if (d.has('+'))
    {
      String ds = d.prefix('+');
      // consume the '+' and build the suffix string
      const char* rest = d.c_str() + static_cast<int>(ds.size()) + 1;
      String hm(rest);
      int hh = 0, mm = 0;
      if (std::sscanf(hm.c_str(), "%d:%d", &hh, &mm) != 2 || !validTime_(hh, mm, 0))
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, d, "Invalid time-of-day after '+'");
      }
      setDate(ds);
      setTime(static_cast<UInt>(hh), static_cast<UInt>(mm), 0);
      return;
    }

    // space separator "yyyy-MM-dd hh:mm:ss" or other
    if (d.has(' '))
    {
      String ds = d.prefix(' ');
      String ts = String(d.c_str() + static_cast<int>(ds.size()) + 1);
      setDate(ds);
      setTime(ts);
      return;
    }

    // Fallback: only date provided (strict formats)
    setDate(d);
    setTime(String("00:00:00"));
  }
} // namespace OpenMS
