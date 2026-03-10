// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/HashUtils.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <string>

namespace OpenMS
{

  /**
      @brief DateTime Class (Qt-free).

      This class implements date/time handling without Qt.
      Import and export to/from both string and integers is possible.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI DateTime
  {
  public:
    DateTime();
    DateTime(const DateTime& date);
    DateTime(DateTime&&) noexcept;
    DateTime& operator=(const DateTime& source);
    DateTime& operator=(DateTime&&) & noexcept;
    ~DateTime();

    bool operator==(const DateTime& rhs) const;
    bool operator!=(const DateTime& rhs) const;
    bool operator<(const DateTime& rhs) const;

    // Date and time setters
    void setDate(const String& date);
    void setTime(const String& time);
    void setDate(UInt month, UInt day, UInt year);
    void setTime(UInt hour, UInt minute, UInt second);
    void set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second);

    // Getters
    void get(UInt& month, UInt& day, UInt& year, UInt& hour, UInt& minute, UInt& second) const;
    void getDate(UInt& month, UInt& day, UInt& year) const;
    String getDate() const;
    void getTime(UInt& hour, UInt& minute, UInt& second) const;
    String getTime() const;

    // Add seconds (can be negative)
    DateTime& addSecs(int s);

    // Now
    static DateTime now();
    static DateTime nowUTC();

    // State
    bool isValid() const;
    bool isNull() const;
    void clear();

    // Formatting
    String toString(const std::string& format = "yyyy-MM-ddThh:mm:ss") const;
    static DateTime fromString(const std::string& date, const std::string& format = "yyyy-MM-ddThh:mm:ss");
    String get() const; // "yyyy-MM-dd hh:mm:ss"

    // Full set from combined string
    void set(const String& date);

  private:
    static bool validDate_(int y, int m, int d);
    static bool validTime_(int hh, int mm, int ss);
    static void parseDate_(const String& s, int& y, int& m, int& d);
    static void parseTime_(const String& s, int& hh, int& mm, int& ss);

    int year_{0};
    int month_{0};
    int day_{0};
    int hour_{0};
    int minute_{0};
    int second_{0};
    bool valid_{false};
  };

} // namespace OpenMS

// Hash function specialization for DateTime
namespace std
{
  template<>
  struct hash<OpenMS::DateTime>
  {
    std::size_t operator()(const OpenMS::DateTime& dt) const noexcept
    {
      // Hash date/time components to match operator==
      std::string datetime_str = dt.toString("yyyy-MM-ddThh:mm:ss");
      return OpenMS::fnv1a_hash_string(datetime_str);
    }
  };
} // namespace std
