// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <cstdint>

namespace OpenMS
{
  /**
      @brief Date Class.

      This class implements date handling without Qt dependencies.
      Import and export to/from both string and integers is possible.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI Date
  {
  public:
    /// Default constructor (undefined date 00/00/0000)
    Date();

    /// Copy/move
    Date(const Date&) = default;
    Date(Date&&) = default;
    Date& operator=(const Date&) = default;
    Date& operator=(Date&&) & = default;

    /**
      @brief sets data from a string

      Supported formats:
      - mm/dd/yyyy
      - dd.mm.yyyy
      - yyyy-mm-dd

      @throws Exception::ParseError if the date is invalid
    */
    void set(const String& date);

    /**
      @brief sets data from three integers

      @throws Exception::ParseError if the date is invalid
    */
    void set(UInt month, UInt day, UInt year);

    /// Returns the current local date
    static Date today();

    /// Returns ISO string 'yyyy-MM-dd' or '0000-00-00' if undefined
    String get() const;

    /// Accessors (return 0 if undefined)
    int year() const;
    int month() const;
    int day() const;

    /// Comparison
    bool operator==(const Date& rhs) const;
    bool operator!=(const Date& rhs) const;

    /// Validity
    bool isValid() const;

    /// Fills the arguments with month/day/year (0 if undefined)
    void get(UInt& month, UInt& day, UInt& year) const;

    /// Sets the undefined date: 00/00/0000
    void clear();

  private:
    static bool isLeap_(int y);
    static bool valid_(int y, int m, int d);

    int year_{0};
    int month_{0};
    int day_{0};
    bool valid_flag_{false};
  };
} // namespace OpenMS
