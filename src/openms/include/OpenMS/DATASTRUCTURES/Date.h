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

namespace OpenMS
{
  /**
      @brief Date Class.

      This class implements date handling.
      Import and export to/from both string and integers is possible.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI Date
  {
public:

    /**
        @brief Default constructor

        Fills the object with an undefined date: 00/00/0000
    */
    Date();

    /// Copy constructor
    Date(const Date& date) = default;

    /// Move constructor
    Date(Date&&) = default;

    /// Assignment operator
    Date& operator=(const Date& source) = default;

    /// Move assignment operator
    Date& operator=(Date&&) & = default;

    /// equal operator
    bool operator==(const Date& rhs) const;

    /**
        @brief sets data from a string

        The following date formats are supported:
        - mm/dd/yyyy
        - dd.mm.yyyy
        - yyyy-mm-dd

        @exception Exception::ParseError is thrown if the date is given in the wrong format
    */
    void set(const String& date);

    /**
        @brief sets data from three integers

        @exception Exception::ParseError is thrown if an invalid date is given
    */
    void set(UInt month, UInt day, UInt year);

    /// Returns the current date
    static Date today();

    /**
        @brief Returns a string representation of the date

        Uses the iso/ansi date format: 'yyyy-mm-dd'
    */
    String get() const;

    /**
        @brief Fills the arguments with the date

        Give the numbers in the following order: month, day and year.
    */
    void get(UInt& month, UInt& day, UInt& year) const;

    ///Sets the undefined date: 00/00/0000
    void clear();

    /// Checks if the date is valid
    bool isValid() const;

private:
    UInt year_;
    UInt month_;
    UInt day_;
  };
} // namespace OPENMS

