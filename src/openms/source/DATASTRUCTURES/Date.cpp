// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/CONCEPT/Exception.h>

using namespace std;

namespace OpenMS
{

  Date::Date(const QDate& date) :
    QDate(date)
  {
  }

  void Date::set(const String& date)
  {
    clear();

    //check for format (german/english)
    if (date.has('.'))
    {
      QDate::operator=(QDate::fromString(date.c_str(), "dd.MM.yyyy"));
    }
    else if (date.has('/'))
    {
      QDate::operator=(QDate::fromString(date.c_str(), "MM/dd/yyyy"));
    }
    else if (date.has('-'))
    {
      QDate::operator=(QDate::fromString(date.c_str(), "yyyy-MM-dd"));
    }

    if (!isValid())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Is no valid german, english or iso date");
    }
  }

  void Date::set(UInt month, UInt day, UInt year)
  {
    if (!setDate(year, month, day))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(year) + "-" + String(month) + "-" + String(day), "Invalid date");
    }
  }

  Date Date::today()
  {
    return QDate::currentDate();
  }

  String Date::get() const
  {
    if (QDate::isValid())
    {
      return toString("yyyy-MM-dd");
    }
    return "0000-00-00";
  }

  void Date::get(UInt& month, UInt& day, UInt& year) const
  {
    day = QDate::day();
    month = QDate::month();
    year = QDate::year();
  }

  void Date::clear()
  {
    QDate::operator=(QDate());
  }

} // namespace OpenMS

