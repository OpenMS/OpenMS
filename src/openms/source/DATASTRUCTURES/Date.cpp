// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <QtCore/QDate>

using namespace std;

namespace OpenMS
{

  Date::Date() :
    date_(make_unique<QDate>())
  {
  }

  Date::Date(const Date& date) :
    date_(date.date_ ? make_unique<QDate>(*date.date_) : make_unique<QDate>())
  {
  }

  Date::Date(const QDate& date) :
    date_(make_unique<QDate>(date))
  {
  }

  Date::Date(Date&& rhs) noexcept :
    date_(rhs.date_ ? std::move(rhs.date_) : make_unique<QDate>())
  {
  }

  Date::~Date() = default;

  Date& Date::operator=(const Date& source)
  {
    if (&source == this)
    {
      return *this;
    }

    if (source.date_ == nullptr)
    {
      // Source is in a 'moved-from' state; create a default date
      date_ = make_unique<QDate>();
    }
    else if (date_ == nullptr)
    { // *this is in a 'moved-from' state; we need to create a date_ object first
      date_ = make_unique<QDate>(*source.date_);
    }
    else
    {
      *date_ = *source.date_;
    }

    return *this;
  }

  Date& Date::operator=(Date&& source) & noexcept
  {
    if (&source == this)
    {
      return *this;
    }

    std::swap(date_, source.date_);

    return *this;
  }

  bool Date::operator==(const Date& rhs) const
  {
    return (*date_ == *rhs.date_);
  }

  bool Date::operator!=(const Date& rhs) const
  {
    return !(*this == rhs);
  }

  bool Date::operator<(const Date& rhs) const
  {
    return (*date_ < *rhs.date_);
  }

  void Date::set(const String& date)
  {
    clear();

    //check for format (german/english)
    if (date.has('.'))
    {
      *date_ = QDate::fromString(date.c_str(), "dd.MM.yyyy");
    }
    else if (date.has('/'))
    {
      *date_ = QDate::fromString(date.c_str(), "MM/dd/yyyy");
    }
    else if (date.has('-'))
    {
      *date_ = QDate::fromString(date.c_str(), "yyyy-MM-dd");
    }

    if (!date_->isValid())
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, date, "Is no valid german, english or iso date");
    }
  }

  void Date::set(UInt month, UInt day, UInt year)
  {
    if (!date_->setDate(year, month, day))
    {
      throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String(year) + "-" + String(month) + "-" + String(day), "Invalid date");
    }
  }

  Date Date::today()
  {
    return Date(QDate::currentDate());
  }

  String Date::get() const
  {
    if (date_->isValid())
    {
      return date_->toString("yyyy-MM-dd");
    }
    return "0000-00-00";
  }

  void Date::get(UInt& month, UInt& day, UInt& year) const
  {
    day = date_->day();
    month = date_->month();
    year = date_->year();
  }

  void Date::clear()
  {
    *date_ = QDate();
  }

  bool Date::isValid() const
  {
    return date_->isValid();
  }

  bool Date::isNull() const
  {
    return date_->isNull();
  }

  int Date::year() const
  {
    return date_->year();
  }

  int Date::month() const
  {
    return date_->month();
  }

  int Date::day() const
  {
    return date_->day();
  }

} // namespace OpenMS

