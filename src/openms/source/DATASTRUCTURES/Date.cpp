// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/Date.h>
#include <ctime>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/config.h>

using namespace std;

namespace OpenMS
{
  Date::Date() :
    QDate()
  {

  }

  Date::Date(const Date & date) :
    QDate(date)
  {
  }

  Date::Date(const QDate & date) :
    QDate(date)
  {
  }

  Date & Date::operator=(const Date & source)
  {
    if (&source == this)
    {
      return *this;
    }
    QDate::operator=(source);

    return *this;
  }

  void Date::set(const String & date)
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
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, date, "Is no valid german, english or iso date");
    }
  }

  void Date::set(UInt month, UInt day, UInt year)
  {
    if (!setDate(year, month, day))
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, String(year) + "-" + String(month) + "-" + String(day), "Invalid date");
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

  void Date::get(UInt & month, UInt & day, UInt & year) const
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
