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
// $Maintainer: Sandro Andreotti $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_DATETIME_H
#define OPENMS_DATASTRUCTURES_DATETIME_H

#include <OpenMS/CONCEPT/Types.h>
#include <QtCore/QDate>

#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  class String;

  /**
      @brief DateTime Class.

      This class implements date handling.
      Import and export to/from both string and integers is possible.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI DateTime :
    public QDateTime
  {
public:

    /**
        @brief Default constructor

        Fills the object with an undefined date: 00/00/0000
    */
    DateTime();
    /// Copy constructor
    DateTime(const DateTime & date);
    /// Copy constructor from Qt base class
    DateTime(const QDateTime & date);

    /// Assignment operator
    DateTime & operator=(const DateTime & source);

    /**
        @brief sets date from a string

        Reads both English, German and iso/ansi date formats: 'MM/dd/yyyy', 'dd.MM.yyyy' or 'yyyy-MM-dd'

        @exception Exception::ParseError
    */
    void setDate(const String & date);

    /**
        @brief sets time from a string

        Reads time format: 'hh:mm:ss'

        @exception Exception::ParseError
    */
    void setTime(const String & date);

    /**
        @brief sets data from three integers

        Give the numbers in the following order: month, day and year.

        @exception Exception::ParseError
    */
    void setDate(UInt month, UInt day, UInt year);

    /**
        @brief sets time from three integers

        Give the numbers in the following order: hour, minute and second.

        @exception Exception::ParseError
    */
    void setTime(UInt hour, UInt minute, UInt second);

    /**
        @brief sets data from six integers

        Give the numbers in the following order: month, day, year, hour, minute, second.

        @exception Exception::ParseError
    */
    void set(UInt month, UInt day, UInt year, UInt hour, UInt minute, UInt second);

    /**
        @brief Fills the arguments with the date and the time

        Give the numbers in the following order: month, day and year, hour minute, second.
    */
    void get(UInt & month, UInt & day, UInt & year, UInt & hour, UInt & minute, UInt & second) const;

    /**
        @brief Fills the arguments with the date

        Give the numbers in the following order: month, day and year.
    */
    void getDate(UInt & month, UInt & day, UInt & year) const;

    /**
        @brief Returns the date as string

        The format of the string is yyyy-MM-dd
    */
    String getDate() const;

    /**
        @brief Fills the arguments with the time

        The arguments are all UInts and the order is hour minute second
    */
    void getTime(UInt & hour, UInt & minute, UInt & second) const;

    /**
        @brief Returns the time as string

        The format of the string is hh:mm:ss
    */
    String getTime() const;

    /// Returns the current date and time
    static DateTime now();

    ///Sets the undefined date: 00/00/0000 00:00:00
    void clear();

    /**
        @brief Returns a string representation of the date and time

        The format of the string will be yyyy-MM-dd hh:mm:ss
    */
    String get() const;

    /**
        @brief Sets date and time

        The following formats are supported:
        - MM/dd/yyyy hh:mm:ss
        - dd.MM.yyyy hh:mm:ss
        - yyyy-MM-dd hh:mm:ss
        - yyyy-MM-ddThh:mm:ss (ISO 8601 format)
        - yyyy-MM-ddZ (ISO 8601 format)
        - yyyy-MM-dd+hh:mm (ISO 8601 format)

        @exception Exception::ParseError
    */
    void set(const String & date);

protected:
  };

} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DATETIME_H
