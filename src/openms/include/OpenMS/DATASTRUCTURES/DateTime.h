// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

#include <QtCore/QDate>

namespace OpenMS
{
  class String;

  /**
      @brief DateTime Class.

      This class implements date handling.
      Import and export to/from both string and integers is possible.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI DateTime
  {
public:

    /**
        @brief Default constructor

        Fills the object with an undefined date: 00/00/0000
    */
    DateTime();

    /// Copy constructor
    DateTime(const DateTime& date) = default;

    /// Move constructor
    DateTime(DateTime&&) = default;

    /// Assignment operator
    DateTime& operator=(const DateTime& source);
	  
    /// Move assignment operator
    DateTime& operator=(DateTime&&) & noexcept;

    /// equal operator
    bool operator==(const DateTime& rhs) const;

    /// not-equal operator
    bool operator!=(const DateTime& rhs) const;

    /// less operator
    bool operator<(const DateTime& rhs) const;

    /**
        @brief sets date from a string

        Reads both English, German and iso/ansi date formats: 'MM/dd/yyyy', 'dd.MM.yyyy' or 'yyyy-MM-dd'

        @exception Exception::ParseError
    */
    void setDate(const String& date);

    /**
        @brief sets time from a string

        Reads time format: 'hh:mm:ss'

        @exception Exception::ParseError
    */
    void setTime(const String& date);

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
    void get(UInt& month, UInt& day, UInt& year, UInt& hour, UInt& minute, UInt& second) const;

    /**
        @brief Fills the arguments with the date

        Give the numbers in the following order: month, day and year.
    */
    void getDate(UInt& month, UInt& day, UInt& year) const;

    /**
        @brief Returns the date as string

        The format of the string is yyyy-MM-dd
    */
    String getDate() const;

    /**
        @brief Fills the arguments with the time

        The arguments are all UInts and the order is hour minute second
    */
    void getTime(UInt& hour, UInt& minute, UInt& second) const;

    // add @param s seconds to date time
    DateTime& addSecs(int s);

    /**
        @brief Returns the time as string

        The format of the string is hh:mm:ss
    */
    String getTime() const;

    /// Returns the current date and time
    static DateTime now();

    /// Returns true if the date time is valid
    bool isValid() const;

    /// return true if the date and time is null 
    bool isNull() const;

    /// Sets the undefined date: 00/00/0000 00:00:00
    void clear();
    
    /* @brief Returns a string representation of the DateTime object.
       @param format "yyyy-MM-ddThh:mm:ss" corresponds to ISO 8601 and should be preferred.	   
	*/
	String toString(std::string format = "yyyy-MM-ddThh:mm:ss") const;

    /* @brief Creates a DateTime object from string representation.
       @param format "yyyy-MM-ddThh:mm:ss" corresponds to ISO 8601 and should be preferred.
	*/
    static DateTime fromString(const std::string& date, std::string format = "yyyy-MM-ddThh:mm:ss");

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
    void set(const String& date);

private:
    QDateTime dt_;
  };

} // namespace OPENMS

