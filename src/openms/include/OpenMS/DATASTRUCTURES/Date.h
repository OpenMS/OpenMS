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

#ifndef OPENMS_DATASTRUCTURES_DATE_H
#define OPENMS_DATASTRUCTURES_DATE_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <QtCore/QDate>

#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  /**
      @brief Date Class.

      This class implements date handling.
      Import and export to/from both string and integers is possible.

      @ingroup Datastructures
  */
  class OPENMS_DLLAPI Date :
    public QDate
  {
public:

    /**
        @brief Default constructor

        Fills the object with an undefined date: 00/00/0000
    */
    Date();
    /// Copy constructor
    Date(const Date & date);
    /// Copy constructor from Qt base class
    Date(const QDate & date);

    /// Assignment operator
    Date & operator=(const Date & source);

    /**
        @brief sets data from a string

        The following date formats are supported:
        - mm/dd/yyyy
        - dd.mm.yyyy
        - yyyy-mm-dd

        @exception Exception::ParseError is thrown if the date is given in the wrong format
    */
    void set(const String & date);

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
    void get(UInt & month, UInt & day, UInt & year) const;

    ///Sets the undefined date: 00/00/0000
    void clear();

protected:
  };
} // namespace OPENMS

#endif // OPENMS_DATASTRUCTURES_DATE_H
