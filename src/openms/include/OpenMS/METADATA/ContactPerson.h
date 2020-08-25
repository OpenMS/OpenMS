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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
      @brief Contact person information

      @ingroup Metadata
  */
  class OPENMS_DLLAPI ContactPerson :
    public MetaInfoInterface
  {
public:

    /// Constructor
    ContactPerson() = default;
    /// Copy constructor
    ContactPerson(const ContactPerson &) = default;
    /// Move constructor
    ContactPerson(ContactPerson&&) = default;
    /// Destructor
    ~ContactPerson() = default;

    /// Assignment operator
    ContactPerson & operator=(const ContactPerson &) = default;
    /// Move assignment operator
    ContactPerson& operator=(ContactPerson&&) & = default;

    /// Equality operator
    bool operator==(const ContactPerson & rhs) const;
    /// Equality operator
    bool operator!=(const ContactPerson & rhs) const;

    /// returns the first name of the person
    const String & getFirstName() const;
    /// sets the first name of the person
    void setFirstName(const String & name);

    /// returns the last name of the person
    const String & getLastName() const;
    /// sets the last name of the person
    void setLastName(const String & name);

    /// sets the full name of the person (gets split into first and last name internally)
    void setName(const String & name);

    /// returns the affiliation
    const String & getInstitution() const;
    /// sets the affiliation
    void setInstitution(const String & institution);

    /// returns the email address
    const String & getEmail() const;
    /// sets the email address
    void setEmail(const String & email);

    /// returns the email address
    const String & getURL() const;
    /// sets the email address
    void setURL(const String & email);

    /// returns the address
    const String & getAddress() const;
    /// sets the address
    void setAddress(const String & email);

    /// returns miscellaneous info about the contact person
    const String & getContactInfo() const;
    /// sets miscellaneous info about the contact person
    void setContactInfo(const String & contact_info);

protected:
    String first_name_;
    String last_name_;
    String institution_;
    String email_;
    String contact_info_;
    String url_;
    String address_;
  };
} // namespace OpenMS

