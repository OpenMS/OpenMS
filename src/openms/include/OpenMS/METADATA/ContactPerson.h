// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

