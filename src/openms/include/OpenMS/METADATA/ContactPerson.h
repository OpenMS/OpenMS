// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
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
    const std::string & getFirstName() const;
    /// sets the first name of the person
    void setFirstName(const std::string & name);

    /// returns the last name of the person
    const std::string & getLastName() const;
    /// sets the last name of the person
    void setLastName(const std::string & name);

    /// sets the full name of the person (gets split into first and last name internally)
    void setName(const std::string & name);

    /// returns the affiliation
    const std::string & getInstitution() const;
    /// sets the affiliation
    void setInstitution(const std::string & institution);

    /// returns the email address
    const std::string & getEmail() const;
    /// sets the email address
    void setEmail(const std::string & email);

    /// returns the email address
    const std::string & getURL() const;
    /// sets the email address
    void setURL(const std::string & email);

    /// returns the address
    const std::string & getAddress() const;
    /// sets the address
    void setAddress(const std::string & email);

    /// returns miscellaneous info about the contact person
    const std::string & getContactInfo() const;
    /// sets miscellaneous info about the contact person
    void setContactInfo(const std::string & contact_info);

protected:
    std::string first_name_;
    std::string last_name_;
    std::string institution_;
    std::string email_;
    std::string contact_info_;
    std::string url_;
    std::string address_;
  };
} // namespace OpenMS

