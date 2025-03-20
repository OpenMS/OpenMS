// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ContactPerson.h>

using namespace std;

namespace OpenMS
{

  bool ContactPerson::operator==(const ContactPerson & rhs) const
  {
    return first_name_ == rhs.first_name_ &&
           last_name_ == rhs.last_name_ &&
           institution_ == rhs.institution_ &&
           email_ == rhs.email_ &&
           contact_info_ == rhs.contact_info_ &&
           url_ == rhs.url_ &&
           address_ == rhs.address_ &&
           MetaInfoInterface::operator==(rhs);
  }

  bool ContactPerson::operator!=(const ContactPerson & rhs) const
  {
    return !(operator==(rhs));
  }

  const String & ContactPerson::getFirstName() const
  {
    return first_name_;
  }

  void ContactPerson::setFirstName(const String & name)
  {
    first_name_ = name;
  }

  const String & ContactPerson::getLastName() const
  {
    return last_name_;
  }

  void ContactPerson::setLastName(const String & name)
  {
    last_name_ = name;
  }

  void ContactPerson::setName(const String & name)
  {
    std::vector<String> tmp;
    if (name.split(',', tmp))
    {
      first_name_ = tmp[1].trim();
      last_name_ = tmp[0].trim();
    }
    else
    {
      if (name.split(' ', tmp))
      {
        first_name_ = tmp[0];
        last_name_ = tmp[1];
      }
      else
      {
        last_name_ = name;
      }
    }
  }

  const String & ContactPerson::getEmail() const
  {
    return email_;
  }

  void ContactPerson::setEmail(const String & email)
  {
    email_ = email;
  }

  const String & ContactPerson::getInstitution() const
  {
    return institution_;
  }

  void ContactPerson::setInstitution(const String & institution)
  {
    institution_ = institution;
  }

  const String & ContactPerson::getContactInfo() const
  {
    return contact_info_;
  }

  void ContactPerson::setContactInfo(const String & contact_info)
  {
    contact_info_ = contact_info;
  }

  const String & ContactPerson::getURL() const
  {
    return url_;
  }

  void ContactPerson::setURL(const String & url)
  {
    url_ = url;
  }

  const String & ContactPerson::getAddress() const
  {
    return address_;
  }

  void ContactPerson::setAddress(const String & address)
  {
    address_ = address;
  }

}

