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

  const std::string & ContactPerson::getFirstName() const
  {
    return first_name_;
  }

  void ContactPerson::setFirstName(const std::string & name)
  {
    first_name_ = name;
  }

  const std::string & ContactPerson::getLastName() const
  {
    return last_name_;
  }

  void ContactPerson::setLastName(const std::string & name)
  {
    last_name_ = name;
  }

  void ContactPerson::setName(const std::string & name)
  {
    std::vector<std::string> tmp;
    if (StringUtils::split(name, ', ', tmp))
    {
      first_name_ = StringUtils::trim(tmp[1]);
      last_name_ = StringUtils::trim(tmp[0]);
    }
    else
    {
      if (StringUtils::split(name, ' ', tmp))
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

  const std::string & ContactPerson::getEmail() const
  {
    return email_;
  }

  void ContactPerson::setEmail(const std::string & email)
  {
    email_ = email;
  }

  const std::string & ContactPerson::getInstitution() const
  {
    return institution_;
  }

  void ContactPerson::setInstitution(const std::string & institution)
  {
    institution_ = institution;
  }

  const std::string & ContactPerson::getContactInfo() const
  {
    return contact_info_;
  }

  void ContactPerson::setContactInfo(const std::string & contact_info)
  {
    contact_info_ = contact_info;
  }

  const std::string & ContactPerson::getURL() const
  {
    return url_;
  }

  void ContactPerson::setURL(const std::string & url)
  {
    url_ = url;
  }

  const std::string & ContactPerson::getAddress() const
  {
    return address_;
  }

  void ContactPerson::setAddress(const std::string & address)
  {
    address_ = address;
  }

}

