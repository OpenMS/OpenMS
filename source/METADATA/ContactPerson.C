// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ContactPerson.h>

using namespace std;

namespace OpenMS
{

  ContactPerson::ContactPerson() :
    MetaInfoInterface(),
    first_name_(),
    last_name_(),
    institution_(),
    email_(),
    contact_info_(),
    url_(),
    address_()
  {

  }

  ContactPerson::ContactPerson(const ContactPerson & source) :
    MetaInfoInterface(source),
    first_name_(source.first_name_),
    last_name_(source.last_name_),
    institution_(source.institution_),
    email_(source.email_),
    contact_info_(source.contact_info_),
    url_(source.url_),
    address_(source.address_)
  {

  }

  ContactPerson::~ContactPerson()
  {

  }

  ContactPerson & ContactPerson::operator=(const ContactPerson & source)
  {
    if (&source == this)
      return *this;

    first_name_ = source.first_name_;
    last_name_ = source.last_name_;
    institution_ = source.institution_;
    email_ = source.email_;
    contact_info_ = source.contact_info_;
    url_ = source.url_;
    address_ = source.address_;
    MetaInfoInterface::operator=(source);

    return *this;
  }

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
