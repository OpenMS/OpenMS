// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ContactPerson.h>

using namespace std;

namespace OpenMS
{

	ContactPerson::ContactPerson():
		MetaInfoInterface()
	{
	  
	}
	
	ContactPerson::ContactPerson(const ContactPerson& source):
		MetaInfoInterface(source),
	  name_(source.name_),
	  institution_(source.institution_),
	  email_(source.email_),
	  contact_info_(source.contact_info_)
	{
	  
	}
	
	ContactPerson::~ContactPerson()
	{
	  
	}
	
	ContactPerson& ContactPerson::operator = (const ContactPerson& source)
	{
	  if (&source == this) return *this;
	  
    name_ = source.name_;
    institution_ = source.institution_;
    email_ = source.email_;
    contact_info_ = source.contact_info_;
    MetaInfoInterface::operator=(source);
	  
	  return *this;
	}
	
  bool ContactPerson::operator== (const ContactPerson& rhs) const
  {
  	return 
	    name_ == rhs.name_ &&
	    institution_ == rhs.institution_ &&
	    email_ == rhs.email_ &&
	    contact_info_ == rhs.contact_info_ &&
  		MetaInfoInterface::operator==(rhs)
  		;
  }
  
  bool ContactPerson::operator!= (const ContactPerson& rhs) const
  {
  	return !(operator==(rhs));
 	}
	
	const String& ContactPerson::getName() const 
	{
	  return name_; 
	}
	
	void ContactPerson::setName(const String& name)
	{
	  name_ = name; 
	}
	
	const String& ContactPerson::getEmail() const 
	{
	  return email_; 
	}
	
	void ContactPerson::setEmail(const String& email)
	{
	  email_ = email; 
	}
	
	const String& ContactPerson::getInstitution() const 
	{
	  return institution_; 
	}
	
	void ContactPerson::setInstitution(const String& institution)
	{
	  institution_ = institution; 
	}
	
	const String& ContactPerson::getContactInfo() const 
	{
	  return contact_info_; 
	}
	
	void ContactPerson::setContactInfo(const String& contact_info)
	{
	  contact_info_ = contact_info; 
	}

}



