// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_CONTACTPERSON_H
#define OPENMS_METADATA_CONTACTPERSON_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS 
{
	/**
		@brief Contact person information
		
		
		
		@ingroup Metadata
	*/
  class OPENMS_DLLAPI ContactPerson
  	: public MetaInfoInterface
  {
    public:
    	/// Constructor
      ContactPerson();
      /// Copy constructor
      ContactPerson(const ContactPerson& source);
      /// Destructor
      ~ContactPerson();
      
      /// Assignment operator
      ContactPerson& operator= (const ContactPerson& source);

      /// Equality operator
      bool operator== (const ContactPerson& rhs) const;
      /// Equality operator
      bool operator!= (const ContactPerson& rhs) const;
			
			/// returns the first name of the person
      const String& getFirstName() const;
      /// sets the first name of the person
      void setFirstName(const String& name);
			
			/// returns the last name of the person
      const String& getLastName() const;
      /// sets the last name of the person
      void setLastName(const String& name);
			
      /// sets the full name of the person (gets split into first and last name internally)
      void setName(const String& name);
			
			/// returns the affiliation
      const String& getInstitution() const;
      /// sets the affiliation
      void setInstitution(const String& institution);
			
			/// returns the email address
      const String& getEmail() const;
      /// sets the email address
      void setEmail(const String& email);

			/// returns the email address
      const String& getURL() const;
      /// sets the email address
      void setURL(const String& email);

			/// returns the address
      const String& getAddress() const;
      /// sets the address
      void setAddress(const String& email);
			
			/// returns miscellaneous info about the contact person
      const String& getContactInfo() const;
      /// sets miscellaneous info about the contact person
      void setContactInfo(const String& contact_info);

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

#endif // OPENMS_METADATA_CONTACTPERSON_H
