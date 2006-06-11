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
// $Id: ContactPerson.h,v 1.1 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
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
  class ContactPerson: public MetaInfoInterface
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
			
			/// returns the name of the person
      const String& getName() const;
      /// sets the name of the person
      void setName(const String& name);
			
			/// returns the affiliation
      const String& getInstitution() const;
      /// sets the affiliation
      void setInstitution(const String& institution);
			
			/// returns the email address
      const String& getEmail() const;
      /// sets the email address
      void setEmail(const String& email);
			
			/// returns miscellaneous info about the contact person
      const String& getContactInfo() const;
      /// sets miscellaneous info about the contact person
      void setContactInfo(const String& contact_info);

    protected:
			String name_;
			String institution_;
			String email_;
			String contact_info_;
  };
} // namespace OpenMS

#endif // OPENMS_METADATA_CONTACTPERSON_H
