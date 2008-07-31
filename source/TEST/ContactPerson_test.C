// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/ContactPerson.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ContactPerson, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ContactPerson* ptr = 0;
CHECK(ContactPerson())
	ptr = new ContactPerson();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ContactPerson())
	delete ptr;
RESULT

CHECK(const String& getContactInfo() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getContactInfo(),"");
RESULT

CHECK(void setContactInfo(const String& contact_info))
  ContactPerson tmp;
	tmp.setContactInfo("bla");
  TEST_EQUAL(tmp.getContactInfo(),"bla");
RESULT

CHECK(const String& getEmail() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getEmail(),"");
RESULT

CHECK(void setEmail(const String& email))
  ContactPerson tmp;
	tmp.setEmail("bla");
  TEST_EQUAL(tmp.getEmail(),"bla");
RESULT

CHECK(const String& getInstitution() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getInstitution(),"");
RESULT

CHECK(void setInstitution(const String& institution))
  ContactPerson tmp;
	tmp.setInstitution("Uni Tuebingen");
  TEST_EQUAL(tmp.getInstitution(),"Uni Tuebingen");
RESULT

CHECK(const String& getFirstName() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getFirstName(),"");
RESULT

CHECK(void setFirstName(const String& name))
  ContactPerson tmp;
	tmp.setFirstName("Meike");
  TEST_EQUAL(tmp.getFirstName(),"Meike");
RESULT

CHECK(const String& getLastName() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getLastName(),"");
RESULT

CHECK(void setLastName(const String& name))
  ContactPerson tmp;
	tmp.setLastName("Meier");
  TEST_EQUAL(tmp.getLastName(),"Meier");
RESULT

CHECK(void setName(const String& name))
  ContactPerson tmp;
	tmp.setName("Diddl Maus");
  TEST_EQUAL(tmp.getFirstName(),"Diddl");
  TEST_EQUAL(tmp.getLastName(),"Maus");
	tmp.setName("Normal, Otto");
  TEST_EQUAL(tmp.getFirstName(),"Otto");
  TEST_EQUAL(tmp.getLastName(),"Normal");
	tmp.setName("Meiser, Hans F.");
  TEST_EQUAL(tmp.getFirstName(),"Hans F.");
  TEST_EQUAL(tmp.getLastName(),"Meiser");
RESULT

CHECK(ContactPerson(const ContactPerson& source))
	ContactPerson tmp;
	tmp.setEmail("ich@du.de");
	tmp.setFirstName("Meike");
	tmp.setLastName("Meier");
	tmp.setInstitution("Uni Tuebingen");
	tmp.setContactInfo("doo");
	tmp.setMetaValue("label",String("label"));
	
	ContactPerson tmp2(tmp);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_EQUAL(tmp2.getFirstName(),"Meike");
  TEST_EQUAL(tmp2.getLastName(),"Meier");
  TEST_EQUAL(tmp2.getEmail(),"ich@du.de");
  TEST_EQUAL(tmp2.getInstitution(),"Uni Tuebingen");
  TEST_EQUAL(tmp2.getContactInfo(),"doo");
RESULT

CHECK(ContactPerson& operator= (const ContactPerson& source))
 	ContactPerson tmp;
	tmp.setEmail("ich@du.de");
	tmp.setFirstName("Meike");
	tmp.setLastName("Meier");
	tmp.setInstitution("Uni Tuebingen");
	tmp.setContactInfo("doo");
	tmp.setMetaValue("label",String("label"));
	
	//normal assignment
	ContactPerson tmp2;
	tmp2 = tmp;
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_EQUAL(tmp2.getFirstName(),"Meike");
  TEST_EQUAL(tmp2.getLastName(),"Meier");
  TEST_EQUAL(tmp2.getEmail(),"ich@du.de");
  TEST_EQUAL(tmp2.getInstitution(),"Uni Tuebingen");
  TEST_EQUAL(tmp2.getContactInfo(),"doo");

	//assignment of empty object
	tmp2 = ContactPerson();
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
  TEST_EQUAL(tmp2.getFirstName(),"");
  TEST_EQUAL(tmp2.getLastName(),"");
  TEST_EQUAL(tmp2.getEmail(),"");
  TEST_EQUAL(tmp2.getInstitution(),"");
  TEST_EQUAL(tmp2.getContactInfo(),"");
RESULT

CHECK(bool operator!= (const ContactPerson& rhs) const)
	ContactPerson tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp.setEmail("ich@du.de");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setFirstName("Meike");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setLastName("Meier");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setInstitution("Uni Tuebingen");
  TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setContactInfo("doo");
  TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
  TEST_EQUAL(tmp==tmp2, false);
RESULT

CHECK(bool operator== (const ContactPerson& rhs) const)
	ContactPerson tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp.setEmail("ich@du.de");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setFirstName("Meike");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setLastName("Meier");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setInstitution("Uni Tuebingen");
  TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setContactInfo("doo");
  TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
  TEST_EQUAL(tmp!=tmp2, true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



