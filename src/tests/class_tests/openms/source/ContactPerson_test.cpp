// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/ContactPerson.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ContactPerson, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ContactPerson* ptr = nullptr;
ContactPerson* nullPointer = nullptr;
START_SECTION(ContactPerson())
	ptr = new ContactPerson();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~ContactPerson())
	delete ptr;
END_SECTION

START_SECTION(const String& getContactInfo() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getContactInfo(),"");
END_SECTION

START_SECTION(void setContactInfo(const String& contact_info))
  ContactPerson tmp;
	tmp.setContactInfo("bla");
  TEST_EQUAL(tmp.getContactInfo(),"bla");
END_SECTION

START_SECTION(const String& getEmail() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getEmail(),"");
END_SECTION

START_SECTION(void setEmail(const String& email))
  ContactPerson tmp;
	tmp.setEmail("bla");
  TEST_EQUAL(tmp.getEmail(),"bla");
END_SECTION

START_SECTION(const String& getURL() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getURL(),"");
END_SECTION

START_SECTION(void setURL(const String& email))
  ContactPerson tmp;
	tmp.setURL("bla");
  TEST_EQUAL(tmp.getURL(),"bla");
END_SECTION

START_SECTION(const String& getAddress() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getAddress(),"");
END_SECTION

START_SECTION(void setAddress(const String& email))
  ContactPerson tmp;
	tmp.setAddress("bla");
  TEST_EQUAL(tmp.getAddress(),"bla");
END_SECTION


START_SECTION(const String& getInstitution() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getInstitution(),"");
END_SECTION

START_SECTION(void setInstitution(const String& institution))
  ContactPerson tmp;
	tmp.setInstitution("Uni Tuebingen");
  TEST_EQUAL(tmp.getInstitution(),"Uni Tuebingen");
END_SECTION

START_SECTION(const String& getFirstName() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getFirstName(),"");
END_SECTION

START_SECTION(void setFirstName(const String& name))
  ContactPerson tmp;
	tmp.setFirstName("Meike");
  TEST_EQUAL(tmp.getFirstName(),"Meike");
END_SECTION

START_SECTION(const String& getLastName() const)
  ContactPerson tmp;
  TEST_EQUAL(tmp.getLastName(),"");
END_SECTION

START_SECTION(void setLastName(const String& name))
  ContactPerson tmp;
	tmp.setLastName("Meier");
  TEST_EQUAL(tmp.getLastName(),"Meier");
END_SECTION

START_SECTION(void setName(const String& name))
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
END_SECTION

START_SECTION(ContactPerson(const ContactPerson& source))
	ContactPerson tmp;
	tmp.setEmail("ich@du.de");
	tmp.setFirstName("Meike");
	tmp.setLastName("Meier");
	tmp.setInstitution("Uni Tuebingen");
	tmp.setContactInfo("doo");
	tmp.setURL("url");
	tmp.setAddress("street");
	tmp.setMetaValue("label",String("label"));
	
	ContactPerson tmp2(tmp);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
  TEST_EQUAL(tmp2.getFirstName(),"Meike");
  TEST_EQUAL(tmp2.getLastName(),"Meier");
  TEST_EQUAL(tmp2.getEmail(),"ich@du.de");
  TEST_EQUAL(tmp2.getInstitution(),"Uni Tuebingen");
  TEST_EQUAL(tmp2.getContactInfo(),"doo");
  TEST_EQUAL(tmp2.getURL(),"url");
  TEST_EQUAL(tmp2.getAddress(),"street");
END_SECTION

START_SECTION(ContactPerson& operator= (const ContactPerson& source))
 	ContactPerson tmp;
	tmp.setEmail("ich@du.de");
	tmp.setFirstName("Meike");
	tmp.setLastName("Meier");
	tmp.setInstitution("Uni Tuebingen");
	tmp.setContactInfo("doo");
	tmp.setURL("url");
	tmp.setAddress("street");
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
  TEST_EQUAL(tmp2.getURL(),"url");
  TEST_EQUAL(tmp2.getAddress(),"street");

	//assignment of empty object
	tmp2 = ContactPerson();
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
  TEST_EQUAL(tmp2.getFirstName(),"");
  TEST_EQUAL(tmp2.getLastName(),"");
  TEST_EQUAL(tmp2.getEmail(),"");
  TEST_EQUAL(tmp2.getInstitution(),"");
  TEST_EQUAL(tmp2.getContactInfo(),"");
  TEST_EQUAL(tmp2.getURL(),"");
  TEST_EQUAL(tmp2.getAddress(),"");
END_SECTION

START_SECTION(bool operator!= (const ContactPerson& rhs) const)
	ContactPerson tmp,tmp2;
	
	TEST_TRUE(tmp == tmp2);
	
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
	tmp.setURL("url");
  TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
	tmp.setAddress("street");
  TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
  TEST_EQUAL(tmp==tmp2, false);
END_SECTION

START_SECTION(bool operator== (const ContactPerson& rhs) const)
	ContactPerson tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp.setEmail("ich@du.de");
	TEST_FALSE(tmp == tmp2);
	
	tmp2 = tmp;
	tmp.setFirstName("Meike");
	TEST_FALSE(tmp == tmp2);
	
	tmp2 = tmp;
	tmp.setLastName("Meier");
	TEST_FALSE(tmp == tmp2);
	
	tmp2 = tmp;
	tmp.setInstitution("Uni Tuebingen");
  TEST_FALSE(tmp == tmp2);
	
	tmp2 = tmp;
	tmp.setContactInfo("doo");
  TEST_FALSE(tmp == tmp2);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
  TEST_FALSE(tmp == tmp2);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



