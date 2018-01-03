// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



