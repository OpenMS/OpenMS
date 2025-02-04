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
#include <OpenMS/METADATA/Acquisition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Acquisition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Acquisition* ptr = nullptr;
Acquisition* nullPointer = nullptr;
START_SECTION(Acquisition())
	ptr = new Acquisition();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~Acquisition())
	delete ptr;
END_SECTION

START_SECTION(const String& getIdentifier() const)
  Acquisition tmp;
  TEST_EQUAL(tmp.getIdentifier(), "");
END_SECTION

START_SECTION(void setIdentifier(const String& identifier))
	Acquisition tmp;
	tmp.setIdentifier("5");
  TEST_EQUAL(tmp.getIdentifier(), "5");
END_SECTION

START_SECTION(Acquisition(const Acquisition& source))
	Acquisition tmp;
	tmp.setIdentifier("5");
	tmp.setMetaValue("label",String("label"));
	Acquisition tmp2(tmp);
	TEST_EQUAL(tmp2.getIdentifier(), "5");
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
END_SECTION

START_SECTION(Acquisition(Acquisition&&) = default)
  Acquisition e, empty;
  e.setIdentifier("Ident");

  Acquisition ef(e);
  Acquisition ef2(e);

  TEST_FALSE(ef == empty)

  // the move target should be equal, while the move source should be empty
  Acquisition ef_mv(std::move(ef));
  TEST_TRUE(ef_mv == ef2)
  TEST_TRUE(ef == empty)
  TEST_EQUAL(ef.getIdentifier().empty(), true)
END_SECTION

START_SECTION(Acquisition& operator= (const Acquisition& source))
	Acquisition tmp,tmp2,tmp3;
	// assignment of a modified object
	tmp2.setIdentifier("5");
	tmp2.setMetaValue("label",String("label"));
	tmp = tmp2;
	TEST_EQUAL(tmp.getIdentifier(), "5");
	TEST_EQUAL((String)(tmp.getMetaValue("label")), String("label"));
	
	// assignment of a default-constructed object
	tmp = tmp3;
	TEST_EQUAL(tmp.getIdentifier(), "");
	TEST_EQUAL(tmp.isMetaEmpty(), true);	
END_SECTION

START_SECTION(bool operator== (const Acquisition& rhs) const)
	Acquisition tmp,tmp2;
	
	TEST_TRUE(tmp == tmp2);
	
	tmp2.setIdentifier("5");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp==tmp2, false);
END_SECTION

START_SECTION(bool operator!= (const Acquisition& rhs) const)
	Acquisition tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setIdentifier("5");
	TEST_FALSE(tmp == tmp2);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_FALSE(tmp == tmp2);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



