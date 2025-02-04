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
#include <OpenMS/METADATA/Software.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Software, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Software* ptr = nullptr;
Software* nullPointer = nullptr;
START_SECTION(Software())
	ptr = new Software();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~Software())
	delete ptr;
END_SECTION

START_SECTION(const String& getName() const)
  Software tmp;
  TEST_EQUAL(tmp.getName(),"");
END_SECTION

START_SECTION(void setName(const String& name))
  Software tmp;
  tmp.setName("name");
  TEST_EQUAL(tmp.getName(),"name");  
END_SECTION

START_SECTION(const String& getVersion() const)
  Software tmp;
  TEST_EQUAL(tmp.getVersion(),"");
END_SECTION

START_SECTION(void setVersion(const String& version))
  Software tmp;
  tmp.setVersion("0.54");
  TEST_EQUAL(tmp.getVersion(),"0.54");
END_SECTION

START_SECTION(Software(const Software& source))
  Software tmp;
  tmp.setVersion("0.54");
  tmp.setName("name");
	
	Software tmp2(tmp);
  TEST_EQUAL(tmp2.getVersion(),"0.54");
  TEST_EQUAL(tmp2.getName(),"name");
END_SECTION


START_SECTION(Software& operator= (const Software& source))
  Software tmp;
  tmp.setVersion("0.54");
  tmp.setName("name");
  
  Software tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getVersion(),"0.54");
  TEST_EQUAL(tmp2.getName(),"name");

  tmp2 = Software();
  TEST_EQUAL(tmp2.getVersion(),"");
  TEST_EQUAL(tmp2.getName(),"");
END_SECTION


START_SECTION(bool operator== (const Software& rhs) const)
  Software edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit = empty;
  edit.setVersion("0.54");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setName("name");
  TEST_EQUAL(edit==empty,false);

END_SECTION

START_SECTION(bool operator!= (const Software& rhs) const)
  Software edit, empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit = empty;
  edit.setVersion("0.54");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setName("name");
  TEST_EQUAL(edit!=empty,true);

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



