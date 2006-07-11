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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/Software.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Software, "$Id: Software_test.C,v 1.3 2006/05/30 15:46:43 marc_sturm Exp $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Software* ptr = 0;
CHECK(Software())
	ptr = new Software();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~Software())
	delete ptr;
RESULT


CHECK(const String& getComment() const)
  Software tmp;
  TEST_EQUAL(tmp.getComment(),"");
RESULT

CHECK(void setComment(const String& comment))
  Software tmp;
  tmp.setComment("comment");
  TEST_EQUAL(tmp.getComment(),"comment");
RESULT

CHECK(const String& getName() const)
  Software tmp;
  TEST_EQUAL(tmp.getName(),"");
RESULT

CHECK(void setName(const String& name))
  Software tmp;
  tmp.setName("name");
  TEST_EQUAL(tmp.getName(),"name");  
RESULT

CHECK(float getCompletionTime() const)
  Software tmp;
  TEST_EQUAL(tmp.getCompletionTime(),0.0);
RESULT

CHECK(void setCompletionTime(float completion_time))
  Software tmp;
  tmp.setCompletionTime(47.11);
  TEST_REAL_EQUAL(tmp.getCompletionTime(),47.11);
RESULT

CHECK(const String& getVersion() const)
  Software tmp;
  TEST_EQUAL(tmp.getVersion(),"");
RESULT

CHECK(void setVersion(const String& version))
  Software tmp;
  tmp.setVersion("0.54");
  TEST_EQUAL(tmp.getVersion(),"0.54");
RESULT

CHECK(Software(const Software& source))
  Software tmp;
  tmp.setVersion("0.54");
  tmp.setName("name");
  tmp.setCompletionTime(1.1);
  tmp.setComment("bla");
  
  Software tmp2(tmp);
  TEST_REAL_EQUAL(tmp.getCompletionTime(),1.1);
  TEST_EQUAL(tmp.getVersion(),"0.54");
  TEST_EQUAL(tmp.getName(),"name");
  TEST_EQUAL(tmp.getComment(),"bla");
RESULT


CHECK(Software& operator= (const Software& source))
  Software tmp;
  tmp.setVersion("0.54");
  tmp.setName("name");
  tmp.setCompletionTime(1.1);
  tmp.setComment("bla");
  
  Software tmp2;
  tmp2 = tmp;
  TEST_REAL_EQUAL(tmp.getCompletionTime(),1.1);
  TEST_EQUAL(tmp2.getVersion(),"0.54");
  TEST_EQUAL(tmp2.getName(),"name");
  TEST_EQUAL(tmp2.getComment(),"bla");

  tmp2 = Software();
  TEST_EQUAL(tmp2.getCompletionTime(),0.0);
  TEST_EQUAL(tmp2.getVersion(),"");
  TEST_EQUAL(tmp2.getName(),"");
  TEST_EQUAL(tmp2.getComment(),"");
RESULT


CHECK(bool operator== (const Software& rhs) const)
  Software edit, empty;
  
  TEST_EQUAL(edit==empty,true);
  
  edit = empty;
  edit.setVersion("0.54");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setName("name");
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setCompletionTime(1.1);
  TEST_EQUAL(edit==empty,false);
  
  edit = empty;
  edit.setComment("bla");
  TEST_EQUAL(edit==empty,false);
RESULT

CHECK(bool operator!= (const Software& rhs) const)
  Software edit, empty;
  
  TEST_EQUAL(edit!=empty,false);
  
  edit = empty;
  edit.setVersion("0.54");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setName("name");
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setCompletionTime(1.1);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setComment("bla");
  TEST_EQUAL(edit!=empty,true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



