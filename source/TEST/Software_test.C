// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

START_TEST(Software, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DateTime time;
time.set("2000-10-09 08:07:40");

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

CHECK(const DateTime& getCompletionTime() const)
  Software tmp;
	String str;
	tmp.getCompletionTime().get(str);
  TEST_EQUAL(str,"0000-00-00 00:00:00");
RESULT

CHECK(void setCompletionTime(const DateTime& completion_time))
  Software tmp;
  tmp.setCompletionTime(time);
  TEST_EQUAL(tmp.getCompletionTime()==time,true);
RESULT

CHECK(void setCompletionTime(const String& completion_time))
  Software tmp;
  tmp.setCompletionTime("2000-10-09 08:07:40");
  TEST_EQUAL(tmp.getCompletionTime()==time,true);
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
  tmp.setCompletionTime(time);
  tmp.setComment("bla");

  Software tmp2(tmp);
  TEST_EQUAL(tmp.getCompletionTime()==time,true);
  TEST_EQUAL(tmp.getVersion(),"0.54");
  TEST_EQUAL(tmp.getName(),"name");
  TEST_EQUAL(tmp.getComment(),"bla");
RESULT


CHECK(Software& operator= (const Software& source))
  Software tmp;
  tmp.setVersion("0.54");
  tmp.setName("name");
  tmp.setCompletionTime(time);
  tmp.setComment("bla");
  
  Software tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp.getCompletionTime()==time,true);
  TEST_EQUAL(tmp2.getVersion(),"0.54");
  TEST_EQUAL(tmp2.getName(),"name");
  TEST_EQUAL(tmp2.getComment(),"bla");

  tmp2 = Software();
  TEST_EQUAL(tmp2.getCompletionTime()==DateTime(),true);
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
  edit.setCompletionTime(time);
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
  edit.setCompletionTime(time);
  TEST_EQUAL(edit!=empty,true);
  
  edit = empty;
  edit.setComment("bla");
  TEST_EQUAL(edit!=empty,true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



