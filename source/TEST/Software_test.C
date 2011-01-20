// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/Software.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Software, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Software* ptr = 0;
START_SECTION(Software())
	ptr = new Software();
	TEST_NOT_EQUAL(ptr, 0)
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



