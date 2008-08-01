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
#include <OpenMS/METADATA/MetaInfoDescription.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MetaInfoDescription, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MetaInfoDescription* ptr = 0;
CHECK(MetaInfoDescription())
	ptr = new MetaInfoDescription();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~MetaInfoDescription())
	delete ptr;
RESULT

CHECK(const SourceFile& getSourceFile() const)
  MetaInfoDescription tmp;
  TEST_EQUAL(tmp.getSourceFile()==SourceFile(),true);
RESULT

CHECK(void setSourceFile(const SourceFile& source_file))
  MetaInfoDescription tmp;
	SourceFile sf;
	sf.setFileType("rm");
	tmp.setSourceFile(sf);
  TEST_EQUAL(tmp.getSourceFile().getFileType(),"rm");
RESULT

CHECK(SourceFile& getSourceFile())
  MetaInfoDescription tmp;
	tmp.getSourceFile().setFileType("rm");
  TEST_EQUAL(tmp.getSourceFile().getFileType(),"rm");
RESULT

CHECK(const String& getComment() const)
  MetaInfoDescription tmp;
  TEST_EQUAL(tmp.getComment(),"");
RESULT

CHECK(void setComment(const String& comment))
  MetaInfoDescription tmp;
  tmp.setComment("comment");
  TEST_EQUAL(tmp.getComment(),"comment");
RESULT

CHECK(const String& getName() const)
  MetaInfoDescription tmp;
  TEST_EQUAL(tmp.getName(),"");
RESULT

CHECK(void setName(const String& name))
  MetaInfoDescription tmp;
  tmp.setName("name");
  TEST_EQUAL(tmp.getName(),"name");
RESULT

CHECK(MetaInfoDescription(const MetaInfoDescription& source))
  MetaInfoDescription tmp;
  tmp.getSourceFile().setFileType("wma");
  tmp.setComment("bla");
  tmp.setName("bla2");
  tmp.setMetaValue("label",String("label"));
  
  MetaInfoDescription tmp2(tmp);
  TEST_EQUAL(tmp2.getSourceFile().getFileType(),"wma");
  TEST_EQUAL(tmp2.getComment(),"bla");
  TEST_EQUAL(tmp2.getName(),"bla2");
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
RESULT

CHECK(MetaInfoDescription& operator= (const MetaInfoDescription& source))
  MetaInfoDescription tmp;
  tmp.getSourceFile().setFileType("wma");
  tmp.setComment("bla");
  tmp.setName("bla2");
  tmp.setMetaValue("label",String("label"));
  
  MetaInfoDescription tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getSourceFile().getFileType(),"wma");
  TEST_EQUAL(tmp2.getComment(),"bla");
  TEST_EQUAL(tmp2.getName(),"bla2");
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");

  tmp2 = MetaInfoDescription();
  TEST_EQUAL(tmp2.getSourceFile().getFileType(),"");
  TEST_EQUAL(tmp2.getComment(),"");
  TEST_EQUAL(tmp2.getName(),"");
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
RESULT

CHECK(bool operator== (const MetaInfoDescription& rhs) const)
  MetaInfoDescription edit, empty;
  
  TEST_EQUAL(edit==empty, true);
  
  edit.getSourceFile().setFileType("wma");
  TEST_EQUAL(edit==empty, false);
  
  edit = empty;
  edit.setComment("bla");
	TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setName("bla2");
	TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setMetaValue("label",String("label"));
  TEST_EQUAL(edit==empty, false);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



