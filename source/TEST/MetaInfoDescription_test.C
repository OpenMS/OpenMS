// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/METADATA/MetaInfoDescription.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MetaInfoDescription, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MetaInfoDescription* ptr = 0;
START_SECTION((MetaInfoDescription()))
	ptr = new MetaInfoDescription();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~MetaInfoDescription()))
	delete ptr;
END_SECTION

START_SECTION((const String& getName() const))
  MetaInfoDescription tmp;
  TEST_EQUAL(tmp.getName(),"");
END_SECTION

START_SECTION((void setName(const String& name)))
  MetaInfoDescription tmp;
  tmp.setName("name");
  TEST_EQUAL(tmp.getName(),"name");
END_SECTION

START_SECTION((const std::vector<DataProcessing>& getDataProcessing() const))
  MetaInfoDescription tmp;
  TEST_EQUAL(tmp.getDataProcessing().size(),0);
END_SECTION

START_SECTION((void setDataProcessing(const std::vector< DataProcessing > &data_processing)))
  MetaInfoDescription tmp;
  std::vector<DataProcessing> dummy;
  dummy.resize(1);
  tmp.setDataProcessing(dummy);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((std::vector<DataProcessing>& getDataProcessing()))
  MetaInfoDescription tmp;
  tmp.getDataProcessing().resize(1);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((MetaInfoDescription(const MetaInfoDescription& source)))
  MetaInfoDescription tmp;
  tmp.setName("bla2");
  tmp.getDataProcessing().resize(1);
  tmp.setMetaValue("label",String("label"));
  
  MetaInfoDescription tmp2(tmp);
  TEST_EQUAL(tmp2.getName(),"bla2");
	TEST_EQUAL(tmp.getDataProcessing().size(),1);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
END_SECTION

START_SECTION((MetaInfoDescription& operator= (const MetaInfoDescription& source)))
  MetaInfoDescription tmp;
  tmp.setName("bla2");
  tmp.getDataProcessing().resize(1);
  tmp.setMetaValue("label",String("label"));
  
  MetaInfoDescription tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getName(),"bla2");
	TEST_EQUAL(tmp.getDataProcessing().size(),1);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	
  tmp2 = MetaInfoDescription();
  TEST_EQUAL(tmp2.getName(),"");
	TEST_EQUAL(tmp2.getDataProcessing().size(),0);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
END_SECTION

START_SECTION((bool operator== (const MetaInfoDescription& rhs) const))
  MetaInfoDescription edit, empty;
  
  TEST_EQUAL(edit==empty, true);
  
  edit = empty;
  edit.setName("bla2");
	TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setMetaValue("label",String("label"));
  TEST_EQUAL(edit==empty, false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



