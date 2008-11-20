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

#include <OpenMS/FORMAT/ControlledVocabulary.h>

///////////////////////////

START_TEST(ControlledVocabulary, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

ControlledVocabulary* ptr = 0;
START_SECTION((ControlledVocabulary()))
	ptr = new ControlledVocabulary();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~ControlledVocabulary()))
	delete ptr;
END_SECTION

START_SECTION(const String& name() const)
	TEST_EQUAL(ControlledVocabulary().name(),"")
END_SECTION

ControlledVocabulary cv;
START_SECTION(void loadFromOBO(const String &name, const String &filename))
	cv.loadFromOBO("bla","data/ControlledVocabulary.obo");
	TEST_EQUAL(cv.name(),"bla")
END_SECTION

START_SECTION(bool exists(const String& id) const)
	TEST_EQUAL(cv.exists("OpenMS:1"),true)
	TEST_EQUAL(cv.exists("OpenMS:2"),true)
	TEST_EQUAL(cv.exists("OpenMS:3"),true)
	TEST_EQUAL(cv.exists("OpenMS:4"),true)
	TEST_EQUAL(cv.exists("OpenMS:5"),true)
	TEST_EQUAL(cv.exists("OpenMS:6"),true)
	TEST_EQUAL(cv.exists("OpenMS:7"),false)
END_SECTION

START_SECTION(const CVTerm& getTerm(const String& id) const)
	const ControlledVocabulary::CVTerm* term=0;
	//Auto
	term = &(cv.getTerm("OpenMS:1"));
	TEST_EQUAL(term->id,"OpenMS:1")
	TEST_EQUAL(term->name,"Auto")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),0)
	TEST_EQUAL(term->unparsed.size(),0)
	//Ford
	term = &(cv.getTerm("OpenMS:2"));
	TEST_EQUAL(term->id,"OpenMS:2")
	TEST_EQUAL(term->name,"Ford")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:1")
	TEST_EQUAL(term->unparsed.size(),0)
	//Mercedes
	term = &(cv.getTerm("OpenMS:3"));
	TEST_EQUAL(term->id,"OpenMS:3")
	TEST_EQUAL(term->name,"Mercedes")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:1")
	TEST_EQUAL(term->unparsed.size(),0)
	//A-Klasse
	term = &(cv.getTerm("OpenMS:4"));
	TEST_EQUAL(term->id,"OpenMS:4")
	TEST_EQUAL(term->name,"A-Klasse")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:3")
	TEST_EQUAL(term->unparsed.size(),3)
	TEST_EQUAL(term->unparsed[0],"xref: unparsed line 1")
	TEST_EQUAL(term->unparsed[1],"xref: unparsed line 2")
	TEST_EQUAL(term->unparsed[2],"xref: unparsed line 3")
	//Mustang
	term = &(cv.getTerm("OpenMS:5"));
	TEST_EQUAL(term->id,"OpenMS:5")
	TEST_EQUAL(term->name,"Mustang")
	TEST_EQUAL(term->obsolete,false)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:2")
	TEST_EQUAL(term->unparsed.size(),0)
	//Ka
	term = &(cv.getTerm("OpenMS:6"));
	TEST_EQUAL(term->id,"OpenMS:6")
	TEST_EQUAL(term->name,"Ka")
	TEST_EQUAL(term->obsolete,true)
	TEST_EQUAL(term->parents.size(),1)
	TEST_EQUAL(*term->parents.begin(),"OpenMS:2")
	TEST_EQUAL(term->unparsed.size(),0)
	
	TEST_EXCEPTION(Exception::InvalidValue , cv.getTerm("OpenMS:7"))
END_SECTION

START_SECTION(bool isChildOf(const String& child, const String& parent) const)
	TEST_EQUAL(cv.isChildOf("OpenMS:6","OpenMS:2"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:5","OpenMS:2"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:2","OpenMS:1"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:3","OpenMS:1"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:4","OpenMS:3"),true)
	TEST_EQUAL(cv.isChildOf("OpenMS:1","OpenMS:6"),false)
	TEST_EQUAL(cv.isChildOf("OpenMS:4","OpenMS:6"),false)
	TEST_EQUAL(cv.isChildOf("OpenMS:2","OpenMS:6"),false)
	TEST_EQUAL(cv.isChildOf("OpenMS:2","OpenMS:3"),false)
	TEST_EXCEPTION(Exception::InvalidValue, cv.isChildOf("OpenMS:7","OpenMS:3"))
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
